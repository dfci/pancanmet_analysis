# Pathway analysis for PanCan project

rm(list=ls())
setwd('/Users/ereznik/Documents/pancanmet/analysis/pancan/')

library(ggplot2)
library(pheatmap)
library(reshape)
library(gridExtra)

# Set parameters here
minpathsize = 4

# Write a function which grabs the legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

load('/Users/ereznik/Documents/pancanmet/data/pathway_v2/allpathwaydata.Rda')

# List of studies below must be identical to that from pathway script
studies = c('BLCA','BRCA','BRCATang','COAD','KICH','KIRC','OV','PAAD','PAADHussain1','PAADHussain2','PRAD','PRADLODA','STAD')
#studies = c('BRCA','KIRC')
# Remove pathways with average size less than X (use X = 5)
p2rm = which( rowMeans(sizemet) < minpathsize )
s2rm = list()
s2rm = c('BRCA','BRCATang')

# Plot the result
ylab = 'Gene Differential Abundance Score'
xlab = 'Metabolite Differential Abundance Score'
metmelt = melt(pathmety[-p2rm,]); colnames(metmelt) = c('Pathway','Cancer')
metmeltx = melt(pathmetx[-p2rm,]); colnames(metmeltx) = c('Pathway','Cancer')
genemelt = melt(pathgeney[-p2rm,]); colnames(genemelt) = c('Pathway','Cancer')
sizemelt = melt(sizemet[-p2rm,]); colnames(sizemelt) = c('Pathway','Cancer')
mergedata = merge(metmelt,genemelt,by=c('Pathway','Cancer')); colnames(mergedata) = c('Pathway','Cancer','Met','Gene')
mergedata = merge(mergedata,metmeltx,by=c('Pathway','Cancer')); colnames(mergedata) = c('Pathway','Cancer','Met','Gene','Offset')
mergedata = merge(mergedata,sizemelt,by=c('Pathway','Cancer')); colnames(mergedata) = c('Pathway','Cancer','Met','Gene','Offset','MetSize')

###############################################################################################
# Make an "averaged" figure over all studies indicating differential abundance score
pathaverage = data.frame( matrix(0,dim(pathmetx)[1],3))
rownames(pathaverage) = rownames(pathmetx)
colnames(pathaverage) = c('Prop','DA','NumStudies')

# Make some temporary data that removes the BRCA studies
if (length(s2rm) > 0){
  s2rm_idx = which(colnames(pathmetx) %in% s2rm)
  temppathmetx = pathmetx[,-s2rm_idx]
  temppathmety = pathmety[,-s2rm_idx]
  tempsizemet = sizemet[,-s2rm_idx]
}else{
  temppathmetx = pathmetx
  temppathmety = pathmety
  tempsizemet = sizemet
}


for (pway in rownames(pathaverage)){
  
  # Get the number of studies with at least 5 measured metabolites in the pathway
  sigstudies = colnames(tempsizemet)[ which(tempsizemet[pway,] > minpathsize) ]
  pathaverage[pway,'Prop'] = mean( temppathmetx[pway,sigstudies] )
  pathaverage[pway,'DA'] = mean( temppathmety[pway,sigstudies] )
  pathaverage[pway,'NumStudies'] = length( sigstudies )
}

# Drop studies with NAN or with too few studies
PA2plot = pathaverage
PA2plot = PA2plot[-which(PA2plot$NumStudies < 5),]
PA2plot$Label = ''
signames = which(abs(PA2plot$DA)>0.2 | PA2plot$Prop > 0.4)
PA2plot[signames,'Label'] = rownames(PA2plot)[signames]

maxx = max(PA2plot$Prop)
maxy = max(abs(PA2plot$DA))
ggplot(PA2plot,aes(Prop,DA,label = Label)) + geom_point(aes(size = NumStudies)) + 
  geom_text(size = 3,position = position_jitter(width=0, height=0.04)) + 
  xlim(0.01,maxx+.2) + ylim(-maxy-0.05,maxy + 0.05) + theme_bw() +
  xlab('Average Proportion of Metabolites Differentially Abundant') + 
  ylab('Average Differential Abundance Score')  + 
  geom_abline(slope = -1) +  geom_abline() + 
  ggsave('../../results/pathway.pdf',height =8, width = 8)


###############################################################################################
# Manipulate the data arrays to get the number of mets/genes which go up/down/dont change
metup = (pathmetx + pathmety)*sizemet/2; 
metdown = (pathmetx - pathmety)*sizemet/2
geneup = (pathgenex + pathgeney)*sizegene/2; geneup[which(is.nan(geneup))] = 0
genedown = (pathgenex - pathgeney)*sizegene/2; genedown[which(is.nan(genedown))] = 0
zeromet = sizemet - metup - metdown
zerogene = sizegene - geneup - genedown; zerogene[which(is.nan(zerogene))] = 0

# Figure out the "length" of each coordinate
lengthmet = sqrt(metup^2 + metdown^2 + zeromet^2)
lengthgene = sqrt(geneup^2 + genedown^2 + zerogene^2)

# Calculate heterogeneity score (a normalized dot product)
antihetscore = (metup*geneup + metdown*genedown + zeromet*zerogene)/(lengthmet*lengthgene)
antihetscore[which(is.nan(antihetscore))] = 1
hetscore = 1-antihetscore
trimhetscore = melt( hetscore[which(rowSums(hetscore)>0),] )

# Now, generate an offset for each study to indicate heterogeneity
mergedata$X_Coord = 0
for (i in 1:length(studies)){
  mergedata[which(mergedata$Cancer == studies[[i]]),'X_Coord'] = 2*i-1
}

mergedata_wgene = mergedata[which(mergedata$Gene!=0),]

# Make scatter plot
scatter = ggplot(mergedata_wgene,aes(Met,Gene)) + geom_point(aes(color = Cancer,size = MetSize))+ 
  xlab(xlab) + ylab(ylab) + theme_classic(base_size = 20) + theme(legend.box = "horizontal") +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),
        axis.text.x = element_blank(),axis.text.y = element_blank())
legend = g_legend( scatter )
scatter = scatter + theme(legend.position = 'none')

scatter + ggsave('../../results/savedfigures/heterog/scatter.pdf',height = 10,width= 10,, useDingbats=FALSE)

# Make a pretty heatmap
mergedata$Pathway = factor(mergedata$Pathway,levels = c(rev(levels(mergedata$Pathway))))
ggplot(mergedata,aes(X_Coord + Offset,Pathway,color = Met)) + geom_point(aes(size = MetSize)) + theme_bw() + 
  scale_colour_gradient2(low = 'blue', mid = 'grey', high ='red',  midpoint = 0) + 
  theme(panel.border = element_blank(), panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = "black")) +
  geom_vline(xintercept = seq(2,length(studies)*2,by=2),linetype = "dotted") + 
  geom_segment(aes(x=X_Coord,y=Pathway,xend = X_Coord+Offset,yend = Pathway)) + xlab('Study') + 
  scale_x_discrete(breaks=seq(1,length(studies)*2,by=2), labels=studies)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  ggsave('../../results/savedfigures/metpathway.pdf',width = 9, height = 9)

# Repeat above for each study

ggplot(mergedata,aes(Cancer,Pathway,color = Gene)) + geom_point(shape = 15,size = 5) + theme_bw() + 
  scale_colour_gradient2(low = 'red', mid = 'white', high ='blue',  midpoint = 0)

# Repeat for each individual study
for (study in unique(mergedata$Cancer)){
  tempdata = mergedata[which(mergedata$Cancer == study),]
  tempdata[,2] = factor(tempdata[,2])
  tempdata = tempdata[which(tempdata$MetSize>=5),]
  ggplot(tempdata,aes(Met,Pathway,color = Met)) + geom_point(aes(size = MetSize)) + theme_bw() + 
    scale_colour_gradient2(low = 'blue', mid = 'grey', high ='red',  midpoint = 0) + 
    theme(panel.border = element_blank(), panel.grid.minor = element_blank(),  
          axis.line = element_line(colour = "black")) +
    geom_vline(xintercept = 0,linetype = "dotted") + 
    geom_segment(aes(x=0,y=Pathway,xend = Met,yend = Pathway)) + xlab('Differential Abundance Score') + 
    scale_x_continuous(limits = c(-1.1, 1.1)) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle(study) +
      ggsave(paste('../../results/savedfigures/metpathway_',study,'.pdf',sep=''),width = 9, height = 9, useDingbats=FALSE)
  
}

for (study in studies){
  plotdata = data.frame(pathmety[,study], pathgeney[,study],sizemet[,study],rownames(pathmety))
  colnames(plotdata) = c('Met','Gene','MetSize','Pathway')
  
  # Remove labels that are not significant
  plotdata$Pathway = as.character(plotdata$Pathway)
  plotdata$Pathway[which(abs(plotdata$Gene) < 0.6 & abs(plotdata$Met) < 0.6)] = ''
  plotdata$Pathway[which(plotdata$MetSize < 5)] = ''
  plotdata$Color = 'Grey'
  plotdata$Color[which(plotdata$Met >= 0.5)] = 'Red'
  plotdata$Color[which(plotdata$Met <= -0.5)] = 'Blue'
  
  plotdata = plotdata[-which(plotdata$MetSize < 5),]
  
  correlation = cor.test(plotdata[,1],plotdata[,2],method = 'spearman')
  print(study)
  print(correlation)
  
  # Write this data to a file
  write.csv(plotdata,paste('../../results/savedfigures/heterog/scatter_',study,'.csv',sep=''))
  
  print( ggplot(plotdata,aes(Met,Gene,label = Pathway, color = Color)) + 
           geom_point(aes(size = MetSize))+ 
           xlab(xlab) + ylab(ylab) + theme_classic() + 
           scale_color_manual(values = c('Blue' = 'blue','Red' = 'red','Grey' = 'grey')) +
           ggtitle(study) + geom_text(hjust=0.4, vjust=0,size = 3) +
           ggsave(paste('../../results/savedfigures/heterog/scatter_',study,'.pdf',sep=''),height = 10,width= 10,, useDingbats=FALSE))
}

###############################################################################################
# Make some pretty density plots

xhist = ggplot(mergedata_wgene,aes(Met,fill = Cancer)) + geom_density(alpha = 0.5) + xlab('') + ylab('') + 
  theme_classic() + facet_grid(Cancer ~ .) + scale_fill_discrete(guide=FALSE)+
  theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank(),
        axis.text.x = element_blank(),axis.text.y = element_blank(), 
        strip.text.x = element_blank(),strip.text.y = element_blank(),strip.background = element_blank())
xhist + ggsave('../../results/savedfigures/heterog/met.pdf',height = 10,width= 10,, useDingbats=FALSE)

yhist = ggplot(mergedata_wgene,aes(Gene,fill = Cancer)) + geom_density(alpha = 0.5) + xlab('') + ylab('') + 
  theme_classic() + facet_grid(~Cancer) + scale_fill_discrete(guide=TRUE) + coord_flip()+
  theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank(),
        axis.text.x = element_blank(),axis.text.y = element_blank(), 
        strip.text.x = element_blank(),strip.text.y = element_blank(),strip.background = element_blank())
yhist + ggsave('../../results/savedfigures/heterog/gene.pdf',height = 10,width= 10,, useDingbats=FALSE)

# Make whole figure
empty =  grid.rect(gp=gpar(col="white"))
pdf('../../results/savedfigures/heterog/fullplot.pdf',height = 10,width= 10,useDingbats=FALSE)
grid.arrange(xhist, legend, scatter, yhist, ncol=2, nrow=2, widths=c(4, 2), heights=c(2, 4))
dev.off()

# # Also calculate enrichment scores for each pathway using fisher tests and odds ratios
# ESmet = matrix(1,dim(pathmetx)[1],dim(pathmetx)[2])
# colnames(ESmet) = colnames(pathmetx)
# rownames(ESmet) = rownames(pathmetx)
# ESgene = ESmet
# 
# for (study in unique(mergedata$Cancer)){
#   m = read.csv(paste('/Users/ereznik/Documents/pancanmet/data/pathview/',study, '_met.csv',sep=''),row.names = 1, header = F)
#   inpath_sigchange = metup[,study] + metdown[,study]
#   inpath_nochange = zeromet[,study]
#   outpath_sigchange = length(which(m[,1]!=0)) - sizemet[,study]
#   outpath_nochange = length(which(m[,1] == 0)) - sizemet[,study]
#   for (i in 1:length(inpath_nochange)){
#     if (inpath_nochange[i] + inpath_sigchange[i] == 0){next}
#     testmat = matrix(c(inpath_sigchange[i], inpath_nochange[i], 
#                        outpath_sigchange[i], outpath_nochange[i]),
#                      nrow = 2,
#                      dimnames = list(Diff = c('SigChange','NoChange'),
#                                      Path = c('InPath','NotinPath')) )
#     ESmet[i,study] = fisher.test( testmat )$p.value
#   }
#   
# }
