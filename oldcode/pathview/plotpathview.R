# Plot integrated metabolite/expression data onto pathways with Pathview
rm(list=ls())
setwd('/Users/ereznik/Documents/reporter/bin/pathview/')
#setwd('/cbio/cslab/home/reznik/reporter/analysis/pathview/') # For the cluster

library(pathview)
library(reshape)
library(ggplot2)
library(KEGGREST)
library(gridExtra)
printpath = '/Users/ereznik/Documents/reporter/results/pathview/'
datapath = '/Users/ereznik/Documents/reporter/data/pathview/'
source('rings2.R')

# Write a function which grabs the legend
g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    legend
}

# Read in KEGG pathway names
maxpath = 92
keggdbpaths = keggList("pathway", "hsa") 
paths = sapply(strsplit(names(keggdbpaths[1:maxpath]),':'), "[", 2)
pathnames = sapply(strsplit(keggdbpaths[1:maxpath],' - '), "[", 1)
keggnames = data.frame(pathnames); rownames(keggnames) = paths

studies = c('OV','BRCApos','BRCAneg','KIRC','ccpap','BLCA','PRAD','STAD','COAD')
#studies = c('BRCApos','BRCAneg')
studies = c('KICH')

pathmetx = matrix(0,dim(keggnames)[1],length(studies)); rownames(pathmetx) = keggnames[paths,1]; colnames(pathmetx) = studies
pathgenex = matrix(0,dim(keggnames)[1],length(studies)); rownames(pathgenex) = keggnames[paths,1]; colnames(pathgenex) = studies
pathmety = matrix(0,dim(keggnames)[1],length(studies)); rownames(pathmety) = keggnames[paths,1]; colnames(pathmety) = studies
pathgeney = matrix(0,dim(keggnames)[1],length(studies)); rownames(pathgeney) = keggnames[paths,1]; colnames(pathgeney) = studies
sizemet = matrix(0,dim(keggnames)[1],length(studies)); rownames(sizemet) = keggnames[paths,1]; colnames(sizemet) = studies

for (study in studies){
  if ( file.exists( paste(datapath,study,'_gene.csv',sep='') ) ){
    isg = TRUE
    g = read.csv(paste(datapath,study,'_gene.csv',sep=''),row.names = 1, header = F)
  } else{isg = FALSE}
  
  # Create the directories we are printing to
  outdirs = c(paste(printpath,study,sep=''),paste(printpath,study,'/metabolograms',sep=''))
  for (outdir in outdirs){
   dir.create(file.path( outdir ), showWarnings = FALSE)
   do.call(file.remove,list(list.files(outdir,full.names=TRUE)))
  }
 
  m = read.csv(paste(datapath,study, '_met.csv',sep=''),row.names = 1, header = F)
  
  for (path in paths){
    setwd(paste(printpath,study,sep=''))
    if (isg){genes = g[,1]; names(genes) = rownames(g)} 
    mets = m[,1]; names(mets) = rownames(m)
    pathname = pathnames[paste('path',path,sep=':')]
    
    # If you want to plot, uncomment below
    if (isg){
      df = pathview(gene.data = genes, cpd.data = mets, pathway.id =  path,
                    species = "hsa", out.suffix = study, kegg.native = T,
                    low = list(gene = 'green',cpd = 'green'),high = list(gene = 'red',cpd = 'red'),
                    limit = list(gene=5,cpd = 3), bins = list(gene=30, cpd=30),
                    kegg.dir = paste(printpath,'keggdir',sep=''))
    } else{
      df = pathview(cpd.data = mets, pathway.id =  path,
                    species = "hsa", out.suffix = study, kegg.native = T,
                    low = list(gene = 'green',cpd = 'green'),high = list(gene = 'red',cpd = 'red'),
                    limit = list(gene=5,cpd = 3), bins = list(gene=30, cpd=30),
                    kegg.dir = paste(printpath,'keggdir',sep=''))
    }
   
    
    # Additionally, write the study and pathway specific genes and metabolites to a table
    if (!is.list(df)){next}
    intm =  intersect(df$plot.data.cpd[,1],rownames(m))
    write.csv(m[intm,],paste(printpath,study,paste(study,path,'met.csv',sep='_'), sep='/' ) )
    
    if (is.list(df)){
      if (isg){
        intg = intersect(df$plot.data.gene[,1],rownames(g))
        write.csv(g[intg,],paste(printpath, study,paste(study,path,'gene.csv',sep='_'), sep='/' ) )
        numgene = length(which(abs(df$plot.data.gene$mol.data)>0)) /
          length(which(!is.na(df$plot.data.gene$mol.data)))
        signgene = ( length(which(df$plot.data.gene$mol.data>0)) - length(which(df$plot.data.gene$mol.data<0))) /
          length(which(!is.na(df$plot.data.gene$mol.data)))
        pathgenex[as.character(keggnames[path,1]),study] = numgene
        pathgeney[as.character(keggnames[path,1]),study] = signgene
        
        # If we have gene data make a ring plot
        if (length(intm)>0){
          rings2( list(data.frame(g[intg,1]),data.frame(m[intm,1])), pathname, paste(printpath,study,
              paste('metabolograms/',study,gsub("[[:punct:][:space:]]", "", pathname),'metabologram.pdf',sep='_'), sep='/' ) )
        }
      }
      
      if (!isg & length(intm)>0){
        # Make a ring plot using only metabolite data
        rings2( list(data.frame(0),data.frame(m[intm,1])), pathname, paste(printpath,study,
        paste('metabolograms/',study,gsub("[[:punct:][:space:]]", "", pathname),'metabologram.pdf',sep='_'), sep='/' ) )
      }
     
      if (length(df$plot.data.cpd) == 0){next}
      nummet =  length(which(abs(df$plot.data.cpd$mol.data)>0)) /
        length(which(!is.na(df$plot.data.cpd$mol.data)))
      signmet = ( length(which(df$plot.data.cpd$mol.data>0)) - length(which(df$plot.data.cpd$mol.data<0))) /
        length(which(!is.na(df$plot.data.cpd$mol.data)))
      pathmetx[as.character(keggnames[path,1]),study] = nummet
      pathmety[as.character(keggnames[path,1]),study] = signmet
      sizemet[as.character(keggnames[path,1]),study] = length(which(!is.na(df$plot.data.cpd$mol.data)))
    }
    
  }
}

#save(pathmetx,pathmety,pathgenex,pathgeney,sizemet,file = '/Users/ereznik/Documents/reporter/data/pathview/allpathviewdata.Rda')

load('/Users/ereznik/Documents/reporter/data/pathview/allpathviewdata.Rda')

# Remove pathways with average size less than X
p2rm = which( rowMeans(sizemet) < 5 )

# Plot the result
ylab = '<-- Genes Lower in Tumor || Genes Higher in Tumor -->'
xlab = '<-- Metabolites Lower in Tumor || Metabolites. Higher in Tumor -->'
metmelt = melt(pathmety[-p2rm,]); colnames(metmelt) = c('Pathway','Cancer')
metmeltx = melt(pathmetx[-p2rm,]); colnames(metmeltx) = c('Pathway','Cancer')
genemelt = melt(pathgeney[-p2rm,]); colnames(genemelt) = c('Pathway','Cancer')
sizemelt = melt(sizemet[-p2rm,]); colnames(sizemelt) = c('Pathway','Cancer')
mergedata = merge(metmelt,genemelt,by=c('Pathway','Cancer')); colnames(mergedata) = c('Pathway','Cancer','Met','Gene')
mergedata = merge(mergedata,metmeltx,by=c('Pathway','Cancer')); colnames(mergedata) = c('Pathway','Cancer','Met','Gene','Offset')
mergedata = merge(mergedata,sizemelt,by=c('Pathway','Cancer')); colnames(mergedata) = c('Pathway','Cancer','Met','Gene','Offset','MetSize')

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
    axis.text.x = element_blank(),axis.text.y = element_blank(), plot.margin = unit(c(3,-5.5,4,3), "mm"))
legend = g_legend( scatter )
scatter = scatter + theme(legend.position = 'none')
    
scatter + ggsave('../../results/savedfigures/heterog/scatter.pdf',height = 10,width= 10,, useDingbats=FALSE)

# Make a pretty heatmap
ggplot(mergedata,aes(X_Coord + Offset,Pathway,color = Met)) + geom_point(aes(size = MetSize)) + theme_bw() + 
  scale_colour_gradient2(low = 'red', mid = 'white', high ='blue',  midpoint = 0) + 
  theme(panel.border = element_blank(), panel.grid.minor = element_blank(),  
          axis.line = element_line(colour = "black")) +
  geom_vline(xintercept = seq(2,length(studies)*2,by=2),linetype = "dotted") + 
  geom_segment(aes(x=X_Coord,y=Pathway,xend = X_Coord+Offset,yend = Pathway)) + xlab('Study') + 
  scale_x_discrete(breaks=seq(1,length(studies)*2,by=2), labels=studies)


ggplot(mergedata,aes(Cancer,Pathway,color = Gene)) + geom_point(shape = 15,size = 5) + theme_bw() + 
  scale_colour_gradient2(low = 'red', mid = 'white', high ='blue',  midpoint = 0)

# Repeat for each individual study

for (study in studies){
  plotdata = data.frame(pathmety[,study], pathgeney[,study],sizemet[,study],rownames(pathmety))
  colnames(plotdata) = c('Met','Gene','MetSize','Pathway')
  print( ggplot(plotdata,aes(Met,Gene,label = Pathway)) + geom_point(aes(size = MetSize))+ 
    xlab(xlab) + ylab(ylab) + theme_classic() + ggtitle(study) + geom_text(hjust=0, vjust=0,aes(size = MetSize)))
}

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
pdf('../../results/savedfigures/heterog/fullplot.pdf',height = 10,width= 10)
grid.arrange(xhist, legend, scatter, yhist, ncol=2, nrow=2, widths=c(4, 2), heights=c(2, 4))
dev.off()
