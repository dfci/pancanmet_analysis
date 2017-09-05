# A script that reads in differential abundance calculations for each study and maps to KEGG pathways for analysis.

rm(list=ls())
setwd('/Users/ereznik/Documents/pancanmet/analysis/pancan/')

#library(reshape)
library(ggplot2)
library(KEGGREST)
library(reshape2)

################
# Set parameters
maxpath = 92
datapath = '/Users/ereznik/Documents/pancanmet/data/pathview/' # Directory where the differential abundances are stored
s2rm = NA#c('BRCA','BRCATang') # Studies to exclude because of biases, use NA to include everything
minpathsize = 5#6 # minimum number of metabolites in a pathway for it to be included
minstudies = 6 # Minimum number of studies with at least minpathsize metabolites in a pathway for pathway to be included

################
# Get KEGG pathway names
keggdbpaths = keggList("pathway", "hsa") 
keggdbpaths = keggdbpaths[1:maxpath]
paths = sapply(strsplit(names(keggdbpaths[1:length(keggdbpaths)]),':'), "[", 2)
pathnames = sapply(strsplit(keggdbpaths[1:length(keggdbpaths)],' - Homo'), "[", 1)
keggnames = data.frame(as.character(pathnames)); rownames(keggnames) = paths

################
# Read the table of KEGG pathway annotations
load('/Users/ereznik/Documents/pancanmet/data/KEGGmetmap.RData')
load('/Users/ereznik/Documents/pancanmet/data/KEGGgenemap.RData')

studies = c('BLCA','BRCA','BRCATang','COAD','KIRC','OV','PAAD','PAADHussain1','PAADHussain2','PRAD','PRADLODA','STAD')
#studies = c('KIRC')

################
# Initialize some data arrays
pathmetx = matrix(0,dim(keggnames)[1],length(studies)); rownames(pathmetx) = keggnames[paths,1]; colnames(pathmetx) = studies
pathmety = matrix(0,dim(keggnames)[1],length(studies)); rownames(pathmety) = keggnames[paths,1]; colnames(pathmety) = studies
sizemet = matrix(0,dim(keggnames)[1],length(studies)); rownames(sizemet) = keggnames[paths,1]; colnames(sizemet) = studies

################
# Run analysis

# Make an array to store all the KEGG IDs and their changes for each cancer type
metres = data.frame(KEGGID = character(0),study = character(0),pathway = character(0),change = numeric(0))

for (study in studies){
  print(paste('Analyzing study', study,'...') )
  
  m = read.csv(paste(datapath,study, '_met.csv',sep=''),row.names = 1, header = F)
  g = read.csv(paste(datapath,study, '_gene.csv',sep=''),row.names = 1, header = F)
  for (path in paths){
    pname = as.character( keggnames[path,1] )
    keggids = names( metmap[[path]] )
    intm = intersect(rownames(m),keggids)
    
    gene_keggids = gmap[[path]]
    
    # Alternatively, you can read in a file of gene IDs
    #gene_keggids = as.character(read.csv('/Users/ereznik/Documents/temp/mitochondrion.csv')[,1])
    
    intg = intersect(g[,2],gene_keggids)
    gidx = which(g[,2] %in% intg)
    dg = g[gidx,]
    
    d = m[intm,]
    
    # Calculate stuff
    numup = length( which(d[,1] > 0) )
    numdown = length( which(d[,1] < 0) )
    size = length(intm)
    
    if (size == 0){next}
    pathmety[pname,study] = (numup - numdown)/size
    pathmetx[pname,study] = (numup + numdown)/size
    sizemet[pname,study] = size
    tempframe = data.frame(rep(study,length(intm)), rep(path,length(intm)),intm,sign(d[intm,1]) )
    metres = rbind(metres, tempframe )
  }
  
}

# Rename colnames of metres
colnames(metres) = c('Study','Pathway','KEGGID','Change')

# Make an "averaged" figure over all studies indicating differential abundance score
pathaverage = data.frame( matrix(0,dim(pathmetx)[1],3))
rownames(pathaverage) = rownames(pathmetx)
colnames(pathaverage) = c('Prop','DA','NumStudies')

# Make some temporary data that removes the BRCA studies
if (!is.na(s2rm[1])){
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
s2drop = which(PA2plot$NumStudies < minstudies)
if (length(s2drop) > 0){
  PA2plot = PA2plot[-s2drop,]
}
PA2plot$Label = ''
signames = which(abs(PA2plot$DA)>0.2 | PA2plot$Prop > 0.35)
PA2plot[signames,'Label'] = rownames(PA2plot)[signames]

# Make a bubble plot
maxx = max(PA2plot$Prop)
maxy = max(abs(PA2plot$DA))
ggplot(PA2plot,aes(Prop,DA,label = Label)) + geom_point(aes(size = NumStudies)) + 
  geom_text(size = 3,position = position_jitter(width=0, height=0.04)) + 
  xlim(0.01,maxx+.2) + ylim(-maxy-0.05,maxy + 0.05) + theme_bw() +
  xlab('Average Proportion of Metabolites Differentially Abundant') + 
  ylab('Average Differential Abundance Score')  + 
  geom_abline(slope = -1) +  geom_abline() + 
  ggsave('../../results/KEGG/pathway_bubbleplot.pdf',height =10, width = 10)

# Make a summarizing heatmap
path2keep = rownames(PA2plot)
xmelt = melt(pathmetx[path2keep,])
rownames(xmelt) = paste(xmelt[,1],':',xmelt[,2])
ymelt = melt(pathmety[path2keep,])
rownames(ymelt) = paste(ymelt[,1],':',ymelt[,2])
sizemelt = melt(sizemet[path2keep,])
rownames(sizemelt) = paste(sizemelt[,1],':',sizemelt[,2])
heatmapdata = sizemelt
heatmapdata$Color = ymelt[rownames(heatmapdata),3]
colnames(heatmapdata) = c('Pathway','Study','Size','Color')
heatmapdata$Pathway = as.character(heatmapdata$Pathway)

# Make an extra column in the in the heatmap that indicates the number of times a pathway is >50% up/down
uqp = unique(heatmapdata$Pathway)
heatmaporder = data.frame(ScoreX = numeric(0),ScoreColor = numeric(0))
for (pathway in uqp){
  sigup = length(which(heatmapdata$Pathway == pathway & heatmapdata$Color > 0))
  sigdown = length(which(heatmapdata$Pathway == pathway & heatmapdata$Color < 0))
  numall = length(which(heatmapdata$Pathway == pathway))
  scorex = (sigup + sigdown )/numall
  pathix = which(heatmapdata$Pathway == pathway)
  scorecolor = sum(heatmapdata[pathix,'Color']*heatmapdata[pathix,'Size'])/sum(heatmapdata[pathix,'Size'])
  heatmaporder[pathway,] = c(scorex,scorecolor)
}

plotorder = rownames(heatmaporder)[order(heatmaporder$ScoreColor)]
heatmapdata$Pathway = factor(heatmapdata$Pathway,levels = plotorder)
heatmaporder$Pathway = factor(rownames(heatmaporder),levels = plotorder)

# Plot the heatmap
ggplot(heatmapdata,aes(Study,Pathway,color = Color)) + geom_point(aes(size = Size),shape=19) + 
  geom_point(data = subset(heatmapdata,Color == 0),aes(size = Size),color = 'grey50',shape = 1)+ 
  geom_point(data = subset(heatmapdata,Size == 0),shape =4,color = 'grey50',size = 5) +
  scale_color_gradient2(low="blue",high="red") + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_size_area() + theme(legend.position="bottom") + 
  ggsave('/Users/ereznik/Documents/pancanmet/results/KEGG/pathway_heatmap.pdf',height = 10,width = 10,useDingbats = FALSE)

# Make an associated plot for each pathway
pplot = ggplot(heatmaporder,aes(ScoreColor,Pathway,color = ScoreColor)) + geom_point(shape = 15,size = 3) + theme_classic(20) + 
  scale_color_gradient2(low = 'blue',mid = 'grey50','high' = 'red') + xlim(-1,1) + 
  geom_segment(aes(x = 0, y = Pathway, xend = ScoreColor, yend = Pathway),color = 'black',linetype = 2) +
  theme(legend.position="none")

pplot + scale_y_discrete(breaks=NULL) + xlab('') + ylab('') + 
  ggsave('/Users/ereznik/Documents/pancanmet/results/KEGG/pathway_bars.pdf',height = 10,width = 4,useDingbats = FALSE)

