# Script to compare metabolomics data from different studies

rm(list = ls())
setwd('/Users/ereznik/Documents/pancanmet/analysis/pancan/')
library(ggplot2)
library(data.table)
library(xlsx)
library(combinat)

# Source multiplot to make a grid of ggplots
source('../rfunc//multiplot_mod.R')

# Set p-value threshold
pthresh = 0.2

temp = fread('../../data/merged_metabolomics/alldata.csv',sep = ',',header = F,skip = 3)
alldata = data.frame( temp, row.names = 1 )
title.line = strsplit( readLines('../../data/merged_metabolomics/alldata.csv', n=1), ',' ) [[1]] # Note that first element is empty ""
colnames(alldata) = title.line[2:length(title.line)]

# Extract relevant data from columns
study = unlist( lapply(strsplit(colnames(alldata),':'), `[[`, 1) )

tmap = read.csv('../../data/merged_metabolomics/merged_tissuetypes.csv', row.names =1, header = T)

# Find the tumor types where we have more than one study
mult = names(which(table(tmap)>1))

pnames = list()
pctr = 1
for (ttype in mult){
  print(ttype)
  
  # Find studies with this tissue type
  s2use = rownames(tmap) [ which(tmap[,1] == ttype) ]
  
  # Reduce the data to these types
  idx = which(study %in% s2use)
  d = alldata[,idx]
  d = d[complete.cases(d),]
  d = as.matrix(d)
  
  # Make a data frame to store results
  fc = data.frame( matrix(0,dim(d)[1],length(s2use)), row.names = rownames(d) ); colnames(fc) = s2use
  p = data.frame( matrix(1,dim(d)[1],length(s2use)), row.names = rownames(d)); colnames(p) = s2use
  padj = data.frame( matrix(1,dim(d)[1],length(s2use)), row.names = rownames(d)); colnames(padj) = s2use
  
  studytype = unlist( lapply(strsplit(colnames(d),':'), `[[`, 1) )
  samptype = unlist( lapply(strsplit(colnames(d),':'), `[[`, 3) )
  
  for ( s in unique(studytype) ){
    # Get normal and tumor idx
    normidx = which(studytype == s & samptype == 'Normal')
    tumidx = which(studytype == s & samptype == 'Tumor')
    
    for (met in rownames(d)){
      fc[met,s] = log2( mean(d[met,tumidx]) / mean(d[met,normidx]) )
      p[met,s] = wilcox.test( d[met,tumidx], d[met,normidx])$p.value
    }
    
    padj[,s] = p.adjust(p[,s],method = 'BH')
  }
  
  write.xlsx2(fc, paste('../../results/sametissue/',ttype,'.xlsx',sep = ''), sheetName='FC',
              col.names=TRUE, row.names=TRUE, append=FALSE)
  
  write.xlsx2(p, paste('../../results/sametissue/',ttype,'.xlsx',sep = ''), sheetName='P',
              col.names=TRUE, row.names=TRUE, append=TRUE)
  
  write.xlsx2(padj, paste('../../results/sametissue/',ttype,'.xlsx',sep = ''), sheetName='PAdj',
              col.names=TRUE, row.names=TRUE, append=TRUE)
  
  # Generate all permutations and plot
  scombs = combn(s2use,2)
  if (length(scombs)==2){scombs = matrix(t(scombs))}
  for (i in 1:dim(scombs)[2]){
    s2plot = scombs[,i]
    plotdata = fc[,s2plot]
    plotdata[which(padj[,s2plot] > pthresh,arr.ind = T)] = 0
    colnames(plotdata) = c('X','Y')
    
    # Also calculate a p-value using fishers exact test
    d2test = plotdata[ -which(abs(plotdata)< 1e-10,arr.ind = T)[,1],]
    if (dim(d2test)[1] >1 ){
      t2test = table(sign(d2test))
      p2print = round( fisher.test(t2test)$p.value, digits = 3)
    }else{
      p2print = 1
    }
    p = ggplot(plotdata,aes(X,Y)) + geom_point() + theme_classic() + xlab(s2plot[1]) + 
      ylab(s2plot[2]) + ggtitle(paste(ttype, 'Fisher Exact P-Value',p2print)) + geom_vline() + geom_hline()
    pname = paste('p',pctr,sep='')
    pctr = pctr + 1
    pnames = c(pnames,pname)
    assign(pname,p)
    print( p )
  }
  
}

# Make a multiplot
pnames = unlist(pnames)
pdf('/Users/ereznik/Documents/pancanmet/results/sametissue.pdf',height = 10,width = 10)
multiplot_mod(pnames,cols = 2)
dev.off()
