# Script to compare metabolomics data from different studies. WE USE UNCORRECTED P-VALUES TO DO THE TESTING!

# Issues:
# 1. Check the multiplot code and the way we are generating combinations.

rm(list = ls())
analysisDir = normalizePath(".")
#setwd(file.path(analysisDir, "analysis"))
setwd('/Users/ereznik/pancanmet_analysis/analysis/')
source('functions/readBigMet.R')
source('functions/multiplot_mod.R')
source('plottingconventions.R')
library(ggplot2)
library(xlsx)
library(combinat)

# Parameters to set
study2drop = c()
pthresh = 0.1

# Read in data
metdata = readBigMet( study2drop )
met = metdata$met; tissuetype = metdata$tissuetype; studytype = metdata$studytype
uqstudy = unique(studytype)

# Read in data on the tissues each study derives from
tmap = read.csv('../data/merged_metabolomics/merged_tissuetypes.csv', row.names =1, header = T)

# Find the tumor types where we have more than one study
mult = names(which(table(tmap)>1))

pnames = list()
pctr = 1
for (ttype in mult){
  print(ttype)
  
  # Find studies with this tissue type
  s2use = rownames(tmap) [ which(tmap[,1] == ttype) ]
  
  # Reduce the data to these types
  idx = which(studytype %in% s2use)
  d = met[,idx]
  d = d[complete.cases(d),]
  d = as.matrix(d)
  
  # Get study and tissuetypes for this data
  d_study = unlist( lapply(strsplit(colnames(d),':'), `[[`, 1) )
  d_tissue = unlist( lapply(strsplit(colnames(d),':'), `[[`, 3) )
  
  # Make a data frame to store results
  fc = data.frame( matrix(0,dim(d)[1],length(s2use)), row.names = rownames(d) ); colnames(fc) = s2use
  p = fc
  padj = fc
  
  for ( s in unique(s2use) ){
    
    # Get normal and tumor idx
    normidx = which(d_study == s & d_tissue == 'Normal')
    tumidx = which(d_study == s & (d_tissue == 'Tumor'  ))
    
    for (metabolite in rownames(d)){
      fc[metabolite,s] = log2( mean(d[metabolite,tumidx]) / mean(d[metabolite,normidx]) )
      p[metabolite,s] = wilcox.test( d[metabolite,tumidx], d[metabolite,normidx])$p.value
    }
    
    padj[,s] = p.adjust(p[,s],method = 'BH')
  }
  
  write.xlsx2(fc, paste('../results/sametissue/',ttype,'.xlsx',sep = ''), sheetName='FC',
              col.names=TRUE, row.names=TRUE, append=FALSE)
  
  write.xlsx2(p, paste('../results/sametissue/',ttype,'.xlsx',sep = ''), sheetName='P',
              col.names=TRUE, row.names=TRUE, append=TRUE)
  
  write.xlsx2(padj, paste('../results/sametissue/',ttype,'.xlsx',sep = ''), sheetName='PAdj',
              col.names=TRUE, row.names=TRUE, append=TRUE)
  
  # Make an FC table with insignificant changes set to zero
  fcadj = fc
  fcadj[p > pthresh] = 0
  
  # Generate all permutations and plot
  scombs = combn(s2use,2)
  if (length(scombs)==2){scombs = matrix(t(scombs))}
  
  for (i in 1:dim(scombs)[2]){
    s2plot = scombs[,i]
    plotdata = fc[,s2plot]
    plotdata[which(padj[,s2plot] > pthresh,arr.ind = T)] = 0
    colnames(plotdata) = c('X','Y')
    
    # Calculate the p value of a fisher test
    tempfc = fcadj[,s2plot]
    tempfc = tempfc[-which(tempfc == 0, arr.ind = TRUE)[,1],]
    tempfc = sign(tempfc)
    
    pvalue = fisher.test(table(tempfc))$p.value
    pvalue = round(pvalue,3)
    
    fignum = ggplot(plotdata,aes(X,Y)) + geom_point() + theme_classic() + 
      xlab(names2plot[s2plot[1]]) + ylab(names2plot[s2plot[2]])  + 
      geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + 
      ggtitle(paste(ttype,', Fisher Test P-Value = ',pvalue,sep = ''))
    pname = paste('fignum',pctr,sep='')
    pctr = pctr + 1
    pnames = c(pnames,pname)
    assign(pname,fignum)
    print( fignum )
  }
  
}

# Make a multiplot
pnames = unlist(pnames)
pdf('../results/sametissue/sametissue.pdf',height = 10,width = 10)
multiplot_mod(pnames,cols = 2)
dev.off()
