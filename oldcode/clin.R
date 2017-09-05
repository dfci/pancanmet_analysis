# Correlate metabolites with stage

# Script to analyze metastases in reporter data and identify common features

rm(list=ls())
setwd('/Users/ereznik/Documents/pancanmet/analysis/pancan/')
library(ggplot2)

# Set significance threshold
thresh = 0.05

# Read in all data, isolate metastases
alldata = read.csv('/Users/ereznik/Documents/pancanmet/data/merged_metabolomics/alldata.csv', row.names = 1,check.names = F)
m = data.frame( sapply(alldata[3:dim(alldata)[1],], function(x) { as.numeric(as.character(x)) } ), check.names = F )
rownames(m) = rownames(alldata)[3:dim(alldata)[1]]
splittypes = strsplit(colnames(alldata),'\\:')
tisstype = sapply(splittypes, "[", 3)
cancertype = alldata[1,1:(dim(alldata)[2])]
cancertype = sapply(unlist(cancertype),as.character)
uqstudy = unique(cancertype)

# Read in the clinical file
clin = read.csv('../../data/merged_metabolomics/clinfeatures.csv',header = T)
for (i in 1:dim(clin)[1]){
  rownames(clin)[i] = paste(clin[i,2],clin[i,1],clin[i,3],sep= ':')
}

s2use = c('BRCA','BRCATang','KIRC','OV','PRAD','PRADLODA')
for (s in s2use){
  
  # Check to see if we have data
  ix = which(cancertype == s)
  ixnames = intersect(colnames(m)[ix],rownames(clin))
  
  if (length(ixnames)<5){next}
  
  d1 = m[,ixnames]
  d1 = d1[complete.cases(d1),]
  d1 = t(d1)
  
  d2 = data.frame( sapply(clin[ixnames,4:5], function(x) { as.numeric(as.character(x)) } ), check.names = F )
  rownames(d2) = ixnames
  
  # Make some arrays
  r = matrix(0,dim(d1)[2],2)
  rownames(r) = colnames(d1)
  colnames(r) = c('Stage','Grade')
  p = matrix(1,dim(d1)[2],2)
  rownames(p) = colnames(d1)
  colnames(p) = c('Stage','Grade')
  
  for (ctype in c('Stage','Grade')){
    if (is.na(d2[1,ctype])){next}
    for (met in 1:dim(d1)[2]){
      temp = cor.test( d1[,met], d2[,ctype],method = 'kendall')
      p[met,ctype] = temp$p.value
      r[met,ctype] = temp$estimate
    }
  }
  padj = p
  padj[,1] = p.adjust(p[,1],method = 'BH')
  padj[,2] = p.adjust(p[,2],method = 'BH')
  
  alldata = cbind(r,padj)
  colnames(alldata) = c('RStage','RGrade','PStage','PGrade')
  write.csv(alldata,paste('../../results/clinassoc/',s,'.csv',sep=''))
}