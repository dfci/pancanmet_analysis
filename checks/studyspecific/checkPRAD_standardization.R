?# Script to compare our data standardization to Sreekumar's in Prostate Data

rm(list = ls())
setwd('/Users/ereznik/pancanmet_analysis/checks/studyspecific/')
library(XLConnect)
library(ggplot2)

old = readWorksheetFromFile('../../data/studies/PRAD/Sreekumar_etal_2009_tissue_metabolites (1).xls',sheet = 'StandardizedData')
uqnames = which(!duplicated(old[,1]))
old = old[uqnames,]
rownames(old) = old[,1]
old = old[,-1]
old = as.matrix(old)

new = read.csv('../../data/metabolomics/PRAD/PRAD_metdata.csv',row.names = 1,header = TRUE)
newrnames = rownames(new)
new = new[-1,]
new = apply(new,2,as.numeric)
rownames(new) = newrnames[2:length(newrnames)]

ixmets = intersect(rownames(new),rownames(old))
res = data.frame()
for (met in ixmets){
  
  # Correlate
  temp = cor.test(old[met,],new[met,],method = 'spearman')
  res[met,'p'] = temp$p.value
  res[met,'r'] = temp$estimate
}

ggplot(res,aes(r)) + geom_histogram()
ggplot(res,aes(-log10(p))) + geom_histogram()

# Plot a sample
met = 'Sarcosine (N-Methylglycine)'
qplot(old[met,],new[met,])
