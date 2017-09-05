# Script to check median absolute deviation code.

rm(list = ls())
setwd('/Users/ereznik/pancanmet_analysis/checks/')
library(ggplot2)

# Set the variatino function we want to use
variation = mad # mad or sd

# Read in the big data matrix and the paired data file
met = as.matrix( read.csv('../data/merged_metabolomics/merged_metabolomics.csv',header = TRUE,row.names = 1,check.names = FALSE) )
pairs = read.csv('../data/merged_metabolomics/tumor_normal_pairs.csv',header = TRUE,stringsAsFactors = FALSE)
pairs$NormalName = paste(pairs$Study,pairs$Normal.Sample,'Normal',sep = ':')
pairs$TumorName = paste(pairs$Study,pairs$Tumor.Sample,'Tumor',sep = ':')

uqstudies = unique(pairs$Study)
metres = data.frame('tumor' = numeric(0),'normal' = numeric(0))
for (s in uqstudies){
  print(s)
  pairix = pairs[which(pairs$Study == s),]
  m2 = met[,pairix$NormalName] # only need this to get the metabolites we should test
  m2 = m2[complete.cases(m2),]
  mets2test = rownames(m2)
  
  for (metname in mets2test){
    
    metres[paste(s,metname,sep=':'),'normal'] = variation( met[metname,pairix$NormalName] )
    metres[paste(s,metname,sep=':'),'tumor'] = variation( met[metname,pairix$TumorName] )
    metres[paste(s,metname,sep=':'),'Study'] = s
  }
  
}

metres$Ratio = log2(metres$tumor/metres$normal)
metres = metres[which(!is.na(metres$Ratio)),]
metres = metres[which(!is.infinite(metres$Ratio)),]

ggplot(metres,aes(Ratio)) + 
  theme_bw() + geom_vline(xintercept = 0) + ylab('Number of Metabolites') + 
  xlab('Log2 Ratio of Median Absolute Deviation, Tumor/Normal') +
  geom_histogram(aes(fill = Study),alpha = 0.5,binwidth = 0.2) + 
  geom_density(aes(y=0.2*..count..,fill = Study),alpha=0.5) +
  facet_grid(Study~.,scales = 'free',drop = FALSE) +
  theme_classic() + geom_vline(xintercept = 0,linetype = "longdash") + 
  theme(legend.position = 'none') +
  theme(strip.background = element_rect(colour="white", fill="white"))

