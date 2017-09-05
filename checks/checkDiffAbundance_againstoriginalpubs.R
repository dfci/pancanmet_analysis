# Check differential abundance against original publications

# Notes:
# 1. For ovarian, they only compare EOC (non-metastatic tumors) to normals, so the fold changes are a little off. 
# If we run diff abundance summary with only tumors and not metastases, we get the exact correct fold change back.


rm(list = ls())
setwd('/Users/ereznik/pancanmet_analysis/checks/')
library(ggplot2)
library(cowplot)
library(XLConnect)

# Set studies to check and studies where we need to log-normalize fold changes
s2check = c('PAADHussain1','PAADHussain2','KIRC','OV')
#s2check = c('PAADHussain1') # for testing
takefold = c('PAADHussain1','PAADHussain2','OV') # for these studies, take the log2 of the fold change

# Read in big table of results
bigres = read.csv('../results/diffabundance/DifferentialAbundanceSummary.csv',header = TRUE,row.names = 1)

# Read in merged mapping
mapping = read.csv('../import/tempdir/merged_mapping.csv',header = TRUE,row.names = 1, stringsAsFactors = FALSE)

# Get names of metabolites
finalnames = read.csv('../import/tempdir/FinalMetaboliteNames.csv',header = TRUE,row.names = 1)
rownames(mapping) = finalnames[rownames(mapping),1]

# Read in original publication results
orig = XLConnect::loadWorkbook(filename = 'checkdata/DiffAbundance_OriginalPubs.xlsx')

for (s in s2check){
  
  # Read in original data
  origd = readWorksheet(orig,sheet = s,rownames = 1,header = TRUE)
  rownames(origd) = sapply(rownames(origd),tolower)
  
  # Make sure it's numeric class
  origd$Q = as.numeric(origd$Q)
  
  if (s %in% takefold){
    origd$FC = log2(origd$FC)
  }
  
  # Get new data
  newd = data.frame( bigres[,c(paste('FC',s,sep = '.'),paste('Padj',s,sep = '.'))])
  colnames(newd) = c('FC','Q')
  newd = newd[complete.cases(newd),]
  
  # Rename rows
  rownames(newd) = mapping[rownames(newd),s]
  rownames(newd) = sapply(rownames(newd),tolower)
  
  # Intersect. The metabolites missing in the intersection were probably those that were dropped because of NA 
  # The drop happens (probably) because Wilcoxen test fails with 100% imputation in either tumor/normal group
  ixmets = intersect(rownames(origd),rownames(newd))
  print(length(ixmets)/dim(newd)[1])
  print(setdiff(rownames(origd),rownames(newd)))
  print(setdiff(rownames(newd),rownames(origd)))
  
  p1 = qplot(newd[ixmets,'FC'],origd[ixmets,'FC']) + theme_bw() + geom_point() + 
    geom_abline() + xlab('New Calc') + ylab('Orig Pub') + ggtitle(paste('FC',s)) + 
    geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
  p2 = qplot( log10(newd[ixmets,'Q']), log10(origd[ixmets,'Q']) ) + theme_bw() + geom_point() + 
    geom_abline() + xlab('New Calc') + ylab('Orig Pub') + ggtitle(paste('Q',s)) + 
    geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
  
  pall = plot_grid(p1,p2)
  print(pall)
}