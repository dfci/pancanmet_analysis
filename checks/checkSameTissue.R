# Script to check RunSameTissue.R is correct. This uses our diff abundance calculations in the saved table. 
# Note that the p-value threshold set in RunSameTissue might be different than that used for diffabundance calculations (e.g. 0.1 vs 0.05)

rm(list = ls())
setwd('/Users/ereznik/pancanmet_analysis/checks/')
library(ggplot2)
library(combinat)
library(cowplot)

# Set p value threshold
pthresh = 5e-2

# Read in data on the tissues each study derives from
tmap = read.csv('../data/merged_metabolomics/merged_tissuetypes.csv', row.names =1, header = T)

# Read in the results of our analysis
da = read.csv('../results/diffabundance/DifferentialAbundanceSummary.csv',header = TRUE,row.names= 1)

# Find the tumor types where we have more than one study
mult = names(which(table(tmap)>1))


for (ttype in mult){
  print(ttype)
  s2use = rownames(tmap)[tmap$Tissue == ttype]
  scombs = combn(s2use,2)
  if (length(scombs)==2){scombs = matrix(t(scombs))}
  
  for (i in 1:dim(scombs)[2]){
   s1 = scombs[1,i]
   s2 = scombs[2,i]
   
   # Read in FC
   fc = da[,c(paste('FC.',s1,sep=''),paste('FC.',s2,sep = ''))]
   p = da[,c(paste('Padj.',s1,sep=''),paste('Padj.',s2,sep = ''))]
   
   fc = fc[which(complete.cases(fc)),]
   p = p[which(complete.cases(p)),]
   p = p[rownames(fc),]
   
   fcadj = fc
   fcadj[p > pthresh] = 0
   
   colnames(fc) = c('S1','S2')
   colnames(fcadj) = c('S1','S2')
   
   p1 = ggplot(fc,aes(S1,S2)) + geom_point() + theme_bw() + 
     xlab(s1) + ylab(s2) + ggtitle(paste(ttype,'FC')) + geom_abline() + 
     geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
   
   p2 = ggplot(fcadj,aes(S1,S2)) + geom_point() + theme_bw() + 
     xlab(s1) + ylab(s2) + ggtitle(paste(ttype,'FCAdj')) + geom_abline() + 
     geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
   
   pall = plot_grid(p1,p2)
   print(pall)
  }
  
}