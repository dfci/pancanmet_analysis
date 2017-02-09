# Script to visualize cancer types by heterogeneity

rm(list = ls())
setwd('/Users/ereznik/pancanmet_analysis/analysis/')
source('plottingconventions.R')
library(ggplot2)
library(ggrepel)

# Read in differential abundance summary
da = read.csv('../results/diffabundance/DifferentialAbundance_SignedChanges.csv',header = TRUE,row.names = 1)

# Read in MAD results
madres = read.csv('../results/variation/MAD_results.csv',row.names = 1,header = TRUE)

# Augment MAD results
madres$NumDiff = NA
madres$DAScore = NA
for (study in rownames(madres)){
  
  temp = da[,study,drop = FALSE]
  temp = temp[which(!is.na(temp)),1,drop = FALSE]
  madres[study,'NumDiff'] = length(which( temp!=0 ))/dim(temp)[1]
  madres[study,'DA'] = (length(which(temp > 0)) - length(which(temp < 0)) )/length(which(!is.na(temp)))
  
}

# Make variation score
madres$Variation = madres$MeanMADRatio
madres[which(madres$P > pthresh),'Variation'] = 0

# Annotate
madres$Tissue = sourcetissue[rownames(madres)]
madres$Name = names2plot[rownames(madres)]

ggplot(madres,aes(Variation,DA,color = Tissue,label = Name)) + geom_point(aes(size = NumDiff)) + 
  ylim(-1,1) + scale_size_continuous(range = c(5,10)) + theme_bw(base_size = 18) + 
  xlab('Metabolic Variation, Tumor/Normal') + ylab('Differential Abundance') + 
  geom_vline(xintercept = 0,linetype = 'dashed') + geom_hline(yintercept = 0,linetype = 'dashed') + 
  scale_size_continuous(name = 'Proportion\nDifferentially\nAbundant\nMetabolites') + 
  geom_text_repel(force = 20,size = 6)
