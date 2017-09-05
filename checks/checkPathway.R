# Script to check pathway calculations

# Result: there are some very minor differences. They probably result from one or two meabolites not being counted the same way between the scripts.
# Qualitatively and quantitatively, the results are nearly identical.

rm(list = ls())
setwd('/Users/ereznik/pancanmet_analysis/analysis/')
library(ggplot2)
library(reshape)
source('plottingconventions.R')

# Read in the pathway map
map = read.csv('../data/KEGG_pathwaymap.csv',header = TRUE,row.names = 1,check.names = FALSE)
map = as.matrix(map)

# Read in the prior results
pathway = read.csv('../results/pathway/AllPathway.csv',header = TRUE,row.names = 1)
dacalc = read.csv('../results/diffabundance/DifferentialAbundanceSummary.csv',header = TRUE,row.names = 1)

# Check each study and pathway
numerror = 0
numcorrect = 0
for (study in colnames(pathway)){
  
  da_study = dacalc[,c(paste0('FC.',study),paste0('Padj.',study))]
  
  for (pway in rownames(pathway)){
    
    # Get the metabolites in the pathway
    mets = colnames(map)[which(map[pway,]==1)]
    temp = dacalc[mets,c(paste0('FC.',study),paste0('Padj.',study))]
    temp2 = temp[complete.cases(temp),]
    
    if (dim(temp2)[1] == 0 & is.na(pathway[pway,study])){
      numcorrect = numcorrect + 1
      next
    }
    
    if (dim(temp2)[1] == 0 & !is.na(pathway[pway,study])){
      print('Error with NAs!')
      print(paste(pway,study,round(pathway[pway,study],2)))
      numerror = numerror + 1
      next
    }
    
    numup = which(temp2[,1] > 0 & temp2[,2] < pthresh)
    numdown = which(temp2[,1] < 0 & temp2[,2] < pthresh)
    psize = dim(temp2)[1]
    
    da = (length(numup) - length(numdown))/psize
    diffda = abs(da - pathway[pway,study])
    if (diffda > 1e-10){
      print('Error!')
      print(paste(study,pway,round(da,2),round(pathway[pway,study],2),'Diff=',diffda))
      numerror = numerror + 1
      #print(temp2)
      #stop()
    }else{
      numcorrect = numcorrect + 1
    }
    # }else{
    #   print(paste(pway,study,da))
    # }
    
  }
}

print(paste0('Numcorrect = ',numcorrect))
print(paste0('Numerror = ',numerror))
