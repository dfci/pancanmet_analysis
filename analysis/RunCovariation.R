# Script to analyze covariation of metabolites

rm(list=ls())
setwd('/Users/ereznik/pancanmet_analysis/analysis/')
source('functions/covariation.R')
source('functions/readBigMet.R')
library(ggplot2)

# Parameters to set
study2drop = c('BLCA')
pcrit = NA#0.05/2500
dropNA = TRUE
minsamples = 20

# Read in data
metdata = readBigMet( study2drop )
met = metdata$met; tissuetype = metdata$tissuetype; studytype = metdata$studytype
uqstudy = unique(studytype)

# Decide whether we will only use metabolites for which we have all measurements or not
if (dropNA){
  met = met[which(complete.cases(met)),]
}

# Initialize a vector to store results in
corrResults = matrix(NA,1,1)

# For each study, calculate the covariation of all metabolite pairs
for (study in uqstudy){
  for (tissue in c('Normal','Tumor')){
   
    studyix = which(studytype == study & tissuetype == tissue)
    if (length(studyix) < minsamples){
      next
    }
    savename = paste(study,':',tissue)
    returnlist = covariation( met[,studyix], pcrit )
    r = returnlist$r; p = returnlist$p; rlist = returnlist$rlist
    
    if (is.na( corrResults[1,1] )){
      corrResults = data.frame( rlist )
      colnames(corrResults) = savename
    }else{
      cnames = colnames(corrResults)
      corrResults = cbind( corrResults, rlist[rownames(corrResults)] )
      colnames(corrResults) = c(cnames,savename)
    }
     
  }
}

# Calculate standard deviation of the results
corrSD = apply(corrResults,1,sd)
numup = apply(corrResults,1,function(x){length(which(x>0))})
numdown = apply(corrResults,1,function(x){length(which(x<0))})
