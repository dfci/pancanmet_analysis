# A script to look at paired correlations in tumor/normal tissues in cross-cancer metabolomics data

rm(list = ls())
setwd('/Users/ereznik/Documents/pancanmet/analysis/pancan/')

library(ggplot2)
library(data.table)
library(pspearman)

thresh = 5e-2

temp = fread('../../data/merged_metabolomics/alldata.csv',sep = ',',header = F,skip = 3)
alldata = data.frame( temp, row.names = 1 )
title.line = strsplit( readLines('../../data/merged_metabolomics/alldata.csv', n=1), ',' ) [[1]] # Note that first element is empty ""
colnames(alldata) = title.line[2:length(title.line)]

# Extract relevant data from columns
study = unlist( lapply(strsplit(colnames(alldata),':'), `[[`, 1) )

# Read in tumor/normal pairs
tnpairs = read.csv('../../data/merged_metabolomics/tumor_normal_pairs4R.csv',header = F,colClasses=c("character","character"))
tnpairs_studies = unlist( lapply( tnpairs[,1],function(x){strsplit(x,'\\:')[[1]][1]}) )

# Initialize a data frame for saving results
uqstudy = unique(study)
savedata = matrix(0,dim(alldata)[1],length(uqstudy))
rownames(savedata) = rownames(alldata)
colnames(savedata) = uqstudy

# For each study, grab the tumor/normal pairs, correlate and plot
for (s in uqstudy){
  print(s)
  
  ix = which(study == s)
  d = alldata[,ix]
  d = d[complete.cases(d),]
  
  # Find the rows which have this study in them
  sidx = which(tnpairs_studies == s)
  ixtum = tnpairs[sidx,1]
  ixnorm = tnpairs[sidx,2]
  
  if (length(ixnorm)==0){next}
  
  dtum = d[,ixtum]
  dnorm = d[,ixnorm]
  
  res = data.frame( matrix(1,dim(d)[1],3) )
  rownames(res) = rownames(dtum)
  colnames(res) = c('R','P','Padj')
  res[,'R'] = 0
  
  for (i in 1:dim(res)[1]){
    
    if (sd(dtum[i,])<1e-10 | sd(dnorm[i,])<1e-10){
      next
    }
    #temp = cor.test(t(dtum[i,]),t(dnorm[i,]),method='spearman')
    temp = spearman.test(t(dtum[i,]),t(dnorm[i,]),approximation = 't-distribution')
    
    # If we have an especially low p-value, confirm it
    if (temp$p.value < 2e-16){
      temp2 = spearman.test(t(dtum[i,]),t(dnorm[i,]))
      if (abs(temp$p.value - temp2$p.value)>1e-16){
          print('Error in computing p values!')
      }
    }
    res[i,1:2] = c(temp$estimate,temp$p.value )
    
  }
  
  # Get rid of zero p-values
  #res[which(res[,2] <= .Machine$double.eps),2] = .Machine$double.eps
  
  res[,3] = p.adjust(res[,2])
  res[,3] = -log10(res[,3])
  
  plotdata = res
  plotdata$Color = 'Not Significant'
  plotdata[which( res[,'Padj']> -log10(thresh)),'Color'] = 'Significant'
  colnames(plotdata) = c('R','P','PAdj','Color')
  
  print( ggplot(plotdata,aes(R,PAdj,color = Color)) + geom_point() + ggtitle(s) + theme_minimal())
  
  # Save
  upidx = which(res[,1] > 0 & res[,3] > -log10(thresh))
  downidx = which(res[,1] < 0 & res[,3] > -log10(thresh))
  
  savedata[upidx,s] = 1
  savedata[downidx,s] = -1
}

rsums = rowSums(savedata)
print(head(sort(rsums, decreasing = T),30))
