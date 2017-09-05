# A script to analyze how different tumor metabolomic profiles are from normal metabolomic profiles.
# The approach is to compute the RMSD for all tumor-normal pairs, and compare to the RMSD for all normal-normal pairs.

rm(list = ls())
setwd('/Users/ereznik/Documents/pancanmet/analysis/pancan/')
source('../rfunc/multiPCAplot.R')
source('../rfunc/multiplot.R')
library(ggplot2)
library(data.table)
library(reshape)
library(NMF)

RMSD <- function(x1,x2) {
  rmsd = sqrt( sum( (x1 - x2)^2 )/ length(x1) )
}


# Some parameters for the code
pthresh = 5e-2
s2rm = 'LGG' # Studies to remove or skip
niter = 1e4 # Number of samples to draw in order to generate a ratio distribution
dim2plot = list(c(1,2),c(2,3),c(3,4),c(4,5))
uselog = FALSE # Log transform data for PCA and heatmap
num2check = 20 # Number of PCA dimensions to correlate
userank = FALSE

# Read in tissue types
ttypes = read.csv('/Users/ereznik/Documents/pancanmet/data/merged_metabolomics/merged_tissuetypes.csv',header = T,row.names = 1)
ttypes[,1] = as.character(ttypes[,1])

# Read in the metabolomics data
temp = fread('../../data/merged_metabolomics/alldata.csv',sep = ',',header = F,skip = 3)
alldata = data.frame( temp, row.names = 1 )
title.line = strsplit( readLines('../../data/merged_metabolomics/alldata.csv', n=1), ',' ) [[1]] # Note that first element is empty ""
colnames(alldata) = title.line[2:length(title.line)]

splittypes = strsplit(colnames(alldata),'\\:')
tissuetype = sapply(splittypes, "[", 1)
uqtis = unique(tissuetype)
uqtis = setdiff(uqtis, s2rm)
#uqtis = 'KIRC' # For testing!

# Make an array to store everything
res = data.frame( matrix(NA,niter,length(uqtis)) )
colnames(res) = uqtis

res2 = data.frame( matrix(NA,niter,length(uqtis)) )
colnames(res2) = uqtis

# For each study, calculate RMSE
for (tis in uqtis){
  ix = which(tissuetype == tis)
  d = alldata[,ix]
  d = d[which(complete.cases(d)),]
  
  # If we are using ranks only, modify here
  if (userank){
    d = t(apply(d,1,rank))
    dropix = which(apply(d,1,sd)==0)
    if (length(dropix)>0){
      d = d[-dropix,]
    }
    d = data.frame(d)
  }else{ # Just scale the data
    cnames = colnames(d)
    d = t(apply(d,1,scale))
    dropix = which(is.na(abs(rowSums(d))))
    if (length(dropix)>0){
      d = d[-dropix,]
    }
    d = data.frame(d)
    colnames(d) = cnames
  }
  
  tumidx = grep('Tumor',colnames(d))
  normidx = grep('Normal',colnames(d))
  
  # Do some PCA and plot
  if(uselog){
    d4pca = log2( d[,c(tumidx,normidx)] )
  }else{
    d4pca = d[,c(tumidx,normidx)]
  }
  
  pca = prcomp(d4pca, scale = TRUE)
  dtype = c(rep('Tumor',length(tumidx)),rep('Normal',length(normidx)))
  if (userank){
    multiPCAplot(pca,dim2plot,dtype,paste('../../results/tum2normvariation/PCA_',tis,'_RANK.pdf',sep=''))
  }else{
    multiPCAplot(pca,dim2plot,dtype,paste('../../results/tum2normvariation/PCA_',tis,'.pdf',sep=''))
  }
  
  
  # For each of the PCA components, calculate the correlation with tumor/normal labels
  pcacor = data.frame(matrix(0,num2check,2))
  for (i in 1:num2check){
    tempd = pca$rotation[,i]
    pcacor[i,1] = -log10( wilcox.test(tempd[which(dtype == 'Tumor')],tempd[which(dtype == 'Normal')])$p.value )
    pcacor[i,2] = round(pca$sdev[i]^2/sum(pca$sdev^2)*100,2)
  }
  cor2plot = data.frame( pcacor )
  cor2plot$Component = 1:num2check
  colnames(cor2plot) = c('P','Variance','Component')
  print( ggplot(cor2plot,aes(Component,P,label = Variance)) + geom_point() + theme_classic() + ggtitle(tis) + 
           geom_text() + xlab('Principal Component Index') + ylab('-Log10 P Value') + 
           ggsave(paste('../../results/tum2normvariation/CorrelationsPCA_',tis,'.pdf',sep='')) )
  
  # Also make a heatmap
  d4heatmap = as.matrix( scale(d4pca,center = TRUE, scale = TRUE) ) 
  aheatmap(d4heatmap,annCol = dtype,Rowv = FALSE,Colv = FALSE, width = 10,height = 10,
           filename = paste('../../results/tum2normvariation/Heatmap_',tis,'.pdf',sep=''))
  
#   # Make a data array to store everything
#   tum2tum = matrix(NA,length(tumidx),length(tumidx))
#   rownames(tum2tum) = colnames(d)[tumidx]
#   colnames(tum2tum) = colnames(d)[tumidx]
#   
#   for (i in 1:length(tumidx)){
#     for (j in 1:length(tumidx)){
#       
#       ix1 = tumidx[i]
#       ix2 = tumidx[j]
#       tum2tum[i,j] = RMSD(d[,ix1],d[,ix2])
#       
#     }
#   }
#   
#   # Make a data array to store everything
#   tum2norm = matrix(NA,length(tumidx),length(normidx))
#   rownames(tum2norm) = colnames(d)[tumidx]
#   colnames(tum2norm) = colnames(d)[normidx]
#   
#   for (i in 1:length(tumidx)){
#     for (j in 1:length(normidx)){
#       
#         ix1 = tumidx[i]
#         ix2 = normidx[j]
#         tum2norm[i,j] = RMSD(d[,ix1],d[,ix2])
#       
#     }
#   }
# 
#   # Make a data array to store everything
#   norm2norm = matrix(NA,length(normidx),length(normidx))
#   rownames(norm2norm) = colnames(d)[normidx]
#   colnames(norm2norm) = colnames(d)[normidx]
#   
#   for (i in 1:length(normidx)){
#     for (j in 1:length(normidx)){
#       
#       ix1 = normidx[i]
#       ix2 = normidx[j]
#       norm2norm[i,j] = RMSD(d[,ix1],d[,ix2])
#       
#     }
#   }
#   
#   # Flatten and merge the data
#   n2nflat = data.frame(as.vector(norm2norm))
#   n2nflat$Type = 'Normal2Normal'
#   colnames(n2nflat) = c('Value','Type')
#   
#   t2nflat = data.frame(as.vector(tum2norm))
#   t2nflat$Type = 'Tumor2Normal'
#   colnames(t2nflat) = c('Value','Type')
#   
#   t2tflat = data.frame(as.vector(tum2tum))
#   t2tflat$Type = 'Tumor2Tumor'
#   colnames(t2tflat) = c('Value','Type')
#   
#   merge = rbind(t2nflat,n2nflat,t2tflat)
#   
#   # Plot the distributions
#   print( ggplot(merge,aes(Type,Value)) + geom_violin() + geom_point() + ggtitle(tis) )
#   
#   # Plot a ratio distribution
#   t2n_sample = sample(t2nflat[,1],niter,replace = TRUE)
#   n2n_sample = sample(n2nflat[,1],niter,replace = TRUE)
#   t2t_sample = sample(t2tflat[,1],niter,replace = TRUE)
#   
#   # Save t2t and t2n comparison to n2n
#   res[,tis] = t2t_sample/n2n_sample
#   res2[,tis] = t2n_sample/n2n_sample
  
  
}
# resmelt = melt(res)
# resmelt[,1] = as.character(resmelt[,1])
# resmelt$Tissue = ttypes[resmelt[,1],1]
# 
# res2melt = melt(res2)
# res2melt[,1] = as.character(res2melt[,1])
# res2melt$Tissue = ttypes[res2melt[,1],1]
# 
# # Plot
# p1 = ggplot(resmelt,aes(variable,log2(value),fill = Tissue)) + geom_boxplot() + 
#   theme_classic() + geom_hline() + ggtitle('T2T/N:N') +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# p2 = ggplot(res2melt,aes(variable,log2(value),fill = Tissue)) + geom_boxplot() + 
#   theme_classic() + geom_hline() + ggtitle('T2N/N:N') +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# if (userank){
#   pdf('../../results/tum2normvariation/ratios_RANK.pdf',height = 10,width = 10)
#   multiplot(p1,p2)
#   dev.off()
# }else{
#   pdf('../../results/tum2normvariation/ratios.pdf',height = 10,width = 10)
#   multiplot(p1,p2)
#   dev.off()
# }
