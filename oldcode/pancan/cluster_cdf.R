# Script to calculate Jaccard distance and do PCA plot from FC data

rm(list=ls())
setwd('/Users/ereznik/Documents/pancanmet/analysis/pancan/')
library(pheatmap)
library(NMF)
library(RColorBrewer)
library(ggplot2)
library(reshape)
library(colorspace)

calcdist = 0
maxval = 3

br = seq( -maxval, maxval,.1)
brb = length(which(br<0)); bra = length(which(br>0))
cl = colorRampPalette(c("blue","white","red"))(length(br)-1)

alldata = read.csv('/Users/ereznik/Documents/pancanmet/data/merged_metabolomics/tnpaired_fc.csv', row.names = 1)
splittypes = strsplit(colnames(alldata),'\\.')
tissuetype = sapply(splittypes, "[", 1)
uqtis = unique(tissuetype)

# For each metabolite, identify the number of studies which it is profiled in
numtis = vector()
for (i in 1:dim(alldata)[1]){
  numtis[i] = length( unique( tissuetype[which(!is.na( alldata[i,] ) ) ] ) )
}

if (calcdist == 1){
  
  # Calculate a Jaccard distance between each sample
  bindata = log2( alldata )
  cdfdata = matrix(0,dim(bindata)[1],dim(bindata)[2])
  for (i in 1:(dim(bindata)[1])){
    notna = which(!is.na(bindata[i,]))
    tempcdf = ecdf( as.numeric( bindata[i,notna] ) )
    cdfdata[i,] = tempcdf( bindata[i,] )
    print(i)
  }
  
  J = matrix(0,dim(cdfdata)[2],dim(cdfdata)[2])
  
  for (i in 1:(dim(cdfdata)[2]-1)){
    print(i)
    for (j in (i+1):dim(cdfdata)[2]){
      
      naidx = which( !is.na(cdfdata[,i]) & !is.na(cdfdata[,j]) )
      d = sum( (cdfdata[naidx,i] - cdfdata[naidx,j])^2 )/sqrt(length(naidx))
      
      J[i,j] = d
      J[j,i] = d
      
    }
  }
  save(J,file = '../../data/clustering/Jdist_cdf.Rdata')
} else(load('../../data/clustering/Jdist_cdf.Rdata'))

# Convert to distance object
d = dist(J)

# Calculate distances
alldata[which(is.na(alldata),arr.ind = T)] = 1

# # Scale and make pretty
# plotdata = t(scale(t(alldata)))

# Use fold change
plotdata = log2(alldata)
plotdata[which(plotdata > maxval,arr.ind = T)] = maxval
plotdata[which(plotdata < -maxval,arr.ind = T)] = -maxval

# Make annotation
annCol = data.frame( tissuetype )
rownames(annCol) = colnames(plotdata)
colnames(annCol) = 'Study'
plotdata[which(is.na(plotdata),arr.ind = T)] = 0
plotdata_final = plotdata[which(numtis>0),]


# Plot
clust = pheatmap(plotdata,clustering_distance_cols = d,show_rownames = F, show_colnames = F,
                 breaks = br, color = cl, annotation =  annCol,
                 filename = '../../results/clustering_cdf.pdf',width = 10,height = 10)

# Separate into two groups of samples, find major differences
s = cutree(clust$tree_col,2)

# Create a matrix to save the differential abundance results from below
diffdata = matrix(0,dim(alldata)[1],length(uqtis))
colnames(diffdata) = uqtis
rownames(diffdata) = rownames(alldata)

for (tistype in uqtis){
  print(tistype)
  c1 = which( tissuetype == tistype & s == 1)
  c1names = colnames(alldata)[c1]
  c2 = which( tissuetype == tistype & s == 2)
  c2names = colnames(alldata)[c2]
  if (length(c1) < 2 | length(c2)< 2){next}
  allidx = c(c1,c2)
  
  # Calculate the differences in FC between these two sample groups
  tempdata = alldata[,allidx]
  diff = vector(); 
  p = vector(); 
  for (i in 1:dim(tempdata)[1]){
    std = sd(tempdata[i,])
    if (std == 0){diff[i] = 0;p[i] = 1;next}
    res = t.test( tempdata[i,c1names],tempdata[i,c2names])
    diff[i] = log2( mean( t( tempdata[i,c1names] ) ) / mean( t( tempdata[i,c2names] ) ) )
    p[i] = res$p.value
  }
  names(diff) = rownames(tempdata)
  names(p) = rownames(tempdata)
  padj = p
  padj = p.adjust(padj,method = 'BH')
  qplot(diff,-log10(padj))
  
  diffdata[,tistype] = diff
  diffdata[which(padj > 1e-2),tistype] = 0
}

plotdata = diffdata
plotdata = plotdata[-which(rowSums(sign(abs(plotdata))) < 2),] # Keep only metabolites changing repeatedly in X number of studies
plotdata[plotdata > maxval] = maxval
plotdata[plotdata < -maxval] = -maxval
pheatmap(plotdata, breaks = br, color = cl,)

