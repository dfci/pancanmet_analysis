# A quick script to visualize the differences in metabolism for "Warburg" vs Non-Warburg Tumors
setwd('/Users/ereznik/Documents/pancanmet/analysis/pancan/viz/')
# Initial pancan metabolomics analysis
rm(list=ls())
library(pheatmap)
library(NMF)
library(RColorBrewer)
library(ggplot2)
library(reshape)
library(colorspace)
library(data.table)

temp = fread('../../../data/merged_metabolomics/alldata.csv',sep = ',',header = F,skip = 3)
alldata = data.frame( temp, row.names = 1 )
title.line = strsplit( readLines('../../../data/merged_metabolomics/alldata.csv', n=1), ',' ) [[1]] # Note that first element is empty ""
colnames(alldata) = title.line[2:length(title.line)]
splittypes = strsplit(colnames(alldata),'\\.')
study = unlist( lapply(strsplit(colnames(alldata),':'), `[[`, 1) )

# Plot just metabolites of interest
hotmets = 1+c(194,13,  94, 101, 160, 165, 184,197, 201, 202, 209, 218, 241,
            254, 269, 291) # get this list from python
notna = which(!is.na(alldata[195,]))

hotdata =  alldata[hotmets,notna]
sortidx = order(hotdata[1,])
hotdata_sort = hotdata[,sortidx]
annCol = study[notna[sortidx]]
plotdata = log2(hotdata_sort)

aheatmap(plotdata,color = diverge_hsv(n=250,p=0.5), breaks = 0,
         Rowv = NA, Colv = NA,annCol = annCol)

########################################################################
# Look for common patterns in metabolite of your choice
hotmet = 'Lactate'
midx = which(rownames(alldata) == hotmet)

s2use = study

p = matrix(1,dim(alldata)[1],length(s2use))
fc = matrix(0,dim(alldata)[1],length(s2use))
colnames(p) = s2use;colnames(fc) = s2use
rownames(p) = rownames(alldata); rownames(fc) = rownames(alldata)

for (s in s2use){
  whichdata = which(study == s)
  tempdata = alldata[,whichdata]
  highlac = which(tempdata[midx,] > 1)
  lowlac = which(tempdata[midx,] < 1)
  if (length(lowlac)<3|length(highlac)<3){next}
  for (i in 1:dim(tempdata)[1]){
    if (length(which(is.na(as.numeric(tempdata[i,]))))>0){next}
    p[i,s] = wilcox.test(as.numeric(tempdata[i,highlac]),as.numeric( tempdata[i,lowlac]) )$p.value
    if (is.nan(p[i,s])){next}
    if (p[i,s]<0.05){fc[i,s] = mean(as.numeric(tempdata[i,highlac])) - mean(as.numeric( tempdata[i,lowlac]) )}
    
  }
}
signfc = sign(fc)
msums = rowSums(abs(signfc))

# Plot the results
notna = which(!is.na(alldata[midx,]))
plotdata = log2(alldata[which(msums>=3),notna])
sortidx = order(alldata[midx,notna])
annCol = study[notna]
aheatmap(plotdata,color = diverge_hsv(n=250,p=0.5), breaks = 0,
         Rowv = NA, Colv = sortidx,annCol = annCol)
