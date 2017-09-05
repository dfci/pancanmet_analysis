# Script to do clustering using only complete data

rm(list=ls())
setwd('/Users/ereznik/Documents/pancanmet/analysis/pancan/')
library(pheatmap)
library(NMF)
library(RColorBrewer)
library(ggplot2)
library(reshape)
library(colorspace)

maxval = 3

br = seq( -maxval, maxval,.1)
brb = length(which(br<0)); bra = length(which(br>0))
cl = colorRampPalette(c("blue","white","red"))(length(br)-1)

alldata = read.csv('/Users/ereznik/Documents/pancanmet/data/merged_metabolomics/tnpaired_fc.csv', row.names = 1)
splittypes = strsplit(colnames(alldata),'\\.')
tissuetype = sapply(splittypes, "[", 1)
uqtis = unique(tissuetype)

fulldata = log2( alldata[which(complete.cases(alldata)),] )
fulldata[which(fulldata>maxval,arr.ind=T)] = maxval
fulldata[which(fulldata< -maxval,arr.ind=T)] = -maxval

# Make annotation
annCol = data.frame( tissuetype )
rownames(annCol) = colnames(fulldata)
colnames(annCol) = 'Study'

pheatmap(fulldata, show_colnames = F,
         breaks = br, color = cl, annotation =  annCol)

