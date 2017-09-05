# Initial pancan metabolomics analysis
rm(list=ls())
library(pheatmap)
library(NMF)
library(RColorBrewer)
library(ggplot2)
library(reshape)
library(colorspace)
library(Hmisc)
cmap = brewer.pal(12,'Set3')
pal <- function(col, border = "light gray", ...){
  n <- length(col)
  plot(0, 0, type="n", xlim = c(0, 1), ylim = c(0, 1),
       axes = FALSE, xlab = "", ylab = "", ...)
  rect(0:(n-1)/n, 0, 1:n/n, 1, col = col, border = border)}


alldata = read.csv('/Users/ereznik/Documents/reporter/data/merged_metabolon/tnpaired_fc.csv', row.names = 1)
splittypes = strsplit(colnames(alldata),'\\.')
tissuetype = sapply(splittypes, "[", 1)

cormat = rcorr(t(alldata),type='spearman')
cormat$r[cormat$P<0.05] = 0
cormat$r[cormat$n<50] = 0
diag(cormat$r) = 0
cormat$r[which(is.na(cormat$r))] = 0
rmrows = which(rowSums(cormat$r) == 0)

pheatmap(cormat$r[-rmrows,-rmrows])
