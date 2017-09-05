# Initial pancan metabolomics analysis
rm(list=ls())
library(pheatmap)
library(NMF)
library(RColorBrewer)
library(ggplot2)
library(reshape)
cmap = brewer.pal(12,'Set3')

alldata = read.csv('/Users/ereznik/Documents/reporter/data/merged_metabolon/alldata.csv', row.names = 1)
metdata = data.frame( sapply(alldata[3:dim(alldata)[1],1:(dim(alldata)[2]-2)], function(x) { as.numeric(as.character(x)) } ), check.names = F )
rownames(metdata) = rownames(alldata)[3:dim(alldata)[1]]
mettype = alldata[3:dim(alldata)[1],dim(alldata)[2]-1]
submettype = alldata[3:dim(alldata)[1],dim(alldata)[2]]
logdata = log2( metdata )
tumortype = alldata[1,1:(dim(alldata)[2]-2)]; names(tumortype) = colnames(metdata)[1:(dim(alldata)[2]-2)]; 
tissuetype = alldata[2,1:(dim(alldata)[2]-2)]; names(tumortype) = colnames(metdata)[1:(dim(alldata)[2]-2)]; 
tissuetype = factor(as.matrix(tissuetype),levels= c('Normal','Tumor','Met'))

# Trim data to include no NA's
nonena = which( apply(logdata, 1,function(x)sum(is.na(x))) == 0 )
trimlogdata = logdata[nonena,]

nonena = which( apply(metdata, 1,function(x)sum(is.na(x))) == 0 )
trimmetdata = metdata[nonena,]

anncol = data.frame('Study' = t(tumortype), 'TissueType' = t(tissuetype))
annrow = data.frame( mettype = mettype[nonena])

# Also set colors
anncolors = list()
ctr = 0
for (i in 1:(dim(anncol)[2]-1)){
  anncolors[[i]] = cmap[1:length(unique(anncol[,i]))]
  ctr = ctr + 1
}
anncolors[[ctr+1]] = c('white','black','red')
ctr = ctr + 1
for (i in 1:dim(annrow)[2]){
  anncolors[[i+ctr]] = cmap[1:length(unique(annrow[,i]))]
}
aheatmap(as.matrix(trimlogdata),annCol = anncol, annRow = annrow, Rowv=FALSE, Colv = FALSE, fontsize = 10, 
         color = '-RdBu:100', cexRow = 0.8, annColors = anncolors)

# # Missing values
# df_merge = melt(as.matrix(logdata)); colnames(df_merge) = c('x','y','z')
# ggplot(df_merge, aes(x = x, y = y)) +
#   geom_tile(data = subset(df_merge, !is.na(z)), aes(fill = z)) +
#   geom_tile(data = subset(df_merge,  is.na(z)), aes(colour = NA),
#             linetype = 0, fill = "pink", alpha = 0.5)
