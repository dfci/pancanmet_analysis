# Initial pancan metabolomics analysis
rm(list=ls())
library(pheatmap)
library(NMF)
library(RColorBrewer)
library(ggplot2)
library(reshape)
library(colorspace)
cmap = brewer.pal(12,'Set3')
pal <- function(col, border = "light gray", ...){
  n <- length(col)
  plot(0, 0, type="n", xlim = c(0, 1), ylim = c(0, 1),
       axes = FALSE, xlab = "", ylab = "", ...)
  rect(0:(n-1)/n, 0, 1:n/n, 1, col = col, border = border)}


alldata = read.csv('/Users/ereznik/Documents/reporter/data/merged_metabolon/tnpaired_fc.csv', row.names = 1)
splittypes = strsplit(colnames(alldata),'\\.')
tissuetype = sapply(splittypes, "[", 1)

# Trim data to include no NA's
nonena = which( apply(alldata, 1,function(x)sum(is.na(x))) == 0 )
trimdata = alldata[nonena,]
logtrimdata = log2(trimdata)
anncol = data.frame('TissueType' = tissuetype); rownames(anncol) = colnames(logtrimdata)

# Also set colors
anncolors = list()
ctr = 0
for (i in 1:dim(anncol)[2] ){
  anncolors[[i]] = cmap[1:length(unique(anncol[,i]))]
  ctr = ctr + 1
}
col_breaks = c(seq(-min(logtrimdata),-1,length=100), # for red
               seq(-.99,1,length=100), # for yellow
               seq(1.01,max(logtrimdata),length=100)) # for green
aheatmap(as.matrix(logtrimdata), color = diverge_hcl(50,p=1), breaks = 0,annCol = anncol)

pca = prcomp(logtrimdata)
plotdata = data.frame( pc1 = pca$rotation[,1], pc2 = pca$rotation[,2], tissue = tissuetype )
ggplot(plotdata,aes(pc1,pc2)) + geom_point(aes(color=tissue, size = 5))

# Repeat process, but now set NA's to 1
morena = which( apply(alldata, 1,function(x)sum(is.na(x))) <= 100  )
moredata = alldata[morena,]
logmoredata = log2(moredata)
logmoredata[is.na(logmoredata)] = 0
aheatmap(as.matrix(logmoredata), color = diverge_hcl(50,p=2), breaks = 0,annCol = anncol)

# Focus in on a metabolite and plot it
hotmet = 'Lactate'
hotdata = t(alldata[hotmet,])
idx = as.vector( order(hotdata) )
plotdata = data.frame(1:length(idx),log2(hotdata[idx]),tissuetype[idx]);colnames(plotdata)=c('x','y','tissue')
ggplot(plotdata,aes(x,y,fill = tissue)) + geom_bar(stat="identity") + theme_classic()


