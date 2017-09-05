# Script to examine covariation of metabolites in select samples

# MAKE SURE TO CHECK THAT THE FISHER P-VALUE CALCULATED BELOW IS CORRECT
rm(list = ls())
library(NMF)
library(RColorBrewer)
library(ggplot2)
library(reshape)
library(colorspace)
library(data.table)
library(pheatmap)
cmap = brewer.pal(12,'Set3')
setwd('/Users/ereznik/Documents/pancanmet/analysis/warburg/')
pthresh = 0.05
ttype = 'Tumor' # 'Tumor' or 'Normal'
hotmet = 'Lactate'

selectfirst <- function(s) strsplit(s, ":")[[1]][1]
selectlast <- function(s) strsplit(s, ":")[[1]][ length( strsplit(s, ":")[[1]] )]

pal <- function(col, border = "light gray", ...){
  n <- length(col)
  plot(0, 0, type="n", xlim = c(0, 1), ylim = c(0, 1),
       axes = FALSE, xlab = "", ylab = "", ...)
  rect(0:(n-1)/n, 0, 1:n/n, 1, col = col, border = border)}

# Read in all the data
temp = fread('../../data/merged_metabolomics/alldata.csv',sep = ',',header = F,skip = 3)
alldata = data.frame( temp, row.names = 1 )
title.line = strsplit( readLines('../../data/merged_metabolomics/alldata.csv', n=1), ',' ) [[1]] # Note that first element is empty ""
colnames(alldata) = title.line[2:length(title.line)]

# Drop empty rowname
emptyix = which(rownames(alldata)=='')
if (length(emptyix) !=0){
  alldata = alldata[-emptyix,]
}

# Extract relevant data from columns
study = unlist( lapply(strsplit(colnames(alldata),':'), `[[`, 1) )
uqtis = unique(study)

################################################################################################################

# For each study, compare the low lactate patients to the high lactate patients
corr = data.frame( matrix(NA,dim(alldata)[1],length(uqtis)) )
rownames(corr) = rownames(alldata)
colnames(corr) = uqtis
corp = corr
corpadj = corr
for (tis in uqtis){
  print(tis)
  
  # Get the samples of this tissue type
  regstr = paste(tis,':.*',ttype,sep = '')
  ix = grep(regstr,colnames(alldata))
  if (length(ix) < 5){
    tisix = which(colnames(corr) == tis)
    corr = corr[,-tisix] 
    corp = corp[,-tisix]
    corpadj = corpadj[,-tisix]
    next
  }
  d = alldata[,ix]
  d = d[which(complete.cases(d)),]
  d = as.matrix(d)
  
  if (!(hotmet %in% rownames(d))){
    tisix = which(colnames(corr) == tis)
    corr = corr[,-tisix] 
    corp = corp[,-tisix]
    corpadj = corpadj[,-tisix]
    next
  }
  
  # Correlate everything with lactate levels
  tempp = data.frame( matrix(NA, dim(d)[1],1) )
  rownames(tempp) = rownames(d)
  for (i in 1:dim(d)[1]){
    mname = rownames(d)[i]
    tempcor = cor.test(d[i,],d[hotmet,])
    corr[mname,tis] = tempcor$estimate
    tempp[mname,1] = tempcor$p.value
  }
  
  # Now we need to adjust p-values
  temppadj = data.frame( p.adjust(tempp[,1],method = 'BH') )
  rownames(temppadj) = rownames(tempp)
 
  # Write in values
  corp[rownames(tempp),tis] = tempp
  corpadj[rownames(temppadj),tis] = temppadj

}

# Now correct p-values equal to zero
corp[which(corp < .Machine$double.eps,arr.ind=T)] = .Machine$double.eps
corpadj[which(corpadj < .Machine$double.eps,arr.ind=T)] = .Machine$double.eps


# Make an aggregated p-value for Fisher's method using uncorrected p-values
fisherp = data.frame( matrix(NA,dim(alldata)[1],8) )
rownames(fisherp) = rownames(alldata)
colnames(fisherp) = c('TestStatistic','P-value','AdjustedP-Value','NumStudies','NumPositiveCor','NumNegativeCor','PctPositiveCor','PctNegativeCor')
for (i in 1:dim(alldata)){
  notNA = which(!is.na(corp[i,]))
  if (length(notNA)==0){next}
  
  fisherp[i,1] = -2*sum(log(corp[i,notNA]))
  fisherp[i,2] = pchisq( fisherp[i,1], 2*length(notNA), lower.tail=FALSE) #NOT SURE IF THIS IS CORRECT!!!!
  
  fisherp[i,4] = length(notNA)
  fisherp[i,5] = length( which(corp[i,notNA]<pthresh & corr[i,notNA] > 0) )
  fisherp[i,6] = length( which(corp[i,notNA]<pthresh & corr[i,notNA] < 0) )
  
  if (fisherp[i,4]>8){
    fisherp[i,7] = fisherp[i,5] / fisherp[i,4]
    fisherp[i,8] = fisherp[i,6] / fisherp[i,4]
  }else{
    fisherp[i,7] = ''
    fisherp[i,8] = ''
  }

}

# Correct p-values
fisherp[,3] = p.adjust(fisherp[,2],method = 'BH')

# Write results to file
writedata = cbind(corp,fisherp)

# Drop all rows with only NA values
NArows = which( apply(fisherp, 1, function(x) all(is.na(x))) )
writedata = writedata[-NArows,]

write.csv(writedata,paste('../../results/metabolitecorrelations/',hotmet,'_',ttype,'_Correlations.csv',sep=''))

# Make a small plot which indicates correlations with glycolytic metabolites
glycmet = read.csv('GlycolysisMetabolites.csv',header = F)
glycmet[,1] = as.character(glycmet[,1])
glycor = corr[glycmet[,1],]
glyp = corp[glycmet[,1],]
glycor[which(glyp > pthresh, arr.ind = T)] = 0
glyplot = melt(as.matrix(glycor))
colnames(glyplot) = c('Metabolite','Study','Value')
ggplot(glyplot,aes(Study,Metabolite,color = Value,size = Value)) + geom_point() + theme_bw() + 
  scale_color_gradient2(low = 'blue',mid = 'white',high = 'red')
