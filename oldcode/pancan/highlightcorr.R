# A quick script to highlight correlations between metabolites

# Initial pancan metabolomics analysis
rm(list=ls())
library(pheatmap)
library(NMF)
library(RColorBrewer)
library(ggplot2)
library(reshape)
library(colorspace)
cl = colorRampPalette(c("blue","white","red"))(15)

temp = read.csv('/Users/ereznik/Documents/pancanmet/data/merged_metabolomics/alldata.csv', row.names = 1,nrows = 2)
splittypes = strsplit(colnames(temp),'\\.')
tissuetype = sapply(splittypes, "[", 1)
alldata = read.csv('/Users/ereznik/Documents/pancanmet/data/merged_metabolomics/alldata.csv', row.names = 1, header = F,skip = 3)
colnames(alldata) = colnames(temp)

########################################################################
# Look for common patterns in metabolite of your choice
hotmet = 'Citrate'
midx = which(rownames(alldata) == hotmet)

s2use = c('BLCA','BRCA','BRCATang','COAD','KICH','KIRC','OV','PRAD','PRADLODA','STAD')

p = matrix(1,dim(alldata)[1],length(s2use))
padj = matrix(1,dim(alldata)[1],length(s2use))
r = matrix(0,dim(alldata)[1],length(s2use))
colnames(p) = s2use;colnames(padj) = s2use;colnames(r) = s2use
rownames(p) = rownames(alldata);rownames(padj) = rownames(alldata); rownames(r) = rownames(alldata)

for (s in s2use){
    
    print(s) 
    
    whichdata = which(tissuetype == s)
    tempdata = alldata[,whichdata]
    whichmet = which(rownames(tempdata) == hotmet)
    if (length(which(is.na(as.numeric(tempdata[whichmet,]))))>0){next}
    
    
    for (i in 1:dim(tempdata)[1]){
        if (length(which(is.na(as.numeric(tempdata[i,]))))>0){next}
        temp = cor.test(as.numeric(tempdata[whichmet,]),as.numeric( tempdata[i,],method = 'spearman') )
        p[i,s] = temp$p.value
        r[i,s] = temp$estimate[[1]]
        if (is.nan(p[i,s])){next}
    }
    padj[,s] = p.adjust(p[,s],method = 'BH')
}

radj = r
radj[padj > 5e-2] = 0
msums = rowSums(abs(sign(radj)))

# Plot the results
aheatmap(radj[which(msums>2),],color = cl,breaks = 0)
