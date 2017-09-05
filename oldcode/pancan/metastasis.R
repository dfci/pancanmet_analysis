# Script to analyze metastases in reporter data and identify common features

rm(list=ls())
setwd('/Users/ereznik/Documents/pancanmet/analysis/pancan/')
library(ggplot2)

# Set significance threshold
thresh = 0.05

# Read in all data, isolate metastases
alldata = read.csv('/Users/ereznik/Documents/pancanmet/data/merged_metabolomics/alldata.csv', row.names = 1,check.names = F)
m = data.frame( sapply(alldata[3:dim(alldata)[1],], function(x) { as.numeric(as.character(x)) } ), check.names = F )
rownames(m) = rownames(alldata)[3:dim(alldata)[1]]
splittypes = strsplit(colnames(alldata),'\\:')
tisstype = sapply(splittypes, "[", 3)

cancertype = alldata[1,1:(dim(alldata)[2])]
cancertype = sapply(unlist(cancertype),as.character)

cancer2include = unique( cancertype[which(tisstype == 'Metastasis')] ) 

# Trim the data down to just these cancer types, drop NA's
idx2keep = which(cancertype%in%cancer2include)
m2 = m[,idx2keep]
m2 = m2[which(complete.cases(m2)),]
m2 = as.matrix(m2)
ttype = tisstype[idx2keep]
ctype = cancertype[idx2keep]

# For each cancer type with metastsases, do a U test
fc = matrix(0,dim(m2)[1],length(cancer2include))
p = matrix(1,dim(m2)[1],length(cancer2include))

ctr = 1
for (c in cancer2include){
  metidx = which(ctype == c & ttype == 'Metastasis')
  otheridx = which(ctype == c & ttype != 'Metastasis')
  
  for (m in 1:dim(m2)[1]){
    p[m,ctr] = wilcox.test( m2[m,metidx],m2[m,otheridx])$p.value
    fc[m,ctr] = log2( mean(m2[m,metidx])/mean(m2[m,otheridx]) )
    
  }
  ctr = ctr + 1
}

# Adjust things
padj = p
for (i in 1:dim(p)[2]){
  padj[,i] = p.adjust(p[,i],method = 'BH')
}
fcadj = fc

fcadj[which(padj>thresh,arr.ind = T)] = 0
qplot(fcadj[,1],fcadj[,2])

# Print metabolites which change in both
print(rownames(m2)[which(fcadj[,1]>0 & fcadj[,2]>0)])
print(rownames(m2)[which(fcadj[,1]<0 & fcadj[,2]<0)])
