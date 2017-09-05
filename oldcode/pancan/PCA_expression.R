# A script to make a PCA plot of metabolic expression data for the pancanmet project.

rm(list = ls())
setwd('/Users/ereznik/Documents/pancanmet/analysis/pancan/')

library(data.table)
library(ggplot2)
library(limma)
source('../rfunc/multiPCAplot.R')

# Set parameters for run
dim2plot = list(c(1,2),c(2,3))

# Set directories
rnadir = '/Users/ereznik/Documents/useful/rnaseq_unnormdata/rnaseq_nov42014_RSEM_normalized/'
savedir = '/Users/ereznik/Documents/pancanmet/results/tum2normvariation/'

# Read in the metabolite gene data
mgenes = read.csv('../../data/Recon2v04_genenames.csv',header = T,row.names = 1)

# Set the studies to read
studies = c('BLCA','BRCA','KIRC','PAAD','PRAD')
#studies = c('BLCA') # For testing

# For each study
for (s in studies){
  print(s)
  
  # Read in the RNA
  temp = fread(paste(rnadir,s,'.txt',sep=''),sep = '\t',header = F,skip = 2)
  rna = data.frame( temp, row.names = 1 )
  title.line = strsplit( readLines(paste(rnadir,s,'.txt',sep=''), n=1), '\t' ) [[1]] # Note that first element is empty ""
  colnames(rna) = title.line[2:length(title.line)]
  
  # Rename the rows
  rnames = unlist( lapply(rownames(rna),function(x){y = strsplit(x,'\\|')[[1]];z=y[2]}))
  rnames = as.character( rnames ) # just in case
  rownames(rna) = rnames
  
  # Get the indices of primary cell lines and metastatic cell lines
  splittypes = strsplit(colnames(rna),'\\-')
  sampletype = sapply(splittypes, "[", 4)
  
  primaryidx=c( grep("01+",sampletype), grep("03+",sampletype) )
  metaidx=c( grep("06+",sampletype), grep("07+",sampletype) )
  normidx=c( grep("11",sampletype) )
  
  # Get rid of genes which have very low read counts (below 16, on average)
  rna = rna[which(apply(rna,1,mean)>16),]
  
  # Restrict ourselves to just tumor and normal cells
  d = rna[,c(primaryidx,normidx)]
  dtype = c( rep('Tumor',length(primaryidx)), rep('Normal',length(normidx)) )
 
  # Get the expected counts, run PCA
  d2plot = log2(d+1)
  pca = prcomp( d2plot, scale. = TRUE)
  multiPCAplot(pca,dim2plot,paste(savedir,dtype,'/PCA_',s,'__AllGenes.pdf',sep = '') ) 
  
  # Repeat, now just using metabolic genes
  dmet = log2( d[unique(mgenes[,1]),] + 1)
  dim1 = dim(dmet)[1]
  
  # Get rid of NAs
  dmet = dmet[which(complete.cases(dmet)),]
  dim2 = dim(dmet)[1]
  
  print(paste('Number of metabolic genes not found:', dim1-dim2) )
  
  # Do PCA again
  pcamet = prcomp( dmet, scale. = TRUE)
  multiPCAplot(pcamet,dim2plot,paste(savedir,dtype,'/PCA_',s,'__MetabolicGenes.pdf',sep = '') ) 
  
}