# Script to compare differential abundance results to those from analysis of individual study data. 

rm(list = ls())
setwd('/Users/ereznik/pancanmet_analysis/checks/')
library(ggplot2)
library(cowplot)

# Read in big table of results
bigres = read.csv('../results/diffabundance/DifferentialAbundanceSummary.csv',header = TRUE,row.names = 1)

# Read in merged mapping
mapping = read.csv('../import/tempdir/merged_mapping.csv',header = TRUE,row.names = 1, stringsAsFactors = FALSE)

# Read in big data file too, to get names for merged mapping
bigmet = read.csv('../data/merged_metabolomics/merged_metabolomics.csv', row.names = 1,check.names = FALSE,header = TRUE)
rownames(mapping) = rownames(bigmet)

calcda = function(tumor,normal){
  # Given two collections of numerical data, calculate differential abundance using variety of tests
  wp = wilcox.test(tumor,normal)$p.value
  
  if (sd(tumor) == 0 | sd(normal) == 0){
    tp = NA
    tfc = NA
  }else{
    ttest = t.test(tumor,normal)
    tp = ttest$p.value
    tfc = log2( ttest$estimate[1] / ttest$estimate[2] )
  }
  
  return(c(tfc,wp,tp))
}

ds = list.dirs('../data/metabolomics/',recursive = FALSE)
for (d in ds){
  if (grepl('oldstudies',d)){
    next
  }
  
  print(d)
  study = strsplit(d,'\\/')[[1]][5]
  met = read.csv(paste(d,'/',study,'_metdata.csv',sep = ''),header = TRUE,row.names = 1)
  mettrim = met[2:dim(met)[1],]
  mettrim = apply(mettrim,2,as.numeric)
  tissueclasses = as.character( t(met[1,]) ) #need to transpose otherwise it throws weird error
  
  tumix = which(tissueclasses %in% c('Tumor','Metastasis'))
  normix = which(tissueclasses %in% 'Normal')
  
  if (length(normix) == 0){
    # Skip, no normals
    next
  }
  
  metres = data.frame( 'fc' = numeric(0), 'wp' = numeric(0), 'tp' = numeric(0) )
  for (i in 1:dim(mettrim)[1]){
    tumtemp = mettrim[i,tumix]
    normtemp = mettrim[i,normix]
    
    metres[i,] = calcda(tumtemp,normtemp )
    
  }
  
  metres$wpadj = p.adjust( metres$wp, method = 'BH')
  metres$tadj = p.adjust( metres$tp, method = 'BH')
  rownames(metres) = sapply(rownames(met)[2:dim(met)[1]], tolower)
  
  # Now, map the results to the common names
  commonmap = data.frame( mapping[,study] )
  commonmap$commonname = rownames(mapping)
  commonmap = commonmap[ -which(is.na(commonmap)),]
  rownames(commonmap) = commonmap[,1]
  
  metres$commonname = commonmap[rownames(metres),'commonname']
  
  # Make plots
  metres$mergefc = bigres[metres$commonname, paste('FC.',study,sep = '')]
  metres$mergep = bigres[metres$commonname, paste('Padj.',study,sep = '')]
  p1 = ggplot(metres,aes(fc,mergefc)) + geom_point() + theme_bw() + 
    xlab('Check') + ylab('Original') + ggtitle(paste(study,'FC')) + geom_abline()
  p2 = ggplot(metres,aes(-log10(wpadj),-log10(mergep))) + geom_point() + 
    theme_bw() + xlab('Check') + ylab('Original') + ggtitle(paste(study,'WilcoxP')) + geom_abline()
  p3 = ggplot(metres,aes(-log10(tadj),-log10(mergep))) + geom_point() + 
    theme_bw() + xlab('Check') + ylab('Original') + ggtitle(paste(study,'TTestP')) + geom_abline()
  print( plot_grid(p1,p2,p3) )
  
  
  
}