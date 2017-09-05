# A script to analyze the output from Irina's analysis and make some plots

rm(list = ls())
setwd('/Users/ereznik/Documents/pancanmet/analysis/pancan/')
library(ggplot2)
library(data.table)

# Set p-value threshold
pthresh = 0.05
minstudies = 2

dtype = 'Stage' # stage or grade, first letter capitalized

# Clear files from the directory we are printing to
printdir = paste('/Users/ereznik/Documents/pancanmet/results/clinassoc/boxplots/',dtype,'/',sep='')
do.call(file.remove,list(list.files(printdir,full.names=TRUE)))

# Read in the metabolomics data
temp = fread('../../data/merged_metabolomics/alldata.csv',sep = ',',header = F,skip = 3)
alldata = data.frame( temp, row.names = 1 )
title.line = strsplit( readLines('../../data/merged_metabolomics/alldata.csv', n=1), ',' ) [[1]] # Note that first element is empty ""
colnames(alldata) = title.line[2:length(title.line)]

# Read in the clinical data itself
clin = read.csv('/Users/ereznik/Documents/pancanmet/data/merged_metabolomics/clinfeatures.csv',header = T)
for (i in 1:dim(clin)[1]){
  rownames(clin)[i] = paste(clin[i,2],clin[i,1],clin[i,3],sep= ':')
}

# Turn the stage and grade variables into character vectors, otherwise this will mess with later analysis
clin$Grade = as.character(clin$Grade)
clin$Stage = as.character(clin$Stage)

# Read in the relevant clinical correlations results
d = read.csv(paste('/Users/ereznik/Documents/pancanmet/results/clinassoc/',dtype,'_Aug152015.csv',sep=''),row.names = 1)

# # Make a QQ plot
# qq = data.frame( matrix(NA, dim(d)[1],2) )
# colnames(qq) = c('Observed','Expected')
# observed_order = order(d$adjusted.combined.p.TUMOR.ONLY)
# observed = d[observed_order,'adjusted.combined.p.TUMOR.ONLY']
# rownames(qq) = rownames(d)[observed_order]
# qq['Observed'] = -(log10(observed))
# 
# expected = c(1:length(observed)) 
# qq[,'Expected'] = -(log10(expected / (length(expected)+1)))
# 
# ggplot(qq,aes(Expected,Observed)) + geom_point() + geom_abline() + theme_bw()


# Make a results data frame
res = data.frame( matrix(NA, 0,3) )
#rownames(res) = rownames(d)
colnames(res) = c('Adj.P','NumStudies','StudyNames')

# Find the columns of d associated with tumors
tumcols = grep('Tumor',colnames(d)) # Careful here that we don't include the p-value columns!

# Find the metabolites which meet the p-value threshold, and have at least minstudies with significant associations
for (met in rownames(d)){
  if (is.na(d[met,'adjusted.combined.p.TUMOR.ONLY'])){next}
  if (d[met,'adjusted.combined.p.TUMOR.ONLY'] > pthresh){next}
  
  res[met,'Adj.P'] = d[met,'adjusted.combined.p.TUMOR.ONLY']
  
  # Count the number of studies using tumor only data that have a significant p-value
  sigstudies = tumcols[which(d[met,tumcols]<pthresh)] # This is a weird line of code
  sigstudynames = colnames(d)[sigstudies]
  shortnames = unlist(lapply(sigstudynames,function(x){y = strsplit(x,'\\.');y[[1]][1]}))
  res[met,'NumStudies'] = length(sigstudies)
  res[met,'StudyNames'] = paste( shortnames, collapse = '|' )
  
  # If there are at least three studies with a correlation, then plot the result
  if (length(sigstudies)>2){
   
    # Get the data
    plotdata = data.frame( matrix(NA,0,3) )
    colnames(plotdata) = c('Metabolite','Study','Clinical')
  
    for (s in shortnames){
      # Find the tumor samples for which we have clinical data
      ix = which(clin$Study == s & clin$Type == 'Tumor')
      ixnames = rownames(clin)[ix]
      
      # Do an intersection because we might have removed some samples (e.g. ccpap)
      ixnames_met = ixnames[which(ixnames%in%colnames(alldata))]
      
      # Store the clinical data
      plotdata[ixnames_met,'Study'] = s
      plotdata[ixnames_met,'Clinical'] = clin[ixnames_met,dtype]
      
      # Get the metabolite data associated with these samples
      plotdata[ixnames_met,'Metabolite'] = t( alldata[met,ixnames_met] )
      
    }
    
    # Make the clinical data a factor and hope that it's in correct order
    plotdata$Clinical = factor( plotdata$Clinical )
    
    # Remove NAs
    plotdata = plotdata[which(complete.cases(plotdata)),]
    if (dim(plotdata)[1] == 0){next}
    
    # Make a plot
    print( ggplot(plotdata,aes(Clinical,Metabolite)) + geom_violin() + geom_point() + 
             facet_wrap(~Study,scales = 'free') + 
             theme_bw() + xlab(dtype) + ylab(met) + 
             ggsave(paste(printdir,met,'_',dtype,'_pthresh',pthresh,'.pdf',sep=''),height = 4,width = 10) )
  }
  
}

# Also make a summary barchart
restable = data.frame( table(res$NumStudies))
colnames(restable) = c('NumStudies','Count')

# Remove row corresponding to 0 numstudies
restable = restable[-which(restable$NumStudies == 0),]
ggplot(restable,aes(NumStudies,Count)) + geom_bar(stat='identity') + 
 xlab('Number of Studies with Significant Association') + ylim(0,80) + theme_classic() + 
  ylab('Number of Metabolites') + ggtitle(paste('Tumor',dtype)) + theme_set(theme_bw(base_size = 20)) + 
  ggsave(paste('/Users/ereznik/Documents/pancanmet/results/clinassoc/boxplots/',dtype,'_summary.pdf',sep=''),height = 10,width = 10)
