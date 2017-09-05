# Script to analyze covariation of metabolites

rm(list=ls())
setwd('/Users/ereznik/pancanmet_analysis/analysis/')
source('functions/covariation.R')
source('functions/readBigMet.R')
source('plottingconventions.R')
library(ggplot2)

# Parameters to set
#study2drop = c('BLCA')
study2drop = c()
pcrit = NA#0.05/2500
dropNA = FALSE
minsamples = 20

# Read in data
metdata = readBigMet( study2drop )
met = metdata$met; tissuetype = metdata$tissuetype; studytype = metdata$studytype
uqstudy = unique(studytype)

# Read in adjacency of metabolites
metadj = read.csv('../results/keggAdjacencyMatrix/keggMetabMetabAdjacencyMatrix.csv',header = TRUE,row.names = 1,check.names = FALSE)
chebidict = read.csv('../import/tempdir/ChEBI_dictionary.csv',header = TRUE,row.names = 2,stringsAsFactors = FALSE)
metadj_rnames = chebidict[rownames(metadj),'FinalName']

# List of metabolites whose correlations we want to know across all studies
pairs2save = read.csv('../data/settings/Covariation_Metabolites2Save.csv',header = TRUE,row.names = 1)

# Don't use duplicated names
uq_rnames = which(!duplicated(metadj_rnames))
metadj2 = metadj[uq_rnames,uq_rnames]

# Set the diagonal to zero
diag(metadj2) = 0
rownames(metadj2) = chebidict[rownames(metadj2),'FinalName']
colnames(metadj2) = chebidict[colnames(metadj2),'FinalName']

# Find adjacent metabolites
adjpairs = which(metadj2 == 1,arr.ind = TRUE)
metpairstrings = paste( rownames(metadj2)[adjpairs[,1]], colnames(metadj2)[adjpairs[,2]],sep = ':' )
  
# Decide whether we will only use metabolites for which we have all measurements or not
if (dropNA){
  met = met[which(complete.cases(met)),]
}

# Initialize a vector to store results in
corrResults = matrix(NA,1,1)

# For each study, calculate the covariation of all metabolite pairs
res = data.frame()
saver = data.frame()
savep = data.frame()
for (study in uqstudy){
  print(study)
  for (tissue in c('Normal','Tumor')){
   
    studyix = which(studytype == study & tissuetype == tissue)
    if (length(studyix) < minsamples){
      next
    }
    savename = paste(study,':',tissue)
    returnlist = covariation( met[,studyix], pcrit )
    r = returnlist$r; p = returnlist$p; rlist = returnlist$rlist
    
    if (is.na( corrResults[1,1] )){
      corrResults = data.frame( rlist )
      colnames(corrResults) = savename
    }else{
      cnames = colnames(corrResults)
      corrResults = cbind( corrResults, rlist[rownames(corrResults)] )
      colnames(corrResults) = c(cnames,savename)
    }
    
    # Save the correlations and p values we are interested in
    saver[rownames(pairs2save),paste(study,tissue,sep = ':')] = apply(pairs2save,1,function(x){r[x[1],x[2]]})
    savep[rownames(pairs2save),paste(study,tissue,sep = ':')] = apply(pairs2save,1,function(x){p[x[1],x[2]]})
    
    # Drop any NA metabolites
    diag(r) = NA
    namet = which( apply(r,1,function(x){all(is.na(x))}) )
    rdrop = r[-namet,-namet]
    
    # Now intersect with the metabolites for which we have calculated adjacency
    ixmets = intersect(rownames(rdrop),rownames(metadj2))
    rdrop_adj = rdrop[ixmets,ixmets]
    metadj3 = metadj2[ixmets,ixmets]
    
    # Now, separate correlations based on whether they are adjacent or not
    corr_adj = rdrop_adj[which(metadj3 == 1,arr.ind = TRUE)]
    corr_notadj = rdrop_adj[which(metadj3 == 0,arr.ind = TRUE)]
    
    # Make a dataframe to plot
    pdata1 = data.frame(corr_adj)
    colnames(pdata1) = c('Correlation')
    pdata1$Adjacency = 'Adjacent'
    
    pdata2 = data.frame(corr_notadj)
    colnames(pdata2) = c('Correlation')
    pdata2$Adjacency = 'NotAdjacent'
    
    pdata = rbind(pdata1,pdata2)
    print( ggplot(pdata,aes(Correlation,fill = Adjacency)) + geom_density(alpha = 0.5) + 
             ggtitle(paste(names2plot[study],tissue)) + theme_bw(base_size = 14) + 
             xlab('Spearman Correlation Coefficient') + ylab('Density') + 
             xlim(-1,1) + 
             ggsave(paste0('../results/covariation/',study,'_',tissue,'.pdf'),height = 6,width = 6))
    
    res[paste(study,tissue),'Mean_NonAdjacent'] = mean(pdata2$Correlation,na.rm = TRUE)
    res[paste(study,tissue),'Mean_Adjacent'] = mean(pdata1$Correlation,na.rm = TRUE)
    res[paste(study,tissue),'Pvalue'] = wilcox.test(pdata1$Correlation,pdata2$Correlation)$p.value
      
  }
  
}

write.csv(res,'../results/covariation/CCM_Covariation_Results.csv')

# Remove insignificant pairwise correlations
saver[which(savep > pthresh,arr.ind = TRUE)] = 0
saver$Metabolite1 = pairs2save$Metabolite1
saver$Metabolite2 = pairs2save$Metabolite2
