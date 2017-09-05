# A script to calculate pathway enrichment scores for meta-analysis of metabolomics data

rm(list = ls())
setwd('/Users/ereznik/Documents/pancanmet/analysis/pancan/')

library(ggplot2)

thresh = 5e-2
datapath = '/Users/ereznik/Documents/pancanmet/data/pathview/'

################
# Read the table of KEGG pathway annotations
keggdata = read.csv('../../data/pathway_v2/kegg_details_munged.csv',header = F,row.names = 1,colClasses='character')

studies = c('BLCA','BRCA','BRCATang','COAD','KICH','KIRC','OV','PAAD','PAADHussein1','PAADHUSSEIN2','PRAD','PRADLODA','STAD')
studies = 'KIRC'
################

# Make a data array to store results
enrich = matrix(0,dim(keggdata)[1],length(studies))
rownames(enrich) = rownames(keggdata)
colnames(enrich) = studies

for (s in studies){
  
  # Read in the metabolomics data
  m = read.csv(paste(datapath,s, '_met.csv',sep=''), header = F)
  rownames(m) = m[,1]
  
  # For each kegg pathway, make a mapping for the metabolites
  keggmap = matrix(0,dim(keggdata)[1],dim(m)[1])
  rownames(keggmap) = rownames(keggdata)
  colnames(keggmap) = m[,1]
  
  for (path in rownames(keggdata)){
    metsinpath = unlist( strsplit(keggdata[path,2],'\\|') )
    intm =  intersect(metsinpath,m[,1])
    if (length(intm)==0){next}
    keggmap[path,intm] = 1
  }
  
  # Now, for each pathway, test for enrichment
  for (path in rownames(keggdata)){
    inidx = which(keggmap[path,] == 1)
    innames = rownames(m)[inidx]
    inpath_sigchange = length( which(m[innames,2] != 0) )
    inpath_nochange = length( which(m[innames,2] == 0) )
    
    outidx = which(keggmap[path,] == 0)
    outnames = rownames(m)[outidx]
    outpath_sigchange = length( which(m[outnames,2] != 0) )
    outpath_nochange = length( which(m[outnames,2] == 0) )
    
    
    testmat = matrix(c(inpath_sigchange, inpath_nochange, 
                       outpath_sigchange, outpath_nochange),
                     nrow = 2,
                     dimnames = list(Diff = c('SigChange','NoChange'),
                                     Path = c('InPath','NotinPath')) )
    enrich[path,s] = fisher.test( testmat )$p.value
  }
}
