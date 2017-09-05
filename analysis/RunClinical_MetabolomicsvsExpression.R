# Script to compare association with stage using metabolomics and expression data

rm(list=ls())
setwd('/Users/ereznik/pancanmet_analysis/analysis/')
source('functions/covariation.R')
source('functions/readBigMet.R')
source('plottingconventions.R')
library(ggplot2)
library(XLConnect)

# Set p value threshold
pthresh = 0.05

# Read in adjacency of metabolites
metadj = read.csv('../results/keggAdjacencyMatrix/keggMetabGeneAdjacencyMatrix.csv',header = TRUE,check.names = FALSE)

# Remove duplicates
duprnames = which(duplicated(metadj[,1]))
if (length(duprnames) > 0){
  metadj = metadj[-duprnames,]
}
rownames(metadj) = metadj[,1]
metadj = metadj[,-1]

# Make a study keep
studykey = c('KIRC' = 'kirc','BRCA' = 'breast','BRCATang' = 'breast')
#studykey = c('KIRC' = 'kirc')

# Read in the clinical data
wb = loadWorkbook('../results/clinical/updated results stage and grade Jan 2017 max1.xlsx')
stage = readWorksheet(wb,header = TRUE,rownames = 1,sheet = 'Stage',startRow = 2)
stagenames = rownames(stage)
stage = sapply(stage[,1:dim(stage)[2]], function(x){as.numeric(as.character(x))} )
rownames(stage) = stagenames

ixmets = intersect(stagenames,rownames(metadj))

for (study in names(studykey)){
  
  print(study)
  res = data.frame()
  resmet = data.frame()
  
  # Read in results
  gres = read.csv(paste0('../results/clinical/stage-limma-diff-analysis/',studykey[study],'-limma-RNAseq-output.txt'),
                  row.names = 1,stringsAsFactors = FALSE,header = TRUE,sep = '\t')
  
  if (study == 'KIRC'){
    # Check against our prior work
    oldfc = read.csv('/Users/ereznik/Documents/ccrc/finaldata/Aug172015/metabologram_stage/CCRCC_gene.csv',header = TRUE,row.names =1 )
    h2e = read.csv('/Users/ereznik/Documents/useful/hugo2entrez.txt',header = TRUE,sep = '\t')
    h2e = h2e[-which(duplicated(h2e$Entrez.Gene.ID)),]
    h2e = h2e[-which(is.na(h2e$Entrez.Gene.ID)),]
    rownames(h2e) = h2e$Entrez.Gene.ID
    oldfc$Hugo = h2e[rownames(oldfc),'Approved.Symbol']
    oldfc$Entrez = rownames(oldfc)
    oldfc = oldfc[-which(duplicated(oldfc$Hugo) | is.na(oldfc$Hugo)),]
    rownames(oldfc) = oldfc$Hugo
    ixgfc = intersect(rownames(oldfc),rownames(gres))
    qplot(oldfc[ixgfc,1],gres[ixgfc,'logFC'])
  }
  
  for (met in ixmets){
    
    # Find the p value and sign
    pdown = stage[met,paste0(study,'.Tumor.greater')]
    pup = stage[met,paste0(study,'.Tumor.less')]
    
    if (is.na(pup) & is.na(pdown)){
      next
    }
    
    if (pup <= pdown){
      metsign = 'Greater in High Stage'
    }else if (pdown < pup){
      metsign = 'Lower in High Stage'
    }else{
      stop('Error in getting p values.')
    }
    
    minmetp = min(c(pup,pdown),na.rm = TRUE)
    if (minmetp < pthresh){
      metsig = TRUE
    }else{
      metsig = FALSE
    }
      
    # Find the adjacent genes
    ixg = colnames(metadj)[which(metadj[met,] == 1)]
    ixg = intersect(ixg,rownames(gres))
    
    for (g in ixg){
      
      resix = paste(met,g,sep = ':')
      res[resix,'Metabolite'] = met
      res[resix,'Gene'] = g
      res[resix,'Study'] = study
      res[resix,'Metabolite.Sig'] = metsig
      res[resix,'Metabolite.P'] = minmetp
      res[resix,'Metabolite.Sign'] = metsign
      res[resix,'Gene.FC'] = gres[g,'logFC']
      res[resix,'Gene.P.Adj'] = gres[g,'adj.P.Val']
      
      if (res[resix,'Gene.P.Adj'] < pthresh){
        res[resix,'Gene.Sig'] = TRUE
      }else{
        res[resix,'Gene.Sig'] = FALSE
      }
      
    }
    
    # Also calculate the "average" score around each metabolite
    resmetix = paste(met,study = ':')
    resmet[resmetix,'Metabolite'] = met
    resmet[resmetix,'Study'] = study
    resmet[resmetix,'Metabolite.Sig'] = metsig
    resmet[resmetix,'Metabolite.P'] = minmetp
    resmet[resmetix,'Metabolite.Sign'] = metsign
    
    gmetix = which(res$Metabolite == met)
    gsig = which(res[gmetix,'Gene.Sig'] == TRUE)
    gnsig = which(res[gmetix,'Gene.Sig'] == FALSE)
    resmet[resmetix,'GeneScore'] = length(gsig)/(length(gnsig) + length(gsig))

  }
  
  # Make sure we don't count genes which are connected to both significant and insignificant metabolites
  sigmets = unique( res[which(res$Metabolite.P <= pthresh),'Metabolite'] )
  notsigmets = unique( res[which(res$Metabolite.P > pthresh),'Metabolite'] )
  sigg_connected = unique(res[which(res$Metabolite %in% sigmets),'Gene'])
  notsigg_connected = unique(res[which(res$Metabolite %in% notsigmets),'Gene'])
  rmgenes = intersect(sigg_connected,notsigg_connected)
  
  res2 = res[-which(res$Gene %in% rmgenes),]
  
  # Analyze the results (choose res or res2)
  sigtable = table(res$Metabolite.Sig,res$Gene.Sig)
  print(sigtable)
  print(chisq.test(sigtable))
  
}
