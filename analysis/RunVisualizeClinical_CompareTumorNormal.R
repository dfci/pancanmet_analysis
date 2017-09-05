# Compare clinical results from tumor and normal samples of same study

rm(list=ls())
setwd('/Users/ereznik/pancanmet_analysis/analysis/')
source('functions/covariation.R')
source('functions/readBigMet.R')
source('plottingconventions.R')
library(ggplot2)
library(XLConnect)
library(reshape)

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}


# Make a study keep, note that each clinical type has only one replicate cancer type each
stagekey = c('KIRC','BRCA')
gradekey = c('KIRC','BRCA','PRAD','PRADLODA')

# Read in the clinical data
wb = loadWorkbook('../results/clinical/updated results stage and grade Jan 2017 max1.xlsx')

stage = readWorksheet(wb,header = TRUE,rownames = 1,sheet = 'stage',startRow = 2)
stagenames = rownames(stage)
stage = sapply(stage[,1:dim(stage)[2]], function(x){as.numeric(as.character(x))} )
rownames(stage) = stagenames

grade = readWorksheet(wb,header = TRUE,rownames = 1,sheet = 'Grade',startRow = 2)
gradenames = rownames(grade)
grade = sapply(grade[,1:dim(grade)[2]], function(x){as.numeric(as.character(x))} )
rownames(grade) = gradenames

res = data.frame()

# for the stage and grade files
for (dtype in c('stage','grade')){
  keyfile = get(paste0(dtype,'key'))
  d = get(dtype)
 
  for (study in keyfile){
    
    # Compare the relevant tumor/normal columns
    tnamepos = paste(study,'Tumor','less',sep = '.')
    nnamepos = paste(study,'Normal','less',sep = '.')
    
    tnameneg = paste(study,'Tumor','greater',sep = '.')
    nnameneg = paste(study,'Normal','greater',sep = '.')
    
    # Find significant metabolites
    posmet = which(d[,tnamepos] < pthresh & d[,nnamepos] < pthresh & (d[,'all.combined.adjusted'] < pthresh | d[,'tumors.combined.adjusted'] < pthresh))
    negmet = which(d[,tnameneg] < pthresh & d[,nnameneg] < pthresh & (d[,'all.combined.adjusted'] < pthresh | d[,'tumors.combined.adjusted'] < pthresh))
    
    posmet = rownames(d)[posmet]
    negmet = rownames(d)[negmet]
    
    # Print
    print(paste(dtype,study))
    print(rownames(d)[posmet])
    print(rownames(d)[negmet])
    
    # Save
    ii1 = paste0(study,dtype,'Positive',sep =':')
    res[ii1,'Study'] = names2plot[study]
    res[ii1,'ClinicalVariable'] = simpleCap(dtype)
    res[ii1,'Type'] = 'Positive Correlation'
    res[ii1,'Number'] = length(posmet)
    res[ii1,'Details'] = paste(posmet,collapse = ';')
    
    ii2 = paste0(study,dtype,'Negative',sep =':')
    res[ii2,'Study'] = names2plot[study]
    res[ii2,'ClinicalVariable'] = simpleCap(dtype)
    res[ii2,'Type'] = 'Negative Correlation'
    res[ii2,'Number'] = length(negmet)
    res[ii2,'Details'] = paste(negmet,collapse = ';')
    
    
  } 
  
}

# Make a label
res$Label = paste0('n=',res$Number)

ggplot(res,aes(Study,Number,fill = Type,label = Label)) + geom_bar(stat = 'identity',position = 'dodge') +
  facet_wrap(~ClinicalVariable,ncol = 1) + theme_bw(base_size = 12) + theme(legend.position  = 'bottom') + 
  xlab(NULL) + ylab('# of Metabolites With Common \nAssociation in Tumor/Normal Samples') + 
  geom_text(position = position_dodge(width = 1)) + 
  ggsave('../results/clinical/TumorNormalSameDirection.pdf',height = 6,width = 6)
