# Compare clinical results across tumor types with replicate studies

rm(list=ls())
setwd('/Users/ereznik/pancanmet_analysis/analysis/')
source('functions/covariation.R')
source('functions/readBigMet.R')
source('plottingconventions.R')
library(ggplot2)
library(XLConnect)
library(reshape)

# Set p value threshold
pthresh = 0.05

# Make a function to annotate significantly changed metabolites
sigchange = function(x,group,string1,string2){
  # x is data frame
  # string1 is string to print if first column is more significant
  # string2 is string to print if second column is more significant
  for (rname in rownames(x)){
    minp = min(x[rname,c(1,2)],na.rm = TRUE)
    x[rname,'p'] = minp
    if (minp > pthresh){
      x[rname,group] = 'NotSignificant'
    }else if(x[rname,1] == minp ){
      x[rname,group] = string1
    }else{
      x[rname,group] = string2
    }
  }
  return(x)
}

# Make a study keep, note that each clinical type has only one replicate cancer type each
stagekey = list('Breast' = c('BRCA','BRCATang'))
gradekey = list('Prostate' = c('PRAD','PRADLODA'))

# Read in the clinical data
wb = loadWorkbook('../results/clinical/updated results stage and grade Jan 2017 max1.xlsx')

stage = readWorksheet(wb,header = TRUE,rownames = 1,sheet = 'Stage',startRow = 2)
stagenames = rownames(stage)
stage = sapply(stage[,1:dim(stage)[2]], function(x){as.numeric(as.character(x))} )
rownames(stage) = stagenames

grade = readWorksheet(wb,header = TRUE,rownames = 1,sheet = 'Grade',startRow = 2)
gradenames = rownames(grade)
grade = sapply(grade[,1:dim(grade)[2]], function(x){as.numeric(as.character(x))} )
rownames(grade) = gradenames

# for the stage and grade files
for (dtype in c('stage','grade')){
  keyfile = get(paste0(dtype,'key'))
  d = get(dtype)
  
  # For each data type with replicate data, 
  # identify metabolites with significant associations
  t1 = keyfile[[1]][1]
  t2 = keyfile[[1]][2]
  dtumor1 = d[,c(paste0(t1,'.Tumor.greater'),paste0(t1,'.Tumor.less'))]
  dtumor2 = d[,c(paste0(t2,'.Tumor.greater'),paste0(t2,'.Tumor.less'))]
  tumormets = intersect(rownames(dtumor1)[which(complete.cases(dtumor1))],rownames(dtumor2)[which(complete.cases(dtumor2))])
  
  dtumor1 = data.frame( dtumor1[tumormets,] )
  dtumor2 = data.frame( dtumor2[tumormets,] )
  
  dtumor1 = sigchange(dtumor1,t1,paste0('Negatively correlated to',dtype),paste0('Positively correlated to ',dtype))
  dtumor2 = sigchange(dtumor2,t2,paste0('Negatively correlated to',dtype),paste0('Positively correlated to ',dtype))
  
  # Tabulate
  restumor = data.frame(dtumor1[,t1],dtumor2[,t2],stringsAsFactors = FALSE)
  colnames(restumor) = c(t1,t2)
  restumor[,t1] = paste(t1,restumor[,t1],sep = '\n')
  restumor[,t2] = paste(t2,restumor[,t2],sep = '\n')
  rownames(restumor) = rownames(dtumor1)
  
  if (dtype == 'grade'){
    dnormal1 = d[,c(paste0(t1,'.Normal.greater'),paste0(t1,'.Normal.less'))]
    dnormal2 = d[,c(paste0(t2,'.Normal.greater'),paste0(t2,'.Normal.less'))]
    normalmets = intersect(rownames(dnormal1)[which(complete.cases(dnormal1))],rownames(dnormal2)[which(complete.cases(dnormal2))])
   
    dnormal1 = data.frame( dnormal1[normalmets,] )
    dnormal2 = data.frame( dnormal2[normalmets,] )
    
    dnormal1 = sigchange(dnormal1,t1,paste0('Negatively correlated to ',dtype),paste0('Positively correlated to ',dtype))
    dnormal2 = sigchange(dnormal2,t2,paste0('Negatively correlated to ',dtype),paste0('Positively correlated to ',dtype))
   
    resnormal = data.frame(dnormal1[,t1],dnormal2[,t2])
    colnames(resnormal) = c(t1,t2)
  }
  
  
  # Identify the significantly changed metabolites
  tumortable = as.data.frame.matrix( table(restumor[,t1],restumor[,t2]) )
  tumortable$Name = rownames(tumortable)
  tumortable2 = melt(tumortable,id.vars = 'Name')
  ggplot(tumortable2,aes(Name,variable,fill = value,label = value)) + geom_tile() + geom_text() + 
    scale_fill_gradient(low = 'white',high = 'lightblue') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    ggsave(paste0('../results/clinical/SameTissue/',t1,'_',t2,'.pdf'),height = 8,width = 8)
  print(tumortable)
  print(chisq.test(table(restumor[,t1],restumor[,t2])))
  
  if(dtype == 'grade'){
    print(table(resnormal[,t1],resnormal[,t2]))
  }
  
}
