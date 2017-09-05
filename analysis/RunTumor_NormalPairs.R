# Script to analyze co-variation of metabolites in tumor/normal samples of common tissues

rm(list = ls())
setwd('/Users/ereznik/pancanmet_analysis/analysis/')
source('functions/readBigMet.R')
source('plottingconventions.R')
library(ggplot2)
library(ggrepel)
library(cowplot)

# Parameters to set
study2drop = c()
pthresh = 0.1

# Read in data
metdata = readBigMet( study2drop )
met = metdata$met; tissuetype = metdata$tissuetype; studytype = metdata$studytype
uqstudy = unique(studytype)

# Change to matrix
met = as.matrix(met)

# Read in tumor/normal paired data
pairs = read.csv('../data/merged_metabolomics/tumor_normal_pairs.csv',header = TRUE,check.names = FALSE)

# Get roots of names
normalIDs = sapply(1:dim(pairs)[1],function(x){paste(pairs[x,3],pairs[x,1],'Normal',sep=':')})
tumorIDs = sapply(1:dim(pairs)[1],function(x){paste(pairs[x,3],pairs[x,2],'Tumor',sep=':')})

# Get the kind of study
uqstudy = unique(pairs$Study)

# Initialize a vector to store results in
corrR = matrix(NA,dim(met)[1],length(uqstudy)); rownames(corrR) = rownames(met); colnames(corrR) = uqstudy
corrP = corrR
corrPadj = corrR
# For each study
for (study in uqstudy){
  studyix = which(pairs$Study == study)
  tum2use = tumorIDs[studyix]
  norm2use = normalIDs[studyix]
  
  d2use = met[,c(tum2use,norm2use)]
  d2use = d2use[which(complete.cases(d2use)),]
  
  # Calculate correlation
  for (metabolite in rownames(d2use) ){
    tempcor = cor.test(d2use[metabolite,tum2use],d2use[metabolite,norm2use],method = 'spearman')
    corrR[ metabolite,study ] = tempcor$estimate
    corrP[ metabolite, study ] = tempcor$p.value
    # qplot(met[metabolite,tum2use],met[metabolite,norm2use])
  }
  
  # Adjust for hypothesis testing
  corrPadj[,study] = p.adjust(corrP[,study],method = 'BH')
}

# Remove insignificant results
corrR[ which(corrPadj > pthresh, arr.ind = TRUE)] = 0
signR = sign(corrR)

# Analyze results

# Count number of NAs in each row, drop any row with more than 3
numNA = apply(corrR,1,function(x){length(which(is.na(x)))})
keepmets = names(numNA)[ which(numNA <= 3) ]

pairscores = data.frame( matrix(NA,length(keepmets),4) )
rownames(pairscores) = keepmets
colnames(pairscores) = c('NumStudies','NumPositiveCor','NumNegativeCor','NumTotalCor')
for (metabolite in keepmets){
  
  pairscores[metabolite,'NumStudies'] = length(uqstudy) - numNA[metabolite]
  pairscores[metabolite,'NumPositiveCor'] = length(which(corrR[metabolite,]>0))
  pairscores[metabolite,'NumNegativeCor'] = length(which(corrR[metabolite,]<0))
  pairscores[metabolite,'NumTotalCor'] = pairscores[metabolite,'NumPositiveCor'] + pairscores[metabolite,'NumNegativeCor']

}

# Sort results and write to file
pairscores = pairscores[order(pairscores$NumTotalCor,decreasing = TRUE),]
write.csv(pairscores,'../results/pairedsamples/pairscores.csv')
pairscores$Index = dim(pairscores)[1]:1
pairscores$Label = rownames(pairscores)
pairscores[which(pairscores$NumTotalCor <= 4),'Label'] = ''
# Plot the results, annotating the top few metabolites
p1 = ggplot(pairscores,aes(Index,NumTotalCor,label = Label)) + 
  ylab('Number of Significant Normal/Tumor Correlations') + 
  theme_bw(base_size = 12) + geom_point() + geom_text_repel(force = 5) + 
  ggsave('../results/pairedsamples/PairedTumorNormal.pdf',height = 6,width = 6)

# Make a bunch of plots
for (metabolite in rownames(pairscores)){
  # skip the boring ones
  if (pairscores[metabolite,'NumTotalCor'] < 3){next}
  
  # Make a temporary data frame of tumor/normal abundance and study
  tempmet = data.frame( met[metabolite,tumorIDs], met[metabolite,normalIDs] )
  tempmet$Study = sapply(rownames(tempmet),function(x){strsplit(x,'\\:')[[1]][1]})
  tempmet = tempmet[which(complete.cases(tempmet)),]
  colnames(tempmet) = c('Tumor','Normal','Study')
  tempmet$StudyName = names2plot[tempmet$Study]
  tempmet$Tissue = sourcetissue[tempmet$Study]
  
  # Plot
  pmet = ggplot(tempmet,aes(log2(Normal),log2(Tumor),color = Tissue)) + geom_point() + 
    theme_bw(base_size = 16) + facet_wrap(~ StudyName) + 
    ggsave(paste0('../results/pairedsamples/TumorNormalplots/',metabolite,'.pdf'),height =8,width = 8)
  
  if (metabolite == 'alanine'){
    pmet2save = pmet
  }
  
}


