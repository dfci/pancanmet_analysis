# Correlate tumor purity with metabolite levels

rm(list = ls())
setwd('/Users/ereznik/pancanmet_analysis/analysis/')
library(ggplot2)
library(XLConnect)
source('plottingconventions.R')
source('functions/readBigMet.R')

# Read in metabolites file for BRCA Tang
d = readBigMet(c())
tangnames = which(d$studytype == 'BRCATang')
met = d$met
met = met[,tangnames]
met = met[which(complete.cases(met)),]
metcnames = sapply(colnames(met),function(x){strsplit(x,'\\:')[[1]][2]})
colnames(met) = metcnames
met = as.matrix(met)

# Read in estimate scores
estimate = readWorksheetFromFile('../data/ESTIMATE_TCGA.xlsx',sheet = 'RNASeqV2',header = TRUE,startRow = 2)
estimate$Sample = sapply(estimate$ID,function(x){y = strsplit(x,'\\.')[[1]];paste(y[2],y[3],sep = '-')})
rownames(estimate) = estimate$Sample

# Intersect samples
ixsamps = intersect(rownames(estimate),colnames(met))

# Correlate
res = data.frame()
for (metabolite in rownames(met)){
  temps = cor.test(met[metabolite,ixsamps],estimate[ixsamps,'Stromal.score'],method = 'spearman')
  tempi = cor.test(met[metabolite,ixsamps],estimate[ixsamps,'Immune.score'],method = 'spearman')
  
  res[metabolite,'StromalR'] = temps$estimate
  res[metabolite,'StromalP'] = temps$p.value
  
  res[metabolite,'ImmuneR'] = tempi$estimate
  res[metabolite,'ImmuneP'] = tempi$p.value
}

res$StromalPadj = p.adjust(res$StromalP,method = 'BH')
res$ImmunePadj = p.adjust(res$ImmuneP,method = 'BH')

res$StromalSig = 'NS'
res$ImmuneSig = 'NS'

res[which(res$StromalP < 0.01),'Significance'] = 'Stroma-Correlated'
res[which(res$ImmuneP < 0.01),'Significance'] = 'Immune-Correlated'
res[which(res$ImmuneP < 0.01 & res$StromalP < 0.01),'Significance'] = 'Stroma and Immune-Correlated'
sigmets = which(res$Significance != 'NS')
res$Label = ''
res[sigmets,'Label'] = rownames(res)[sigmets]

ggplot(res,aes(StromalR,ImmuneR,color = Significance,label = Label)) + 
  geom_point() + theme_bw(base_size = 16) + geom_text_repel(show.legend = FALSE) + 
  theme(legend.position = 'bottom',legend.direction = "vertical",nco) +
  ylab('Correlation, Immune Score') + xlab('Correlation, Stromal Score') + 
  ggsave('../results/validation/TumorPurity.pdf',height = 8,width = 8)
