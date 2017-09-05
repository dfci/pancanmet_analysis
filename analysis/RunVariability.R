# Analyze variability in paired samples of tumor/normal

# Issues: 1. Do we need to worry about imputed data biasing MAD comparison? Doesn't seem like it from testing.

rm(list=ls())
setwd('/Users/ereznik/pancanmet_analysis/analysis/')
library(ggplot2)
library(reshape2)
source('plottingconventions.R')
source('functions/readBigMet.R')

# Read in the paired metabolomics data
metdata = readBigMet(c())
met = metdata[[1]]
tissuetype = metdata[[2]]
tumnormtype = metdata[[3]]

uqtis = unique(tissuetype)

# Read in the pairs
pairs = read.csv('../data/merged_metabolomics/tumor_normal_pairs.csv',check.names = FALSE,header = TRUE,stringsAsFactors = FALSE)

# Compare the results of MAD analysis to differential abundance analysis
da = read.csv('../results/diffabundance/DifferentialAbundanceSummary.csv',
              header = TRUE,row.names = 1,check.names = FALSE,stringsAsFactors = FALSE)

firstflag = TRUE
statres = data.frame()

for (study in uqtis){
  
  # Check if this study appears in the paired data
  if (!(study %in% pairs$Study)){
    next
  }
  
  # Make the paired data
  ix = pairs[pairs$Study == study,]
  tumnames = paste(ix$Study,ix$`Tumor Sample`,'Tumor',sep= ':')
  normnames = paste(ix$Study,ix$`Normal Sample`,'Normal',sep =':')
  
  # Check that all the names are in the met data
  if (length(setdiff(c(tumnames,normnames),colnames(met)))!=0){
    stop('Error in naming pairs.')
  }
  
  studydata = met[,c(tumnames,normnames)]
  studydata = studydata[complete.cases(studydata),]
  tumor = studydata[,tumnames]
  normal = studydata[,normnames]
  
  # Normalize tumor and normal data by the within-cohort median
  tumor_norm = t( apply(tumor,1,function(x){x/median(x)}) )
  normal_norm = t( apply(normal,1,function(x){x/median(x)}) )
  
  # Take the log of this
  tumor_norm = log2(tumor_norm)
  normal_norm = log2(normal_norm)
  
  tumdisp = apply(tumor_norm,1,mad) #only this calculation needs to be run with tumor_norm and normal_norm
  normdisp = apply(normal_norm,1,mad)
  
  # Also calculate for each metabolite the fraction of imputed samples
  tumimp = apply(tumor,1,function(x){ length(which(abs(x - min(x))<1e-3))/length(x) })
  normalimp = apply(normal,1,function(x){ length(which(abs(x - min(x))<1e-3))/length(x) })
  
  tempdisp = data.frame(tumdisp,normdisp,names(tumdisp),rep(study,length(tumdisp)),tumimp,normalimp)
  colnames(tempdisp) = c('Tumor','Normal','Metabolite','Study','TumorImputed','NormalImputed')
  tempdisp$Ratio = log2(tempdisp$Tumor/tempdisp$Normal)
  
  # Add in differential abundance data
  tempdisp[,'FC'] = da[rownames(tempdisp),paste('FC',study,sep = ' ')]
  tempdisp[,'Padj'] = da[rownames(tempdisp),paste('Padj',study,sep = ' ')]
  
  # Statistical test
  statres[study,'GreaterP'] = wilcox.test(tempdisp$Ratio,alternative = 'greater')$p.value
  statres[study,'LessP'] =  wilcox.test(tempdisp$Ratio,alternative = 'less')$p.value
  statres[study,'P'] = wilcox.test(tempdisp$Ratio)$p.value
  
  # Calculate mean on non-infinite cases and non-NA cases
  keepix = which(!is.na(tempdisp$Ratio) & !is.infinite(tempdisp$Ratio))
  statres[study,'MeanMADRatio'] = mean(tempdisp[keepix,'Ratio'],na.rm = TRUE)
  
  # Modify rownames of tempdisp
  rownames(tempdisp) = paste(study,rownames(tempdisp),sep = ':')
  if (firstflag){
    res = tempdisp
    firstflag = FALSE
  }else{
    res = rbind(res,tempdisp)
  }
  
}

res$Name = names2plot[as.character(res$Study)]
res$Tissue = sourcetissue[as.character(res$Study)]

# Repeat plot using only metabolite which are not heavily imputed
res_lowimp = res[which(res$TumorImputed < 0.2 & res$NormalImputed < 0.2),]

# # Cast the data
# rescast = dcast(res_lowimp,Metabolite ~ Study, value.var = 'Ratio')
# rownames(rescast) = rescast$Metabolite
# rescast = rescast[,-which(colnames(rescast) == 'Metabolite')]
# 
# # Count the number of NAs and the average of the ratio for each metabolite
# rescast$notNA = NA
# rescast$AverageMAD = NA
# for (metabolite in rownames(rescast)){
#   naix = which(!is.na(rescast[metabolite,]))
#   rescast[metabolite,'notNA'] = length(naix)
#   if (length(naix) == 0){
#     next
#   }
#   rescast[metabolite,'AverageMAD'] = mean( as.matrix(rescast[metabolite,naix]), na.rm = TRUE )
#   
# }
#
# # Analyze rescast
# rescast_trim = rescast[which(rescast$notNA > 5),]
# rescast_trim = rescast_trim[order(rescast_trim$AverageMAD),]

# For each study, calculate the correlation between MAD and FC
tempcor = list()
for (study in unique(res$Study)){
  tempres = res[which(res$Study == study),]
  tempcor[[study]] = cor.test(tempres$Ratio,tempres$FC,method = 'spearman')
}

for (dname in list('res','res_lowimp')){
  if (dname == 'res'){
    d = res
    file2save = '../results/variation/MAD.pdf'
    file2save2 = '../results/variation/MAD_vs_DA.pdf'
  }
  
  if (dname == 'res_lowimp'){
    d = res_lowimp
    file2save = '../results/variation/MAD_lowimputation.pdf'
    file2save2 = '../results/variation/MAD_vs_DA_lowimputation.pdf'
  }
  
  ggplot(d,aes(Ratio)) +
    theme_bw(base_size = 20) + geom_vline(xintercept = 0) + ylab('Number of Metabolites') + 
    xlab('Log2 Ratio of Median Absolute Deviation, Tumor/Normal') +
    geom_histogram(aes(fill = Tissue),alpha = 0.5,binwidth = 0.2) + 
    geom_density(aes(y=0.2*..count..,fill = Tissue),alpha=0.5) +
    facet_grid(Name~.,scales = 'free',drop = FALSE) +
    theme_classic() + geom_vline(xintercept = 0,linetype = "longdash") + 
    theme(legend.position = 'none') + xlim(-4,4) + 
    theme(strip.background = element_rect(colour="white", fill="white") ) + 
    scale_fill_manual(values = sourcetissue_color) + 
    ggsave(file2save,height = 10,width = 6)
    
  
  ggplot(d,aes(FC,Ratio)) + geom_point() + facet_wrap(~Name) + theme_bw() +
    xlab('Log2 Ratio, Tumor:Normal') + ylab('Normalized MAD Ratio') +
    ggsave(file2save2,height = 10,width = 10) 

}

# Write out statres to a file
write.csv(statres,'../results/variation/MAD_results.csv')

  