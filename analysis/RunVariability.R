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
  
  tumdisp = apply(tumor,1,mad)
  normdisp = apply(normal,1,mad)
  
  # Also calculate for each metabolite the fraction of imputed samples
  tumimp = apply(tumor,1,function(x){ length(which(abs(x - min(x))<1e-3))/length(x) })
  normalimp = apply(normal,1,function(x){ length(which(abs(x - min(x))<1e-3))/length(x) })
  
  tempdisp = data.frame(tumdisp,normdisp,names(tumdisp),rep(study,length(tumdisp)),tumimp,normalimp)
  colnames(tempdisp) = c('Tumor','Normal','Metabolite','Study','TumorImputed','NormalImputed')
  tempdisp$Ratio = log2(tempdisp$Tumor/tempdisp$Normal)
  
  
  # Statistical test
  statres[study,'GreaterP'] = wilcox.test(tempdisp$Ratio,alternative = 'greater')$p.value
  statres[study,'LessP'] =  wilcox.test(tempdisp$Ratio,alternative = 'less')$p.value
  statres[study,'P'] = wilcox.test(tempdisp$Ratio)$p.value
  
  # Calculate mean on non-infinite cases and non-NA cases
  keepix = which(!is.na(tempdisp$Ratio) & !is.infinite(tempdisp$Ratio))
  statres[study,'MeanMADRatio'] = mean(tempdisp[keepix,'Ratio'],na.rm = TRUE)
  
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

for (dname in list('res','res_lowimp')){
  if (dname == 'res'){
    d = res
    file2save = '../results/Variation/MAD.pdf'
  }
  
  if (dname == 'res_lowimp'){
    d = res_lowimp
    file2save = '../results/Variation/MAD_lowimputation.pdf'
  }
  
  ggplot(d,aes(Ratio)) + facet_grid(Name~.) + 
    theme_bw() + geom_vline(xintercept = 0) + ylab('Number of Metabolites') + 
    xlab('Log2 Ratio of Median Absolute Deviation, Tumor/Normal') +
    geom_histogram(aes(fill = Tissue),alpha = 0.5,binwidth = 0.2) + 
    geom_density(aes(y=0.2*..count..,fill = Tissue),alpha=0.5) +
    facet_grid(Name~.,scales = 'free',drop = FALSE) +
    theme_classic() + geom_vline(xintercept = 0,linetype = "longdash") + 
    theme(legend.position = 'none') +
    theme(strip.background = element_rect(colour="white", fill="white") ) + 
    ggsave(file2save,height = 10,width = 10)

}


# Write out statres to a file
write.csv(statres,'../results/variation/MAD_results.csv')
  