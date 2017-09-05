# Script to compare Celeste's data to ours

rm(list = ls())
setwd('/Users/ereznik/pancanmet_analysis/analysis/')

library(ggplot2)

#############################################################################
# Functions
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}
#############################################################################

# Read in Simon data
simon = read.csv('../data/validation/Simon_Metabolon_FC.csv',header = T,row.names = 1)
rownames(simon) = unlist(lapply(rownames(simon),simpleCap))
simon$Color = 'Significant'
simon[which(simon$Qval>0.05),'Color'] = 'NS'
simon$FC = log2(simon$FC)
rownames(simon) = tolower(rownames(simon))

# Read in our data
pancan = read.csv('../results/diffabundance/DifferentialAbundanceSummary.csv',header = TRUE,row.names = 1,stringsAsFactors = FALSE)
pancan = pancan[,c('FC.KIRC','Padj.KIRC')]
pancan$Color = 'Significant'
pancan[which(pancan[,2]>0.05),'Color'] = 'NS'
ixnames = intersect(rownames(pancan),rownames(simon))

fc = data.frame( pancan[ixnames,1], simon[ixnames,1],pancan[ixnames,'Color'],simon[ixnames,'Color'])
colnames(fc) = c('MSKCC','UPenn','MSKColor','UPennColor')
rownames(fc) = ixnames
fc$Significance = 'NS'
fc[which(fc$MSKColor == 'Significant' & fc$UPennColor == 'Significant'),'Significance'] = 'Significant'

fc$MSKSign = NA
fc$PennSign = NA

fc[which(fc$MSKCC > 0 & fc$MSKColor == 'Significant'),'MSKSign'] = 1
fc[which(fc$MSKCC < 0 & fc$MSKColor == 'Significant'),'MSKSign'] = -1

fc[which(fc$UPenn > 0 & fc$UPennColor == 'Significant'),'PennSign'] = 1
fc[which(fc$UPenn < 0 & fc$UPennColor == 'Significant'),'PennSign'] = -1
fisherp = fisher.test(table(fc$MSKSign,fc$PennSign))$p.value
corval = cor.test(fc$MSKCC,fc$UPenn,method = 'spearman')

titlestr = paste0('Kidney, Fisher P-value <1e-16\nSpearman rho ',round(corval$estimate,3),' p-value <1e-16')
ggplot(fc,aes(MSKCC,UPenn,color = Significance)) + geom_point() + theme_classic(20) + 
  xlab('Hakimi et al, 2016') + ylab('Li et al, 2014') + 
  scale_color_manual(values = c('NS' = 'gray','Significant' = 'red')) + 
  ggtitle(titlestr) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  ggsave('../results/validation/ClearCellValidation.pdf',height = 8,width = 8)
print(cor.test(fc$MSKCC,fc$UPenn),method= 'spearman')
