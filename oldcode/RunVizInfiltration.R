# Script to visualize metabolomics correlations with RNA-Seq estimate scores

rm(list = ls())
setwd('/Users/ereznik/pancanmet_analysis/analysis/')
library(ggplot2)
library(ggrepel)

cvals = c('Stromal:NS\nImmune:NS' = 'gray','Stromal:NS\nImmune:Significant' = 'blue',
          'Stromal:Significant\nImmune:NS' = 'red')

r = read.csv('../results/immuneInfiltrationAllColsCors.csv',header = TRUE,row.names = 1)
p = read.csv('../results/immuneInfiltrationAllColsPVals.csv',header = TRUE,row.names = 1)

# Get data to plot
pdata = r[,1:2]
pdata$StromalP = p$StromalScore
pdata$ImmuneP = p$ImmuneScore

pdata$StromalColor = 'Significant'
pdata$ImmuneColor = 'Significant'
pdata[which(pdata$StromalP > 0.01),'StromalColor'] = 'NS'
pdata[which(pdata$ImmuneP > 0.01),'ImmuneColor'] = 'NS'
pdata$Fill = paste0('Stromal:',pdata$StromalColor,'\n','Immune:',pdata$ImmuneColor)
pdata$Label = NA
sigix = which(pdata$StromalColor == 'Significant' | pdata$ImmuneColor == 'Significant')
pdata[sigix,'Label'] = rownames(pdata)[sigix]

ggplot(pdata,aes(StromalScore,ImmuneScore,color = Fill,label = Label)) + geom_point() + theme_bw() +
  scale_color_manual(values = cvals,name = 'Statistical\nSignificance') + 
  geom_text_repel(show.legend = FALSE) +
  theme(legend.position = 'bottom')
  
