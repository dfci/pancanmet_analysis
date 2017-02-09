# Script to visualize cancer types by heterogeneity

rm(list = ls())
setwd('/Users/ereznik/pancanmet_analysis/analysis/')
source('plottingconventions.R')
library(ggplot2)
library(reshape2)

# Make empty theme
new_theme_empty <- theme_bw()
new_theme_empty$line <- element_blank()
new_theme_empty$rect <- element_blank()
new_theme_empty$strip.text <- element_blank()
new_theme_empty$axis.text <- element_blank()
new_theme_empty$plot.title <- element_blank()
new_theme_empty$axis.title <- element_blank()
new_theme_empty$plot.margin <- structure(c(0, 0, -1, -1), unit = "lines", valid.unit = 3L, class = "unit")

# Read in the signed changes file
d = read.csv('../results/diffabundance/DifferentialAbundance_SignedChanges.csv',header = TRUE,row.names= 1)
d$Metname = rownames(d)
d2 = melt(d,id.vars = 'Metname')
d2$Tissue = sourcetissue[as.character(d2$variable)]

# Order the studies
d2$variable = factor(d2$variable,levels = c('KIRC','BLCA','BRCA','BRCATang','OV','PAAD','PAADHussain1',
                                            'PAADHussain2','PRAD','PRADLODA'))

# Set tissue colors
tissuecols = c('Kidney' = 'gray','Bladder' = 'green','Breast' = 'blue','Ovary' = 'red',
               'Pancreas' = 'green','Prostate' = 'purple','Zero' = 'blue')

ccmmetabolites = c('glucose',
                   'glucose-6-phosphate (g6p)',
                   'fructose-6-phosphate',
                   'fructose-1,6-bisphosphate',
                   'dihydroxyacetone phosphate (dhap)',
                   'dihydroxyacetonephosphate (dhap)',
                   'glyceraldehyde-3-phosphate',
                   '3-phosphoglycerate',
                   '2-phosphoglycerate',
                   'phosphoenolpyruvate (pep)',
                   'lactate',
                   'citrate',
                   'cis-aconitate',
                   'succinate',
                   'fumarate',
                   'malate',
                   'alpha-ketoglutarate',
                   'glutamate',
                   'glutamine')

# For each metabolite in central carbon metabolism, make a barplot
for (metname in ccmmetabolites){
  dplot = d2[which(d2$Metname == metname),]
  zerodata = dplot[which(dplot$value == 0),]
  
  # Set NAs to zero in dplot so that they are plotted
  dplot[which(is.na(dplot$value)),'value'] = 0
  
  if(length(which(is.na(dplot$value))) == dim(dplot)[1]){
    print(paste0('Skipping ',metname))
    next
  }
  
  if (dim(zerodata)[1] > 0){
    
    ggplot(dplot,aes(variable,value)) + geom_bar(stat = 'identity',aes(fill = Tissue)) + new_theme_empty + 
      geom_bar(data = zerodata,stat = 'identity',aes(x = variable,y = .35,color = Tissue,fill = 'none'),size = 10) +
      geom_bar(data = zerodata,stat = 'identity',aes(x = variable,y = -.35, color = Tissue,fill = 'none'),size = 10) +
      scale_fill_manual(values = tissuecols) + theme(legend.position = 'none') + 
      scale_color_manual(values = tissuecols) + 
      geom_hline(yintercept = 0,size = 6) + ylim(-1,1) + 
      ggsave(paste0('../results/barplots_bymetabolite/',metname,'.pdf'))
  }else{
    ggplot(dplot,aes(variable,value,fill = Tissue)) + geom_bar(stat = 'identity') + new_theme_empty + 
      scale_fill_manual(values = tissuecols) + theme(legend.position = 'none') + 
      geom_hline(yintercept = 0,size = 6) + ylim(-1,1) + 
      ggsave(paste0('../results/barplots_bymetabolite/',metname,'.pdf'))
}
  
}