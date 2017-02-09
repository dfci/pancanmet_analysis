# Visualize differential abundance results

rm(list = ls())
setwd('/Users/ereznik/pancanmet_analysis/analysis/')
source('../analysis/plottingconventions.R')
source('../import/useful_metimport.R')
library(ggplot2)
library(reshape2)

d = read.csv('../results/diffabundance/DifferentialAbundanceSummary.csv',header = TRUE,row.names = 1)
dsign = read.csv('../results/diffabundance/DifferentialAbundance_SignedChanges.csv',header = TRUE,row.names = 1)
d$FracDiff = d$TotalDiff/d$notNA
d = d[order(d$FracDiff,decreasing = TRUE),]
d$Index = 1:dim(d)[1]

######
# Plot 1

# Make inset plot
ggplot(d,aes(Index,FracDiff)) + geom_point() + theme_bw()

######
# Plot 2
histdata = as.data.frame( table(d$TotalDiff) )
ggplot(histdata,aes(Var1,Freq)) + geom_bar(stat='identity') + xlab('Frequency of Differential Abundance') + 
  ylab('Number of Metabolites') + ggsave('../results/diffabundance/diffabundancehistogram.pdf',height = 8, width = 8)



######
# Plot 3

# Make main plot
rsums = rowSums(abs(dsign),na.rm = TRUE)
sigmets = names(rsums)[which(rsums > 5)]
rsums2 = rsums[sigmets]
dsign2 = dsign[sigmets,]
dsign2$Name = rownames(dsign2)
dsign2 = melt(dsign2,id.vars = 'Name')
dsign2$value = factor(dsign2$value)
dsign2$Name = sapply(dsign2$Name,simpleCap)

# Make nice annotation, first capitalize metclasses
names(metclasses) = sapply(names(metclasses),simpleCap)
metclasses = metclasses[order(metclasses)]
metorder = names(metclasses)

dsign2$Name2Plot = names2plot[as.character(dsign2$variable)]
dsign2$MetClass = factor(metclasses[as.character(dsign2$Name)])
dsign2$MetName = factor(dsign2$Name,levels = metorder)

dsign2$Tissue = sourcetissue[as.character(dsign2$variable)]
dsign2$ValueNoFactor = as.numeric(as.character(dsign2$value))

print(unique(as.character(dsign2[which(is.na(dsign2$MetClass)),'Name'])))


# ggplot(dsign2,aes(Name2Plot,MetName,shape = value, color = value)) + geom_point(size = 4,alpha = 0.75) + 
#   theme_bw() + 
#   scale_shape_manual(name = 'Direction of Change',
#                      values = c('-1' = '\u25BC','0' = '\u25CF','1' = '\u25B2')) + 
#   scale_color_manual(values = c('-1' = 'blue','0' = 'gray','1' = 'red')) + 
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
#   facet_wrap(~Tissue,nrow = 1,scale = 'free_x')

ggplot(dsign2,aes(Name2Plot,MetName,color = value)) + geom_point(size = 4,alpha = 0.75) + 
  theme_bw() + geom_point()+
  scale_color_manual(values = c('-1' = 'blue','0' = 'gray','1' = 'red')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.text.y = element_text(angle = 45)) + 
  facet_grid(~Tissue,scale = 'free_x', space = 'free_x') + 
  xlab('') + ylab('') + theme(legend.position = 'none') + 
  ggsave('../results/diffabundance/DifferentialAbundance.pdf',height = 8,width = 8)

