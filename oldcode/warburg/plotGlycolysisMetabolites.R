# Glycolysis plot

rm(list = ls())
library(NMF)
library(RColorBrewer)
library(ggplot2)
library(reshape)
library(colorspace)
library(data.table)
library(pheatmap)
cmap = brewer.pal(12,'Set3')
setwd('/Users/ereznik/Documents/pancanmet/analysis/warburg/')

# Set p-value threshold
pthresh = 0.05

# Read in the metabolites
glyc = read.csv('GlycolysisMetabolites.csv',header = 0)
glyc[,1] = as.character(glyc[,1])

# Read the metabolite results
results = read.csv('../../results/DifferentialAbundanceSummary.csv',header = 1,row.names = 1)
results = results[glyc[,1],]

studies_temp = colnames(results)[ grep('FC',colnames(results)) ]
studies = unlist(lapply(studies_temp,function(x){strsplit(x,'\\.')[[1]][2]}))

# Make a dataframe for plotting
res = data.frame(matrix(NA,length(glyc[,1]),length(studies)))
colnames(res) = studies
rownames(res) = glyc[,1]

for (study in studies){
  for (met in glyc[,1]){
    
    tempsign = sign( results[met,paste('FC',study,sep='.')] )
    tempp = results[met,paste('Padj',study,sep='.')]
   
    if(is.na(tempsign)){next}
    
    if (tempp > pthresh){
      res[met,study] = 0
    } else{
      res[met,study] = tempsign
    }
  }
  
}

# Plot a pie chartres
resmelt = melt(as.matrix(res))
colnames(resmelt) = c('Metabolite','Study','Value')
ggplot(resmelt,aes(Study,Value,fill = Value)) + geom_bar(stat = 'identity') + coord_polar("y", start=0) + facet_wrap(~Metabolite) + 
  theme_classic() + scale_color_discrete(na.value = "grey50")
