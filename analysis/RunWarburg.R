# Script to compare levels of lactate in paired tumor/normal tissues

rm(list=ls())
setwd('/Users/ereznik/pancanmet_analysis/analysis/')
library(ggplot2)
library(ggrepel)
source('functions/readBigMet.R')

# Choose the metabolite to focus on, by default should be lactate
hotmet = 'lactate'

# Read in the paired metabolomics data
d = read.csv('../data/merged_metabolomics/tnpaired_fc.csv',header = TRUE,row.names = 1,check.names = FALSE)
splittypes = strsplit(colnames(d),'\\:')
studytype = sapply(splittypes, "[", 1)
names(studytype) = colnames(d)
tissuemap = read.csv('../data/merged_metabolomics/merged_tissuetypes.csv',header = TRUE,row.names = 1)

pdata = data.frame( t(d[hotmet,]) )
splittypes = sapply(rownames(pdata),function(x){ strsplit(x,'\\:') })
pdata$Study = sapply(splittypes,function(x){x[[1]]})
pdata$Tissue = tissuemap[pdata$Study,1]
colnames(pdata) = c('Metabolite','Study','Tissue')
pdata = pdata[order(pdata$Metabolite),]
pdata$Metabolite = log2(pdata$Metabolite)
pdata$IX = (1:dim(pdata)[1])/dim(pdata)[1]
pdata$Elevation = NA
pdata[which(pdata$Metabolite>0),'Elevation'] = 'Elevated'
pdata[which(pdata$Metabolite<0),'Elevation'] = 'Depleted'

# Make a waterfall plot
xinterceptidx = length(which(pdata$Metabolite < 0))/dim(pdata)[1]
ggplot(pdata,aes(IX,Metabolite,color = Tissue)) + geom_point(size =3)  +
  annotate("rect", xmin = 0, xmax = xinterceptidx, ymin = -Inf, ymax = 0, fill= "blue",alpha = 0.1) +
  annotate("rect", xmin = xinterceptidx, xmax = Inf, ymin = 0, ymax = Inf, fill= "red",alpha = 0.1) +
  xlab('Percentile') + ylab(paste(expression(Log2),'Tumor Lactate:Normal Lactate')) + theme_bw() + 
  scale_size_continuous(guide=FALSE) + scale_color_discrete(name = 'Study') + 
  theme( legend.position=c(0.75, .2),text = element_text(size = 20) ) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = xinterceptidx) 
  ggsave(filename = '../results/warburg/warburg_lactateplot.pdf',width = 10,height = 10,useDingbats = FALSE)

# For all metabolites, go through and calculate the proportion of "elevated" values
elevated = apply(d,1,function(x){ length(which(x>1))/length(x[!is.na(x)] ) } )
uqelevated = apply(d,1,function(x){length(unique(studytype[which(!is.na(x))]))})
elev = data.frame('Proportion' = elevated, 'Num' = uqelevated)
elev = elev[elev$Num > 6,]
elev = elev[order(elev$Proportion),]
elev$Name = ''
sigelev = which(elev$Proportion < 0.3 | elev$Proportion > 0.7)
elev[sigelev,'Name'] = rownames(elev)[sigelev]
elev$IX = (1:dim(elev)[1])/dim(elev)[1]

# Plot
ggplot(elev,aes(IX,Proportion,label = Name)) + geom_point() + 
  xlab('Quantile') + ylab('Proportion of Samples with\nHigher Abundance in Tumor') + 
  geom_text_repel(data = subset(elev,elev$Name!=''),angle = 0,size = 6,
                  box.padding = unit(0.1,'lines'),point.padding = unit(0.1,'lines')) + 
  theme_bw() + theme(text = element_text(size=40)) +
  ggsave('../results/warburg/allmetaboltes_paired.pdf',height = 10,width = 10,useDingbats = FALSE)
  