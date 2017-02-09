# Script to calculate fold changes

rm(list=ls())
setwd('/Users/ereznik/pancanmet_analysis/analysis/')
source('plottingconventions.R')
library(RColorBrewer)
library(ggplot2)
library(reshape)
library(colorspace)
library(scales)
library(cowplot)

calcdist = 1
maxval = 3
textsize = 20

br = seq( -maxval, maxval,.1)
brb = length(which(br<0)); bra = length(which(br>0))
cl = colorRampPalette(c("blue","white","red"))(length(br)-1)

met = read.csv('../data/merged_metabolomics/merged_metabolomics.csv', row.names = 1,check.names = FALSE,header = TRUE)
splittypes = strsplit(colnames(met),'\\:')
tissuetype = sapply(splittypes, "[", 1)
tumnormtype = sapply(splittypes, "[", 3)
uqtis = unique(tissuetype)

# Remove LGG, no normals
uqtis = setdiff(uqtis,'LGG')

# Read in the tissue type
tmap = read.csv('../data/merged_metabolomics/merged_tissuetypes.csv', row.names =1, header = T)
lggidx = which(rownames(tmap) == 'LGG')
rnamestmap = rownames(tmap)
tmap = data.frame( tmap[-lggidx,] )
rownames(tmap) = rnamestmap[-lggidx]

# For each tumor type and metabolite, calculate fold changes
fc = matrix(NA,dim(met)[1],length(uqtis)); rownames(fc) = rownames(met); colnames(fc) = uqtis;
p = matrix(1,dim(met)[1],length(uqtis)); rownames(p) = rownames(met); colnames(p) = uqtis;

for (j in 1:length(uqtis)){
  print(uqtis[j])
  
  tumidx = which(tissuetype == uqtis[j] & (tumnormtype == 'Tumor' | tumnormtype == 'Metastasis') ) 
  normidx = which(tissuetype == uqtis[j] & tumnormtype == 'Normal')
  
  if (length(normidx) < 2 | length(tumidx) < 2){next}
  
  for (i in 1:dim(met)[1]){
    
    if ( any(is.na(met[i,tumidx])) ){
      p[i,j] = NA
      fc[i,j] = NA
    }else{
      p[i,j] = wilcox.test( as.numeric(met[i,tumidx]) , as.numeric(met[i,normidx]) )$p.value
      fc[i,j] = log2( mean(as.numeric(met[i,tumidx])) / mean(as.numeric(met[i,normidx]))  ) 
    }
    
  }
  
}

# Correct p-values
padj = p
for (i in 1:dim(p)[2]){
  padj[,i] = p.adjust(p[,i],method = 'BH')
}

fcadj = fc
fcadj[which(padj > pthresh & !is.na(fc),arr.ind = TRUE)] = 0

# Also make a matrix just indicating whether metabolites go up/down in each study, this will be useful for pathway analysis
fcsign = sign( fcadj )

# Before writing to a file, lets calculate some statistics
metstats = data.frame( matrix(0,dim(met)[1],11) )
rownames(metstats) = rownames(met); 
colnames(metstats) = c('SigStudies','SigTissues','NumUp','PctUp','NumDown','PctDown',
                       'PctUpFC','PctDownFC','TotalDiff','TotalDiffFC','notNA');
    
for (i in 1:dim(met)[1]){ #for each metabolite
  
  notNA = length(which(!is.na(fcadj[i,])))
  
  sigix = which(fcadj[i,]!=0)
  metstats[i,'SigStudies'] = paste(colnames(fc)[sigix],collapse = ';')
  metstats[i,'SigTissues'] = paste( unique(sourcetissue[colnames(fc)[sigix]]), collapse = ';')
  metstats[i,'notNA'] = notNA
  
  metstats[i,'NumUp'] = length(which(fcadj[i,]>0))
  metstats[i,'NumDown'] = length(which(fcadj[i,]<0))
  
  metstats[i,'TotalDiff'] = metstats[i,'NumUp'] + metstats[i,'NumDown']
  
  metstats[i,'PctUp'] = metstats[i,'NumUp']/notNA
  metstats[i,'PctDown'] = metstats[i,'NumDown']/notNA
  
  # Repeat, now considering which metabolites have at least a 2-fold change
  metstats[i,'PctUpFC'] = length(which(fcadj[i,] > 1))/notNA
  metstats[i,'PctDownFC'] = length(which(fcadj[i,] < -1))/notNA
  metstats[i,'TotalDiffFC'] = length(which(fcadj[i,] > 1)) + length(which(fcadj[i,] < -1))
  
}

# Bind all data together
colnames(fc) = lapply(colnames(fc),function(x){paste('FC',x)})
colnames(p) = lapply(colnames(p),function(x){paste('P',x)})
colnames(padj) = lapply(colnames(padj),function(x){paste('Padj',x)})
writedata = cbind(metstats,fc,padj)
write.csv(writedata,'../results/diffabundance/DifferentialAbundanceSummary.csv')

# Write the sign of each metabolite change to a file as well
write.csv(fcsign,'../results/diffabundance/DifferentialAbundance_SignedChanges.csv')

##########################################################################################
# Make bar plots to indicate the relative frequency of imputation and differential abundance
##########################################################################################

# Make a bar plot that indicates the number of metabolites tested and the number with sig fold changes for all studies
bardata = data.frame(matrix(0,length(uqtis),10),row.names = uqtis)
colnames(bardata) = c('Total','SigDiff','SigUp','SigDown','SigUpFC',
                      'SigDownFC','TumorImputed','NormalImputed','NumTumors','NumNormals')

for (j in 1:length(uqtis)){
  tis2use = uqtis[j]
  
  tumidx = which(tissuetype == tis2use & tumnormtype == 'Tumor')
  normidx = which(tissuetype == tis2use & tumnormtype == 'Normal')
  tempdata = met[,c(tumidx,normidx)]
  completedata = tempdata[complete.cases(tempdata),]
  
  # Number of metabolites measured
  bardata[uqtis[j],'Total'] = dim(completedata)[1]
  
  # Number of differentially abundant metabolites
  bardata[uqtis[j],'SigDiff'] = length(which( fcadj[,tis2use]!= 0 & !is.na(fcadj[,tis2use]) ))
  
  # Number of differentially increased metabolites
  bardata[uqtis[j],'SigUp'] = length(which(fcadj[,tis2use] > 0))
  
  # Number of differentially decreased metabolites
  bardata[uqtis[j],'SigDown'] = length(which(fcadj[,tis2use] < 0))
  
  # Number of metabolites with >20% imputation for tumors
  rowmins = apply(completedata,1,min)
  completeplus = cbind(completedata,rowmins)
  
  tumimput = apply( completeplus,1,function(x){length( which( abs(x[1:length(tumidx)] - x[length(x)] ) < 1e-5) ) } ) # Check only number of imputed TUMORS
  bardata[uqtis[j],'TumorImputed'] = length(which(tumimput > length(tumidx)*0.2))
  
  normstart = length(tumidx)+1
  normend = length(tumidx) + length(normidx)
  normimput = apply( completeplus,1,function(x){length( which( abs( x[normstart:normend] - x[length(x)] ) < 1e-5) ) } ) # Check only number of imputed NORMALS
  bardata[uqtis[j],'NormalImputed'] = length(which(normimput > length(normidx)*0.2))
  
  # Count number of metabolites with abolute fold change > 2
  bardata[uqtis[j],'SigUpFC'] = length(which(fcadj[,tis2use] > 1))
  bardata[uqtis[j],'SigDownFC'] = length(which(fcadj[,tis2use] < -1))
  
  # Number of tumor and normal samples
  bardata[uqtis[j],'NumTumors'] = length(tumidx)
  bardata[uqtis[j],'NumNormals'] = length(normidx)
}

bardata$Study = names2plot[ rownames(bardata) ]
bardatamelt = melt(bardata,id.vars = 'Study')

ggplot(bardatamelt,aes(Study,value,fill = variable)) + geom_bar(width = 0.7,stat = 'identity',position = 'dodge') + 
  theme_minimal() +  theme(text = element_text(size=textsize)) + 
  scale_fill_brewer(palette='Spectral',name = '',labels = colnames(bardata)) +
  theme(axis.text.x = element_text(angle = 45, hjust = .5),legend.position=c(.8, .8),text = element_text(size=textsize)) + 
  ylab('Number of Metabolites') + ggsave('../results/barplots/alldata_barplot.pdf',height = 10,width = 10)


########################################################################################################
# Start making special plots

# Repeat bar plot using proportions 
cols2use = c("SigDiff","SigUp","SigDown","SigUpFC","SigDownFC","TumorImputed","NormalImputed")
barplot = bardata[,cols2use]/bardata[,1] 
barplot$Study = names2plot[ as.character(rownames(barplot)) ]
barplotmelt = melt(barplot,id.vars = 'Study')
barplotmelt$variable = factor(as.character(barplotmelt$variable), levels = cols2use)

ggplot(barplotmelt,aes(Study,value,fill = variable)) + geom_bar(stat = 'identity',position = 'dodge',width = 0.6,alpha = 0.8) + 
  theme_minimal() +  theme(text = element_text(size=textsize)) + 
  scale_fill_brewer(palette='Dark2',name = '',labels = 
                      c('% Metabolites Changed','% Metabolites Increased','% Metabolites Decreased',
                        '% Metabolites Increased > 2-fold','% Metabolites Decreased > 2-fold',
                        'Metabolites Imputed > 20%,Tumor','Metabolites Imputed > 20%, Normal')) +
  theme(legend.text=element_text(size=10)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = .5),legend.position=c(.85, .85),text = element_text(size=textsize)) + 
  ylab('Proportion') + ggsave('../results/barplots/alldata_proportions_barplot.pdf',height = 10,width = 10)

# Repeat one more time but removing imputation data
cols2use = c("SigDiff","SigUp","SigDown","SigUpFC","SigDownFC")
barplot = bardata[,cols2use]/bardata[,1] 
barplot$Study = names2plot[ as.character(rownames(barplot)) ]
barplotmelt = melt(barplot,id.vars = 'Study')
barplotmelt$variable = factor(as.character(barplotmelt$variable), levels = cols2use)

ggplot(barplotmelt,aes(Study,value,fill = variable)) + geom_bar(stat = 'identity',position = 'dodge',width = 0.6,alpha = 0.8) + 
  theme_minimal() +  theme(text = element_text(size=textsize)) + ylab('Proportion') +
  scale_fill_brewer(palette='Dark2',name = '',labels = 
                      c('% Metabolites Different',
                        '% Metabolites Increased','% Metabolites Decreased',
                        '% Metabolites Increased > 2-fold','% Metabolites Decreased > 2-fold')) +
  theme(axis.text.x = element_text(angle = 45, hjust = .5),legend.position=c(.85, .85),text = element_text(size=textsize)) + 
  ggsave('../results/barplots/FC_proportions_barplot.pdf',height = 10,width = 10)

# One more time using just imputed data
cols2use = c("TumorImputed","NormalImputed")
barplot = bardata[,cols2use]/bardata[,1] 
barplot$Study = names2plot[ as.character(rownames(barplot)) ]
barplotmelt = melt(barplot,id.vars = 'Study')
barplotmelt$variable = factor(as.character(barplotmelt$variable), levels = cols2use)
ggplot(barplotmelt,aes(Study,value,fill = variable)) + geom_bar(stat = 'identity',position = 'dodge',width = 0.6,alpha = 0.8) + 
  theme_minimal() + scale_fill_brewer(palette='Dark2',name = '',labels = 
                                        c('% Metabolites Imputed, Tumor','% Metabolites Imputed, Normal')) +
  theme(axis.text.x = element_text(angle = 45, hjust = .2),legend.position=c(.85, .85),text = element_text(size=textsize)) + 
  ylab('Proportion of Metabolites') + ggsave('../results/barplots/imput_barplot.pdf',height = 10,width = 10)

# One more time using just the differences with fold greater than 2
cols2use = c("SigUpFC","SigDownFC")
barplot = bardata[,cols2use]/bardata[,1] 
barplot$Study = names2plot[ as.character(rownames(barplot)) ]
barplotmelt = melt(barplot,id.vars = 'Study')
barplotmelt$variable = factor(as.character(barplotmelt$variable), levels = cols2use)
ggplot(barplotmelt,aes(Study,value,fill = variable)) + geom_bar(stat = 'identity',position = 'dodge',width = 0.6,alpha = 0.8) + 
  theme_minimal() + scale_fill_brewer(palette='Dark2',name = '',labels = 
                                        c('% Metabolites Increased > 2-fold','% Metabolites Decreased > 2-fold')) +
  theme(axis.text.x = element_text(angle = 45, hjust = .2),legend.position=c(.85, .85),text = element_text(size=textsize)) + 
  ylab('Proportion of Metabolites') + ggsave('../results/barplots/fc2_barplot.pdf',height = 10,width = 10)

# Make a plot of the proportion of differentially abundant metabolites versus sample size
barplot = bardata
barplot[,1:(dim(bardata)[2]-3)] = barplot[,1:(dim(bardata)[2]-3)]/bardata[,1]
barplot$SourceTissue = sourcetissue[ as.character(rownames(barplot)) ]
p1 = ggplot(barplot,aes(NumTumors + NumNormals,SigDiff,color = SourceTissue,label = Study)) + geom_point() + 
  scale_color_discrete(name = 'Cancer Type') + 
 geom_text(vjust="inward",hjust="inward",show.legend = FALSE) + xlab('') + ylab('\n') + ggtitle('BH-Corrected P-Value < 0.05')

p2 = ggplot(barplot,aes(NumTumors + NumNormals,SigUpFC + SigDownFC,color = SourceTissue,label = Study)) + geom_point() + 
  geom_text(vjust="inward",hjust="inward",show.legend = FALSE) +
  scale_color_discrete(name = 'Cancer Type') + 
  xlab('Number of Samples (Tumor + Normal)') + ylab('Proportion of Metabolites\nDifferentially Abundant') + 
  ggtitle('BH-Corrected P-Value < 0.05 and > 2-Fold Change in Abundance')

p1p2 = plot_grid(p1,p2,ncol = 1,labels = 'AUTO')
save_plot('../results/diffabundance/samples_vs_diffabundance.pdf',p1p2,base_height = 8,base_width = 8)

