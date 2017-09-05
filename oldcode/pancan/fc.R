# Script to calculate fold changes

rm(list=ls())
setwd('/Users/ereznik/Documents/pancanmet/analysis/pancan/')
library(pheatmap)
library(NMF)
library(RColorBrewer)
library(ggplot2)
library(reshape)
library(colorspace)
library(scales)

calcdist = 1
maxval = 3
pthresh = 0.05
textsize = 20

br = seq( -maxval, maxval,.1)
brb = length(which(br<0)); bra = length(which(br>0))
cl = colorRampPalette(c("blue","white","red"))(length(br)-1)

alldata = read.csv('/Users/ereznik/Documents/pancanmet/data/merged_metabolomics/alldata.csv', row.names = 1,check.names = F)
met = data.frame( sapply(alldata[3:dim(alldata)[1],], function(x) { as.numeric(as.character(x)) } ), check.names = F )
rownames(met) = rownames(alldata)[3:dim(alldata)[1]]
splittypes = strsplit(colnames(alldata),'\\:')
tissuetype = sapply(splittypes, "[", 1)
uqtis = unique(tissuetype)

# Remove LGG, no normals
uqtis = setdiff(uqtis,'LGG')

# Read in the tissue type
tmap = read.csv('../../data/merged_metabolomics/merged_tissuetypes.csv', row.names =1, header = T)
lggidx = which(rownames(tmap) == 'LGG')
rnamestmap = rownames(tmap)
tmap = data.frame( tmap[-lggidx,] )
rownames(tmap) = rnamestmap[-lggidx]

tnvalue = sapply(splittypes, "[", 3)

# For each tumor type and metabolite, calculate fold changes
fc = matrix(NA,dim(met)[1],length(uqtis)); rownames(fc) = rownames(met); colnames(fc) = uqtis;
p = matrix(1,dim(met)[1],length(uqtis)); rownames(p) = rownames(met); colnames(p) = uqtis;

for (j in 1:length(uqtis)){
  print(uqtis[j])
  
  tumidx = which(tissuetype == uqtis[j] & tnvalue == 'Tumor')
  normidx = which(tissuetype == uqtis[j] & tnvalue == 'Normal')
  
  if (length(normidx) < 2 | length(tumidx) < 2){next}
  
  for (i in 1:dim(met)[1]){
    
    if (is.na(met[i,tumidx[1]])){next}
    
    p[i,j] = wilcox.test( as.numeric(met[i,tumidx]) , as.numeric(met[i,normidx]) )$p.value
    fc[i,j] = log2( mean(as.numeric(met[i,tumidx])) / mean(as.numeric(met[i,normidx]))  )
    
  }
  
}

# Correct p-values
padj = p
for (i in 1:dim(p)[2]){
  padj[,i] = p.adjust(p[,i],method = 'BH')
}

fcadj = fc
fcadj[which(padj > pthresh & !is.na(fc),arr.ind = T)] = 0

# Before writing to a file, lets calculate some statistics
metstats = data.frame( matrix(0,dim(met)[1],6) )
rownames(metstats) = rownames(met); 
colnames(metstats) = c('NumUp','PctUp','NumDown','PctDown','TotalDiff','notNA');
    
for (i in 1:dim(met)[1]){
  
  notNA = length(which(!is.na(fcadj[i,])))
  
  metstats[i,'notNA'] = notNA
  
  metstats[i,'NumUp'] = length(which(fcadj[i,]>0))
  metstats[i,'NumDown'] = length(which(fcadj[i,]<0))
  
  metstats[i,'TotalDiff'] = metstats[i,'NumUp'] + metstats[i,'NumDown']
  
  metstats[i,'PctUp'] = metstats[i,'NumUp']/notNA
  metstats[i,'PctDown'] = metstats[i,'NumDown']/notNA
  
}

# Make a table showing the frequency of metabolites changing versus their measurement frequency
freqtable = data.frame( table(metstats[,c('notNA','TotalDiff')]) )
freqtable[-which(freqtable$Freq == 0),]
ggplot(freqtable,aes(notNA,TotalDiff)) + geom_point(aes(size= log10(Freq))) + theme_classic(20) + theme(text = element_text(size=textsize))

# Make a quantile plot
quantileplot = data.frame(Diff = metstats$NumUp - metstats$NumDown, TotalDiff = metstats$TotalDiff, 
                          notNA = metstats$notNA, pct = metstats$PctUp + metstats$PctDown)

# Reorder things
quantileplot = quantileplot[order(quantileplot[,'Diff'],quantileplot[,'notNA']),]
quantileplot$CDF = rank(quantileplot$Diff,ties.method = 'first')/dim(quantileplot)[1]
ggplot(quantileplot,aes(CDF,Diff,size = notNA,color = pct)) + geom_point() + theme_classic(20) +
  scale_color_gradient2(low = 'green',mid = 'yellow',high='red',midpoint = 0.5,name = '% Differentially Abundant') + xlab('Quantile') + 
  scale_size_continuous(name = '# of Times Measured') +
  ylab('# Studies Increased - # Studies Decreased') +
  scale_x_continuous(breaks=pretty_breaks(10)) +
  scale_y_continuous(breaks=pretty_breaks(10)) +
  ggsave('/Users/ereznik/Documents/pancanmet/results/QuantilePlot.pdf',height = 10,width =10,useDingbats = FALSE)

# Bind all data together
colnames(fc) = lapply(colnames(fc),function(x){paste('FC',x)})
colnames(p) = lapply(colnames(p),function(x){paste('P',x)})
colnames(padj) = lapply(colnames(padj),function(x){paste('Padj',x)})
writedata = cbind(fc,padj,metstats)
write.csv(writedata,'../../results/DifferentialAbundanceSummary.csv')

# Make a histogram indicating the frequency of differential abundance across studies
ggplot(metstats,aes(NumUp + NumDown)) + geom_histogram(binwidth=0.5, origin = -0.3) + theme_classic() +  # Change origin value to align bins
  xlab('Number of Studies with Significant Differential Abundance') + 
  ylab('Number of Metabolites') +
  scale_x_continuous(breaks=c(0:max(metstats$NumUp + metstats$NumDown + 1)) ) + 
  theme(text = element_text(size=textsize)) + 
  ggsave('/Users/ereznik/Documents/pancanmet/results/Histogram_NumDiffAbundant.pdf',height = 10,width =10)

# Make an extra fc dataframe for plotting
fcadj = fc
fcadj[which(p > pthresh,arr.ind = T)] = 0

##########################################################################################
# Make bar plots to indicate the relative frequency of imputation and differential abundance
##########################################################################################

# Make a bar plot that indicates the number of metabolites tested and the number with sig fold changes for all studies
bardata = data.frame(matrix(0,length(uqtis),8),row.names = uqtis)
colnames(bardata) = c('Total','SigDiff','SigUp','SigDown','TumorImputed','NormalImputed','NumTumors','NumNormals')
for (j in 1:length(uqtis)){
  tumidx = which(tissuetype == uqtis[j] & tnvalue == 'Tumor')
  normidx = which(tissuetype == uqtis[j] & tnvalue == 'Normal')
  tempdata = met[,c(tumidx,normidx)]
  completedata = tempdata[complete.cases(tempdata),]
  
  # Number of metabolites measured
  bardata[uqtis[j],'Total'] = dim(completedata)[1]
  
  # Number of differentially abundant metabolites
  bardata[uqtis[j],'SigDiff'] = length(which(abs(fcadj[, paste('FC',uqtis[j]) ])>1))
  
  # Number of differentially increased metabolites
  bardata[uqtis[j],'SigUp'] = length(which(fcadj[,paste('FC',uqtis[j])] > 1))
  
  # Number of differentially decreased metabolites
  bardata[uqtis[j],'SigDown'] = length(which(fcadj[,paste('FC',uqtis[j])]< -1))
  
  # Number of metabolites with >20% imputation for tumors
  rowmins = apply(completedata,1,min)
  completeplus = cbind(completedata,rowmins)
  tumimput = apply( completeplus,1,function(x){length( which( abs(x[1:length(tumidx)] - x[length(x)] ) < 1e-5) ) } ) # Check only number of imputed TUMORS
  bardata[uqtis[j],'TumorImputed'] = length(which(tumimput > length(tumidx)*0.2))
  
  normimput = apply( completeplus,1,function(x){length( which( abs( x[(length(tumidx)+1):(length(x)-1)] - x[length(x)] ) < 1e-5) ) } ) # Check only number of imputed NORMALS
  bardata[uqtis[j],'NormalImputed'] = length(which(normimput > length(normidx)*0.2))
  
  # Number of tumor and normal samples
  bardata[uqtis[j],'NumTumors'] = length(tumidx)
  bardata[uqtis[j],'NumNormals'] = length(normidx)
}

bardata$Study = rownames(bardata)
bardatamelt = melt(bardata,id.vars = 'Study')

ggplot(bardatamelt,aes(Study,value,fill = variable)) + geom_bar(width = 0.7,stat = 'identity',position = 'dodge') + 
  theme_minimal() +  theme(text = element_text(size=textsize)) + 
  scale_fill_brewer(palette='Dark2',name = '',labels = colnames(bardata)) +
  theme(axis.text.x = element_text(angle = 45, hjust = .5),legend.position=c(.8, .8),text = element_text(size=textsize)) + 
  ylab('Proportion')

# Repeat bar plot using proportions 
barplot = bardata
barplot[,2:6] = barplot[,2:6]/barplot[,1] 
barplot = barplot[,-c(1,2,7,8)]
barplotmelt = melt(barplot,id.vars = 'Study')

ggplot(barplotmelt,aes(Study,value,fill = variable)) + geom_bar(stat = 'identity',position = 'dodge',width = 0.6,alpha = 0.8) + 
  theme_minimal() +  theme(text = element_text(size=textsize)) + 
  scale_fill_brewer(palette='Dark2',name = '',labels = 
                      c('% Metabolites Increased','% Metabolites Decreased','% Metabolites Imputed,Tumor','% Metabolites Imputed, Normal')) +
  theme(axis.text.x = element_text(angle = 45, hjust = .5),legend.position=c(.85, .85),text = element_text(size=textsize)) + 
  ylab('Proportion')

# Repeat one more time but removing imputation data
barplot = bardata
barplot[,2:6] = barplot[,2:6]/barplot[,1] 
barplot = barplot[,-c(1,2,5,6,7,8)]
barplotmelt = melt(barplot,id.vars = 'Study')

ggplot(barplotmelt,aes(Study,value,fill = variable)) + geom_bar(stat = 'identity',position = 'dodge',width = 0.6,alpha = 0.8) + 
  theme_minimal() +  theme(text = element_text(size=textsize)) +
  scale_fill_brewer(palette='Dark2',name = '',labels = 
                      c('% Metabolites Increased','% Metabolites Decreased')) +
  theme(axis.text.x = element_text(angle = 45, hjust = .5),legend.position=c(.85, .85),text = element_text(size=textsize)) + 
  ylab('Proportion of Metabolites') + ggsave('../../results/FC_barplot.pdf',height = 10,width = 10)

# One more time using just imputed data
imputplot = bardata
imputplot[,2:6] = imputplot[,2:6]/imputplot[,1]
imputplot = imputplot[,c(5,6,9)]
imputplotmelt = melt(imputplot,id.vars = 'Study')
ggplot(imputplotmelt,aes(Study,value,fill = variable)) + geom_bar(stat = 'identity',position = 'dodge',width = 0.6,alpha = 0.8) + 
  theme_minimal() + scale_fill_brewer(palette='Dark2',name = '',labels = 
                      c('% Metabolites Imputed, Tumor','% Metabolites Imputed, Normal')) +
  theme(axis.text.x = element_text(angle = 45, hjust = .5),legend.position=c(.85, .85),text = element_text(size=textsize)) + 
  ylab('Proportion of Metabolites') + ggsave('../../results/imput_barplot.pdf',height = 10,width = 10)

##########################################################################################
# Make a heatmap of metabolites found at least n times
##########################################################################################
nzidx = rownames(metstats)[which((metstats[,'NumUp'] + metstats[,'NumDown']) > 5)]

plotdata = fcadj[nzidx,]
colnames(plotdata) = unlist(lapply(colnames(plotdata),function(x){y = strsplit(x,' ');y[[1]][[2]]}))
diffscore = (rowSums(sign(plotdata)))/dim(plotdata)[2]
plotdata = plotdata[order(diffscore),]
pheatmap(plotdata,color = cl,breaks = br,cluster_col = FALSE,annotation_col = tmap,
         filename = '../../results/FC_heatmap.pdf',height = 10,width = 10)

# Get data for a cool stacked bar plot
meltmets = melt(sign(plotdata))

# Re-order the levels of the metabolite names
metord = rank(diffscore,ties.method = 'first')
meltmets$Xpos = metord[as.character( meltmets$X1 ) ]

ggplot(meltmets,aes(Xpos,value, fill = X2)) + geom_bar(stat= 'identity',alpha = 0.8) + 
  theme_minimal() + coord_flip() + geom_hline()
