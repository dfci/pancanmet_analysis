# A script to analyze tumors in the PanCan study which do not exhibit the Warburg effect
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
pthresh = 0.05
textsize = 28

selectfirst <- function(s) strsplit(s, ":")[[1]][1]
selectlast <- function(s) strsplit(s, ":")[[1]][ length( strsplit(s, ":")[[1]] )]

pal <- function(col, border = "light gray", ...){
    n <- length(col)
    plot(0, 0, type="n", xlim = c(0, 1), ylim = c(0, 1),
         axes = FALSE, xlab = "", ylab = "", ...)
    rect(0:(n-1)/n, 0, 1:n/n, 1, col = col, border = border)}

# Import the fold change data
fc = fread('/Users/ereznik/Documents/pancanmet/data/merged_metabolomics/tnpaired_fc.csv',sep = ',',header = T,skip = 1)
fc = data.frame(fc,row.names = 1)
title.line = strsplit( readLines('/Users/ereznik/Documents/pancanmet/data/merged_metabolomics/tnpaired_fc.csv', n=1), ',' ) [[1]] # Note that first element is empty ""
colnames(fc) = title.line[2:length(title.line)]
splittypes = strsplit(colnames(fc),'\\:')
tissuetype = sapply(splittypes, "[", 1)

# Read in all the data
temp = fread('../../data/merged_metabolomics/alldata.csv',sep = ',',header = F,skip = 3)
alldata = data.frame( temp, row.names = 1 )
title.line = strsplit( readLines('../../data/merged_metabolomics/alldata.csv', n=1), ',' ) [[1]] # Note that first element is empty ""
colnames(alldata) = title.line[2:length(title.line)]

# Drop empty rowname
emptyix = which(rownames(alldata)=='')
if (length(emptyix) !=0){
  alldata = alldata[-emptyix,]
  fc = fc[-emptyix,]
}

# Extract relevant data from columns
study = unlist( lapply(strsplit(colnames(alldata),':'), `[[`, 1) )

# Focus in on a metabolite and plot it
hotmet = 'Lactate'
hotdata = t(fc[hotmet,])
idx = as.vector( order(hotdata) )
plotdata = data.frame(1:length(idx),log2(hotdata[idx]),tissuetype[idx])
colnames(plotdata)=c('x','y','tissue')
rownames(plotdata) = rownames(hotdata)[idx]
plotdata = plotdata[complete.cases(plotdata),]

# Make a waterfall plot
xinterceptidx = length(which(plotdata$y < 0))
ggplot(plotdata,aes(x,y,color = tissue)) + geom_point(aes(size =5))  +
  annotate("rect", xmin = 0, xmax = xinterceptidx, ymin = -Inf, ymax = 0, fill= "blue",alpha = 0.1) +
  annotate("rect", xmin = xinterceptidx, xmax = Inf, ymin = 0, ymax = Inf, fill= "red",alpha = 0.1) +
  xlab('') + ylab(paste(expression(Log[2]),'Tumor Lactate:Normal Lactate')) + theme_bw() + 
  scale_size_continuous(guide=FALSE) + scale_color_discrete(name = 'Study') + 
  theme( legend.position=c(0.75, .2),text = element_text(size = textsize) ) +
  geom_hline() + geom_vline(xintercept = xinterceptidx) + 
  ggsave(filename = '../../results/warburg_lactateplot.pdf',width = 10,height = 10)

# Make an additional column to indicate whether the tumor is Warburgian or not
plotdata$Warburg = '[Lactate] >2-Fold Higher in Tumor'
plotdata[which(plotdata$y< -1),'Warburg'] = '[Lactate] >2-Fold Lower in Tumor'
plotdata[which(abs(plotdata$y)< 1),'Warburg'] = 'No Substantial Change'

plotdata$UpDown = 'Lactate Up'
plotdata[which(plotdata$y<0),'UpDown'] = 'Lactate Down'

# Make a bar plot which indicates, for each study, what proportion of samples were "non-Warburgian"
warburgtable = table(plotdata$tissue,plotdata$Warburg)
warburgtable = warburgtable/rowSums(warburgtable)
warburgtable = data.frame(warburgtable)
colnames(warburgtable) = c('Study','WarburgEffect','Frequency')
warburgtable = warburgtable[complete.cases(warburgtable),]

ggplot(warburgtable,aes(Study,Frequency,fill = WarburgEffect)) + geom_bar(width = 0.5,stat = 'identity',position = 'dodge') + 
  theme_bw() + scale_fill_grey(name = 'Warburg Effect') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position=c(.3, .92)) + 
  theme(text = element_text(size = textsize) ) +
  scale_fill_manual(values=c('red','blue','grey'), name='') + xlab('') +
  ylab('Proportion') + ggsave('../../results/warburg_barplot.pdf',width = 10,height = 10)

updowntable = table(plotdata$tissue,plotdata$UpDown)
updowntable = updowntable/rowSums(updowntable)
updowntable = data.frame(updowntable)
colnames(updowntable) = c('Study','UpDown','Frequency')

ggplot(updowntable,aes(Study,Frequency,fill = UpDown)) + geom_bar(width = 0.5,stat = 'identity',position = 'dodge') + 
  theme_bw() + scale_fill_grey(name = 'Warburg Effect') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position=c(.3, .9)) + 
  theme(text = element_text(size = textsize) ) +
  ylab('Proportion') + ggsave('../../results/warburg_updown_barplot.pdf',width = 10,height = 10)

# Print how many tissues have less lactate than normal tissue
numlow = length(which(plotdata$y < 0))/dim(plotdata)[1]
print(paste('Proportion of samples with low lactate levels compared to normal:',numlow))

################################################################################################################
# Make a plot which compares the fold change of each paired tissue to the fold change versus the median normal tissue

# For each tissue type, compute the median lactate level
uqtis = as.character( unique( plotdata$tissue ) )
normallactate = data.frame( matrix(NA, length(uqtis),1) )
rownames(normallactate) = uqtis
for (tis in uqtis){
  regstr = paste(tis,'.*Normal',sep = '')
  ix = grep(regstr,colnames(alldata),perl = TRUE)
  normallactate[tis,1] = median( as.matrix( alldata['Lactate',ix] ) )
}

# For each sample, compare the tumor's lactate level to the median for a normal tissue
plotdata$vsMedianNormal = NA
for (pat in rownames(plotdata)){
  plotdata[pat,'vsMedianNormal'] = log2( alldata['Lactate',pat] / normallactate[as.character(plotdata[pat,'tissue']),1] )
}

# Plot
xmax = max(plotdata$y) + 0.5
xmin = min(plotdata$y) - 0.5
ymax = max(plotdata$vsMedianNormal) + 0.5
ymin = min(plotdata$vsMedianNormal) - 0.5
ggplot(plotdata,aes(y,vsMedianNormal,color = tissue)) + geom_point(aes(shape=tissue)) + theme_classic() + 
  xlab('Ratio, Lactate Tumor:Lactate Normal, Paired Tissue') + 
  ylab('Ratio, Lactate Tumor:Median of All Normal Lactate') + geom_hline() + geom_vline() + 
  scale_shape_manual(values=seq(0,length(uqtis))) + 
  #geom_segment(x = 1, y = 1, xend = 1, yend = ymax, color = 'gray') + 
  #geom_segment(x = 1, y = 1, xend = xmax, yend = 1, color = 'gray') + 
  #geom_segment(x = -1, y = -1, xend = -1, yend = ymin, color = 'gray') + 
  #geom_segment(x = -1, y = -1, xend = xmin, yend = -1, color = 'gray')  +
  facet_wrap(~tissue,nrow = 3) + geom_abline(linetype=2) + 
  theme(legend.position="none",text = element_text(size=textsize)) + 
  ggsave('../../results/warburg_scatter.pdf',height = 10,width = 10)
################################################################################################################

# Find the patients which show very decreases in hot metabolite, and check if this is due to abnormally high normal levels
lowpats = rownames(hotdata)[which(hotdata <= 1)]

# Get all of the alldata
tumlac = data.frame( t( alldata[hotmet,rownames(plotdata)] ) )
rownames(tumlac) = rownames(plotdata)

# Repeat for the normal data
tum2norm = read.csv( '/Users/ereznik/Documents/pancanmet/data/merged_metabolomics/tumor_normal_pairs4R.csv',header = F,row.names = 1)
normlac = data.frame( t( alldata[hotmet,as.character( tum2norm[rownames(tumlac),1 ] ) ] ) )

# Bind all data
plotdata2 = rbind(tumlac,normlac)
plotdata2$Study = sapply(rownames(plotdata2),selectfirst)
plotdata2$Type = sapply(rownames(plotdata2),selectlast)

# If there are any duplicates (ie 2 tumors for 1 normal), remove
plotdata2[which(plotdata2$Type=='Normal.1'),'Type'] = 'Normal'

# Indicate whether hot metabolite goes up or down
plotdata2$DeltaMet = paste('Log2', hotmet, 'Tumor:Normal > 0')
plotdata2[lowpats,'DeltaMet'] = paste('Log2', hotmet, 'Tumor:Normal < 0')
plotdata2[as.character(tum2norm[lowpats,1]),'DeltaMet'] = paste('Log2', hotmet, 'Tumor:Normal < 0')

# Get rid of bladder cancer, no hotmet
plotdata2 = plotdata2[-which(plotdata2$Study %in% c('BLCA','COAD','STAD')),] #### CORRECT THIS!!!!! WHEN HOTMET CHANGES, NOT CORRECT!

# Make faceted histogram
ggplot( plotdata2, aes_string(x = hotmet,fill = 'DeltaMet')) + geom_density(adjust = .5) + 
    facet_grid(Study ~ Type,scales = 'free') +
    theme_classic() + theme(text = element_text(size=textsize)) + 
    ggsave('../../results/warburg_source.pdf',width = 10,height = 10)

################################################################################################################

# For each study, compare the low lactate patients to the high lactate patients
highvslow = data.frame( matrix(NA,dim(alldata)[1],length(uqtis)) )
rownames(highvslow) = rownames(alldata)
colnames(highvslow) = uqtis
highvslowp = highvslow
highvslowpadj = highvslow
for (tis in uqtis){
  
  # Get the samples of this tissue type
  regstr = paste(tis,':.*Tumor',sep = '')
  ix = grep(regstr,colnames(alldata))
  d = alldata[,ix]
  d = d[which(complete.cases(d)),]
  d = as.matrix(d)
  
  # NOTE: some of the columns in d don't have paired normal tissue. We get around that using lesslac and morelac below
  
  # Make sure lactate is in this data
  if (!(hotmet %in% rownames(d))){
    print(paste(hotmet),'not in study',tis,',dropping from analysis...')
    highvslow = highvslow[,-which(colnames(highvslow) == tis)] 
    highvslowp = highvslowp[,-which(colnames(highvslowp) == tis)]
    highvslowpadj = highvslowpadj[,-which(colnames(highvslowpadj) == tis)]
    next
  
  }
  
  # Find the samples which have less lactate than matched normal tissue
  lesslac = grep(regstr, colnames(fc)[ which(fc[hotmet,] < 1) ], perl = TRUE, value = TRUE )
  morelac = grep(regstr, colnames(fc)[ which(fc[hotmet,] > 1) ], perl = TRUE, value = TRUE )
  
  if (length(lesslac) < 2){
    print(paste('Inadequate samples of ',hotmet,'in study',tis,',dropping from analysis...'))
    highvslow = highvslow[,-which(colnames(highvslow) == tis)] 
    highvslowp = highvslowp[,-which(colnames(highvslowp) == tis)]
    highvslowpadj = highvslowpadj[,-which(colnames(highvslowpadj) == tis)]
    next
  }
  if (length(morelac) < 2){
    print(paste('Inadequate samples of ',hotmet,'in study',tis,',dropping from analysis...'))
    highvslow = highvslow[,-which(colnames(highvslow) == tis)] 
    highvslowp = highvslowp[,-which(colnames(highvslowp) == tis)]
    highvslowpadj = highvslowpadj[,-which(colnames(highvslowpadj) == tis)]
    next
  }
  
  # Do a diff abundance test
  tempp = data.frame( matrix(NA, dim(d)[1],1) )
  rownames(tempp) = rownames(d)
  for (i in 1:dim(d)[1]){
    mname = rownames(d)[i]
    highvslow[mname,tis] = log2( mean(d[mname,morelac])/mean(d[mname,lesslac]) )
    tempp[mname,1] = wilcox.test( d[mname,morelac], d[mname,lesslac] )$p.value
  }
  
  # We need to correct the p-values in a smart way, because they are stored in an array that is larger than it should be
  temppcor = p.adjust(tempp[,1],method= 'BH')
  temppcor = data.frame(temppcor)
  rownames(temppcor) = rownames(d)
  
  # Correct
  highvslowp[rownames(tempp),tis] = tempp
  highvslowpadj[rownames(temppcor),tis] = temppcor
  highvslow[which(highvslowpadj[,tis] > pthresh),tis] = 0
}

highvslowsums = rowSums(sign(abs(highvslow)))
sigidx = which(abs(highvslowsums) > 1)
pheatmap(highvslow[sigidx,])
