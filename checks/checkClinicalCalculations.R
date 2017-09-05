# Script to check clinical associations

rm(list = ls())
setwd('/Users/ereznik/pancanmet_analysis/checks/')
library(ggplot2)
library(cowplot)
library(XLConnect)
library(metap)

# Determine whether we want to use stage or grade
f2use = 'grade'

# Read in big data file too, to get names for merged mapping
bigmet = read.csv('../data/merged_metabolomics/merged_metabolomics.csv', row.names = 1,check.names = FALSE,header = TRUE)

# Read in clinical dta
clin = read.csv('../data/merged_metabolomics/clinfeatures.csv',row.names =  1, stringsAsFactors = FALSE)

# Read in Irina's results
irina = readWorksheetFromFile('../results/clinical/updated results stage and grade Jan 2017 max1.xlsx',sheet = f2use,
                              header = TRUE,rownames = 1,startRow = 2)
irina_cnames = colnames(irina)
irina_rnames = rownames(irina)
irina = apply(irina,2,as.numeric)
rownames(irina) = irina_rnames
colnames(irina) = irina_cnames

# For each study in clin, if it has grade information, correlate
uqstudy = unique(clin$study)

# Make a table for all results
resexist = FALSE

for (study in uqstudy){
  print(study)
  
  clinstudy = clin[clin$study == study,]
  clinstudy = clinstudy[,c('study','type',f2use)]
  clinstudy = clinstudy[complete.cases(clinstudy),]
  
  if (dim(clinstudy)[1] == 0){
    print(paste('Skipping this study for no samples with grade info:',study))
    next
  }
  
  # Rename the rownames
  clinstudy$sample = rownames(clinstudy)
  rownames(clinstudy) = paste(study,clinstudy$sample,clinstudy$type,sep = ':')
  
  # Calculate correlation with each metabolite separately in tumor and normal samples
  tumsamples = rownames(clinstudy)[clinstudy$type=='Tumor']
  tumsamples = intersect(tumsamples,colnames(bigmet)) # Added Jan 2017
  if (length(which(clinstudy$type == 'Normal'))>0){
    normsamples = rownames(clinstudy)[clinstudy$type=='Normal']
    normsamples = intersect(normsamples,colnames(bigmet)) # Added Jan 2017
  }else{
    normsamples = c()
  }
  
  
  # Restrict to metabolite data
  smallmet = bigmet[,c(tumsamples,normsamples)]
  smallmet = smallmet[complete.cases(smallmet),]
  smallmet = as.matrix(smallmet)
  
  # for each metabolite, calculate p-value of correlation
  res = data.frame('PositivePTumor' = numeric(0),'NegativePTumor' = numeric(0),
                   'PositivePNormal' = numeric(0),'NegativePNormal' = numeric(0) )
  
  for (metname in rownames(smallmet)){
    res[metname,'PositivePTumor'] = cor.test(clinstudy[tumsamples,f2use],smallmet[metname,tumsamples],
                                             method = 'spearman',alternative = 'greater')$p.value  
    
    res[metname,'NegativePTumor'] = cor.test(clinstudy[tumsamples,f2use],smallmet[metname,tumsamples],
                                             method = 'spearman',alternative = 'less')$p.value  
  }
  
  if (length(normsamples)>0){
    
    for (metname in rownames(smallmet)){
      res[metname,'PositivePNormal'] = cor.test(clinstudy[normsamples,f2use],smallmet[metname,normsamples],
                                               method = 'spearman',alternative = 'greater')$p.value  
      
      res[metname,'NegativePNormal'] = cor.test(clinstudy[normsamples,f2use],smallmet[metname,normsamples],
                                               method = 'spearman',alternative = 'less')$p.value  
    }
    
  }
  
  # Plot the results against Irina's findings
  ixmets = intersect(rownames(irina),rownames(res))
  irina_positivecor = paste(study,'.Tumor.less',sep = '')
  irina_negativecor = paste(study,'.Tumor.greater',sep = '')
  p1 = qplot(irina[ixmets,irina_positivecor],res[ixmets,'PositivePTumor']) + scale_x_log10() + scale_y_log10() + 
    ggtitle(study) + geom_abline()
  p2 = qplot(irina[ixmets,irina_negativecor],res[ixmets,'NegativePTumor']) + scale_x_log10() + scale_y_log10() + 
    ggtitle(study) + geom_abline()
  
  print(p1)
  print(p2)
  
  # Save the results to allres
  if (!(resexist)){
    allres = data.frame(res[ixmets,c('PositivePTumor','NegativePTumor')])
    colnames(allres) = c(irina_positivecor,irina_negativecor)
    resexist = TRUE
  }else{
    tempdf = data.frame(res[ixmets,c('PositivePTumor','NegativePTumor')])
    colnames(tempdf) = c(irina_positivecor,irina_negativecor)
    allres = merge(allres,tempdf,by = 'row.names',all = TRUE)
    rownames(allres) = allres$Row.names
    
    # Drop the row.names column
    allres = allres[,which(colnames(allres) != 'Row.names')]
  }
  
  if (length(normsamples)>0){
    
    irina_positivecornormal = paste(study,'.Normal.less',sep = '')
    irina_negativecornormal = paste(study,'.Normal.greater',sep = '')
    p3 = qplot(irina[ixmets,irina_positivecornormal],res[ixmets,'PositivePNormal']) + scale_x_log10() + scale_y_log10() + 
      ggtitle(paste(study,'Normal')) + geom_abline()
    p4 = qplot(irina[ixmets,irina_negativecornormal],res[ixmets,'NegativePNormal']) + scale_x_log10() + scale_y_log10() + 
      ggtitle(paste(study,'Normal')) + geom_abline()
    
    print(p3)
    print(p4)
    
    tempdf = data.frame(res[ixmets,c('PositivePNormal','NegativePNormal')])
    colnames(tempdf) = c(irina_positivecornormal,irina_negativecornormal)
    allres = merge(allres,tempdf,by = 'row.names',all = TRUE)
    rownames(allres) = allres$Row.names
    
    # Drop the row.names column
    allres = allres[,which(colnames(allres) != 'Row.names')]
    
  }  
  
}

# Combine the p-values for tumors and for all samples
negcolumns = grep('greater',colnames(allres))
negtumorcolumns = grep('Tumor.greater',colnames(allres))
poscolumns =  grep('less',colnames(allres))
postumorcolumns = grep('Tumor.less',colnames(allres))
allres$all.combined = NA
allres$tumors.combined = NA 

for (metname in rownames(allres)){
  
  for (c2use in c('all.combined','tumors.combined')){
    
    if (c2use == 'all.combined'){
      posc2use = poscolumns
      negc2use = negcolumns
    }else if(c2use == 'tumors.combined'){
      posc2use = postumorcolumns
      negc2use = negtumorcolumns
    }else{
      stop('Error in choosing columns to combine')
    }
    
    # Calculate one-sided p values for all data
    posvalues = na.omit( t(allres[metname,posc2use] ) )
    if (length(posvalues) > 1){
      posp = sumlog( na.omit( posvalues ) )$p
    }else if (length(posvalues) == 1){
      posp = posvalues  
    }else{
      posp = NA
    }
    
    negvalues = na.omit( t(allres[metname,negc2use] ) )
    if (length(negvalues) > 1){
      negp = sumlog( na.omit( negvalues ) )$p
    }else if (length(negvalues) == 1){
      negp = negvalues  
    }else{
      negp = NA
    }
    
    if (is.na(posp) & is.na(negp)){
      allres[metname,c2use] = NA
    }else{
      allres[metname,c2use] = 2*min(posp,negp,na.rm = TRUE)
    } 
    
  }
  
}

# Adjust p-values 
allres$all.combined.adjusted = p.adjust(allres$all.combined,method = 'BH')
allres$tumors.combined.adjusted = p.adjust(allres$tumors.combined,method = 'BH')

# Compare to Irina
ixmetsall = intersect(rownames(irina),rownames(allres))
allcompare = qplot(irina[ixmetsall,'all.combined.adjusted'],allres[ixmetsall,'all.combined.adjusted']) + 
  scale_x_log10() + scale_y_log10() + ggtitle('AllData') + geom_abline()

alltumor = qplot(irina[ixmetsall,'tumors.combined.adjusted'],allres[ixmetsall,'tumors.combined.adjusted']) + 
  scale_x_log10() + scale_y_log10() + ggtitle('TumorData') + geom_abline()

print(allcompare)
print(alltumor)

alldiff = abs( log10(irina[ixmetsall,'all.combined.adjusted']) - log10(allres[ixmetsall,'all.combined.adjusted']) )
