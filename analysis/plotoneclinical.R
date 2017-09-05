# A script to visualize a clinical association and plot it

setwd('/Users/ereznik/pancanmet_analysis/analysis/')
library(ggplot2)
library(data.table)
source('plottingconventions.R')

# Remove all files from results from file
do.call(file.remove,list(list.files("../results/clinical/plots/",full.names = TRUE)))

# Open file
res = read.csv('../results/clinical/Visualized_grade.csv',header = TRUE, row.names = 1)
for (met in rownames(res)){
  print(met)
  
  pall = res[met,'all.combined.adjusted']
  ptumor = res[met,'tumors.combined.adjusted']
  
  if (pall > 5e-2 & ptumor > 5e-2){next} # skip insignificant results
  
  # Choose the metabolite and study (or studies) you would like to plot
  #met = 'asparagine' # for testing
  studysplit = unlist( strsplit( as.character(res[met,'SigNames'] ),';') )
  studies = sapply( studysplit, function(x){strsplit(x,'\\.')[[1]][1]})
  studies = unique(studies)
  
  # Read in big data and clinical data
  temp = fread('../data/merged_metabolomics/merged_metabolomics.csv',data.table = FALSE,header = TRUE)
  d = temp[,2:dim(temp)[2]]
  rownames(d) = temp[,1]
  clin = read.csv('../data/merged_metabolomics/clinfeatures.csv',header = TRUE)
  rownames(clin) = paste(clin[,2],clin[,1],clin[,3],sep = ':')
  
  # Trim the data
  clin = clin[clin$study %in% studies,]
  ixpats = intersect(rownames(clin),colnames(d))
  clin = clin[ixpats,]
  d = d[met,ixpats]
  
  # Merge the data
  pdata  = clin
  pdata[,'Metabolite'] = t(d)
  colnames(pdata) = c(colnames(clin),'Metabolite')
  pdata$grade = factor(pdata$grade)
  
  # Drop any NA's
  naix = which(is.na(pdata$grade))
  if (length(naix)!=0){
    pdata = pdata[-naix,]
  }
  
  # Put everything in the appropriate order
  xlevels = c( paste('Normal: Grade',1:4), paste('Tumor: Grade',1:4) )
  pdata$XLevels = paste(pdata$type,': Grade ',pdata$grade,sep = '')
  pdata$XLevels = factor(pdata$XLevels,levels = xlevels)
  
  # Convert y to log scale
  pdata$Metabolite = log2(pdata$Metabolite)
  
  # Convert study names
  pdata$StudyName = names2plot[as.character(pdata$study)]
  
  # Decide on number of columns
  if (length(unique(pdata$study)) <= 3){
    #numcols = 1;width = 5;
    numcols = 2;width = 10; # make everything 2x2
  }else{numcols = 2;width = 10}
  
  # Plot
  ggplot(pdata,aes(XLevels,Metabolite,fill = type)) + geom_boxplot(outlier.shape = NA) + theme_bw() + 
    geom_jitter() + 
    xlab('Grade') + ylab(paste('Log2',met)) + 
    facet_wrap(~StudyName,ncol = numcols,scale = 'free_y') + theme(legend.position = 'None') + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    theme(strip.background = element_rect(colour = 'white', fill = 'white', size = 3, linetype='dashed')) +
    scale_fill_manual(values = c('Normal' = 'blue','Tumor' = 'red')) + 
    ggsave(paste('../results/clinical/plots/',met,'.pdf',sep = ''),height = 10,width = width, useDingbats = FALSE)
    
}