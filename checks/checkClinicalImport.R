# Script to compare newly imported clinical data to old clinical data

rm(list = ls())
analysisDir = normalizePath(".")
shinyDir = '../pancanmet/'
setwd('/Users/ereznik/pancanmet_analysis/checks/')
#setwd(file.path(analysisDir, "checks"))

met = read.csv('../data/merged_metabolomics/merged_metabolomics.csv',header = TRUE,row.names = 1,check.names = FALSE)

logfile = 'ClinicalCheck.txt'
write('Starting to compare clinical files...',logfile,append = FALSE)

# Read in old and new clinical data
#old = read.csv('checkdata/clinfeatures_DEPRECATED_April2016.csv',header = TRUE)
old = read.csv('checkdata/clinfeatures_DeprecatedNov2016.csv',header = TRUE)
colnames(old)[1:5] = c('Sample','Study','Type','Stage','Grade') # capitalization of column names matters

new = read.csv('../data/merged_metabolomics/clinfeatures.csv', header = TRUE)

# Make good rownames
rownames(old) = paste(old$Study,old$Sample,old$Type,sep = ':')
rownames(new) = paste(new$study,new$sample,new$type,sep = ':')

# Make sure that patients in metabolomics data match those in clinical data
diffpatsmet = setdiff(rownames(new),colnames(met))
if (length(diffpatsmet)!=0){
  print('Error in patient names in clinical/metabolomics file.')
  print(diffpatsmet)
}

ixpats = intersect(rownames(old),rownames(new))
diffpats = c( setdiff(rownames(old),rownames(new)), setdiff(rownames(new),rownames(old)) )

# For each study, make a table of the intersecting grades and stages
ixold = old[ixpats,]
ixnew = new[ixpats,]

ixold$Stage = as.character(ixold$Stage)
ixold$Grade = as.character(ixold$Grade)

ixnew$stage = as.character(ixnew$stage)
ixnew$grade = as.character(ixnew$grade)

for (study in unique(ixold$Study)){
  
  studypats = rownames(ixold)[which(ixold$Study == study) ]
  stagetab = table( ixold[studypats,'Stage'], ixnew[studypats,'stage'] )
  gradetab = table( ixold[studypats,'Grade'], ixnew[studypats,'grade'] )
  
  #write('---------',logfile,append = TRUE)
  #write(paste('Working on',study),logfile,append = TRUE)
  #write(stagetab,logfile,append = TRUE)
  #write(gradetab,logfile,append = TRUE)
  
  print(study)
  print(stagetab)
  print(gradetab)
  print('Number of mismatching stages:')
  print(length(which(ixold[studypats,'Stage']!=ixnew[studypats,'stage'])))
  print('Number of mismatching grades:')
  print(length(which(ixold[studypats,'Grade']!=ixnew[studypats,'grade'])))
}

setwd(analysisDir)
