# Script to compare the association between kidney metabolites and clinical data using Irina's calculations and ours from Cancer Cell paper

rm(list = ls())
setwd('/Users/ereznik/pancanmet_analysis/checks/')
library(ggplot2)
library(cowplot)
library(XLConnect)
library(metap)

# Determine whether we want to use stage or grade
f2use = 'stage'

# Read in Irina's results
irina = readWorksheetFromFile('../results/clinical/updated results stage and grade Jan 2017 max1.xlsx',
                              sheet = f2use,
                              header = TRUE,rownames = 1,startRow = 2)
irinak = irina[,c("KIRC.Tumor.less","KIRC.Tumor.greater")]
irinak[,1] = as.numeric(irinak[,1])
irinak[,2] = as.numeric(irinak[,2])
irinak[which(irinak == 'NA',arr.ind = TRUE)] = NA
irinak = irinak[which(complete.cases(irinak)),]

# Read in big data file too, to get names for merged mapping
bigmet = read.csv('../data/merged_metabolomics/merged_metabolomics.csv', 
                  stringsAsFactors = FALSE,row.names = 1,check.names = FALSE,header = TRUE)


# Read in the mapping to the original metabolites
mapping = read.csv('../import/tempdir/merged_mapping.csv',header = TRUE,row.names = 1,stringsAsFactors = FALSE)
rownames(mapping) = rownames(bigmet)
rownames(irinak) = mapping[rownames(irinak),'KIRC']

# Set significant and insignificant
irinak$Change = 'None'
irinak[which(irinak$KIRC.Tumor.less < 5e-2),'Change'] = 'Higher in High Stage'
irinak[which(irinak$KIRC.Tumor.greater < 5e-2),'Change'] = 'Higher in Low Stage'

# Read in data from Cancer Cell paper
ccell = read.csv('/Users/ereznik/Documents/ccrc/savedresults/ccrc_Utest_stage.csv',header = TRUE,row.names= 1)
rownames(ccell) = sapply(rownames(ccell),tolower)
ccell$Change = 'None'
ccell[which(ccell$Adjusted.P.Value < 5e-2 & ccell$Log2.FC > 0),'Change'] = 'CCell:Higher in High Stage'
ccell[which(ccell$Adjusted.P.Value < 5e-2 & ccell$Log2.FC < 0),'Change'] = 'CCell:Higher in Low Stage'

irinak$CCell = ccell[rownames(irinak),'Change']
print( table(irinak$Change,irinak$CCell))
