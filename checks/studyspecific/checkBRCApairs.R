# A script to confirm that we have the correct pairs of breast tumors/normals

rm(list = ls())
setwd('/Users/ereznik/pancanmet_analysis/checks/studyspecific/')
library(XLConnect)

# The data is formatted really weirdly but we can still work with it
sidata = readWorksheetFromFile('../checkdata/BRCApairs_fromSI.xlsx',sheet = 'Sheet1')
sidata = sidata[complete.cases(sidata),]
rownames(sidata) = sidata$Metabolon.Tumor.ID.

tnpairs = read.csv('../../data/merged_metabolomics/tumor_normal_pairs.csv',header = TRUE)
tnpairs = tnpairs[tnpairs$Study == 'BRCA',]
rownames(tnpairs) = tnpairs$Tumor.Sample

for (rname in rownames(tnpairs)){
  # Get the pair in the sidata and confirm its the same
  tnvalue = as.character( tnpairs[rname,'Normal.Sample'] )
  sivalue = as.character( sidata[rname,'Metabolon.Normal.ID.'] )
  
  if (tnvalue != sivalue){
    stop('Error in matching!')
  }else{
    print(c(tnvalue,sivalue,'Perfect Match'))
  }
  
}

# Also, read in the BRCA data and confirm that the age of the samples is the same
metdata = readWorksheetFromFile('../../data/studies/BRCA/BRCA_alldata.xlsx',sheet = 'OrigData',
                                startRow = 2, startCol = 10,endRow = 18,check.names = FALSE,
                                header = TRUE,rownames = 1)

metdata = apply(metdata,2,as.character)
for (rname in rownames(tnpairs)){
  # Check that all the data for them is the same
  tumsample = rname
  normsample = as.character( tnpairs[rname,'Normal.Sample'] )
  
  tdata = metdata[,tumsample]
  ndata = metdata[,normsample]
  
  diffix = which(tdata != ndata)
  print(c(tumsample,normsample,tdata[diffix],ndata[diffix]))
  
}
