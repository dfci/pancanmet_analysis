# Check that tumor/normal assignments are correct. 
# An easy negative control is just changing any of the studies so that normals get labeled as tumors

rm(list = ls())
setwd('/Users/ereznik/pancanmet_analysis/checks/')
library(ggplot2)
library(XLConnect)

logfile = 'Excel2CSVcheck_tumornormal.txt'
write('Beginning to check importing of data...\n',logfile,append = FALSE)

comptypes = function(study,alltypes,excel){
  
  diffix = which(excel$Type != alltypes[rownames(excel),'Type'])
  sameix = which(excel$Type == alltypes[rownames(excel),'Type'])
  
  if (length(diffix)!=0){
    stop('Error in matching!')
  }else{
    write(paste(study,'is ok!'),logfile,append = TRUE)
    write(paste('Number of samples that match:',length(sameix)),logfile,append = TRUE)
    write(paste('Number of samples in big data:',length(which(alltypes$Study == study))),logfile,append = TRUE)
  }
}

########################################################################################
# Read in the big data file
bigd = read.csv('../data/merged_metabolomics/merged_metabolomics.csv',header = TRUE,row.names = 1, check.names = FALSE)
alltypes = data.frame('FullName' = colnames(bigd), stringsAsFactors = FALSE)
alltypes$Short = sapply(alltypes$FullName,function(x){
  y = strsplit(x,'\\:')[[1]];
  z = paste(y[1:2],collapse = ':')
})
alltypes$Type = sapply(alltypes$FullName,function(x){
  y = strsplit(x,'\\:')[[1]][3];
})
alltypes$Study = sapply(alltypes$FullName,function(x){
  y = strsplit(x,'\\:')[[1]][1];
})
alltypes$JustID = sapply(alltypes$FullName,function(x){
  y = strsplit(x,'\\:')[[1]][2];
})
rownames(alltypes) = alltypes$Short

########################################################################################
# 1. Bladder
# Read the excel file
BLCAwb = loadWorkbook('../data/studies/BLCA/PutlurietalData_profiling (1).xlsx')
BLCA = readWorksheet(BLCAwb,'Profiling_allcpds_processed',header = FALSE, check.names = FALSE,
                     rownames = 1,endRow = 2)
# Make values
BLCA_types = data.frame(t(BLCA['Cancer Status',]))
rownames(BLCA_types) = paste('BLCA:',BLCA['CPD corrected',],sep = '')
BLCA_types$Type = NA
BLCA_types[which(BLCA_types$Cancer.Status == 'Cancer'),'Type'] = 'Tumor'
BLCA_types[which(BLCA_types$Cancer.Status == 'Benign'),'Type'] = 'Normal'

# Compare
comptypes('BLCA',alltypes,BLCA_types)

########################################################################################
# 2. Breast Terunama
BRCAwb = loadWorkbook('../data/studies/BRCA/BRCA_alldata.xlsx')
BRCA = readWorksheet(BRCAwb,'ScaledImpData',header = FALSE, check.names = FALSE,startRow = 2,endRow = 3,
                     startCol = 10,rownames = 1)

# Make values
BRCA_types = data.frame(t(BRCA['DIAG',]))
rownames(BRCA_types) = paste('BRCA:',BRCA['SAMPLE_ID',],sep = '')
BRCA_types$Type = NA
BRCA_types[which(BRCA_types == 'TUMOR'),'Type'] = 'Tumor'
BRCA_types[which(BRCA_types == 'NORMAL'),'Type'] = 'Normal'

# Compare
comptypes('BRCA',alltypes,BRCA_types)

########################################################################################
#3. BRCATang
# This one is simple, all tumors have TCGA-style barcode, and all normals are of length 2 and have a leading N
BTtypes = alltypes[which(alltypes$Study == 'BRCATang'),]
for (row in rownames(BTtypes)){
  if (grepl('\\-',BTtypes[row,'Short']) & BTtypes[row,'Type'] == 'Tumor'){
    write(paste('This should be a tumor:',BTtypes[row,'FullName']),logfile,append = TRUE)
  }else if (substr(BTtypes[row,'JustID'],1,1) == 'N' & BTtypes[row,'Type'] == 'Normal'){
    write(paste('This should be a normal:',BTtypes[row,'FullName']),logfile,append = TRUE)
  }else{
    stop()
  }
}
write('Breast Tang is OK!',logfile,append = TRUE)

########################################################################################
#4. KIRC
KIRCwb = loadWorkbook('../data/studies/KIRC/Supp Table 2_MedianNormalized.xlsx')
KIRC = readWorksheet(KIRCwb,sheet = 'Median Normalized',header = TRUE,check.names = FALSE,
                     startRow = 2, endRow = 3,startCol = 12, rownames = 1)
KIRC_types = data.frame( t(KIRC), stringsAsFactors = FALSE )
rownames(KIRC_types) = paste('KIRC:',rownames(KIRC_types),sep = '')
KIRC_types$Type = NA
KIRC_types[which(KIRC_types$TISSUE.TYPE == 'TUMOR'),'Type'] = 'Tumor'
KIRC_types[which(KIRC_types$TISSUE.TYPE == 'NORMAL'),'Type'] = 'Normal'
comptypes('KIRC',alltypes,KIRC_types)

########################################################################################
#5. PAADHussain1
PH1wb = loadWorkbook('../data/studies/PAADHussain1/NCIA-08-10VW CDT (final).xlsx')
PH1 = readWorksheet(PH1wb,sheet = 'ScaledImpData',header = TRUE,check.names = FALSE,
                    startRow = 2, startCol = 11, endRow = 6,rownames = 1)
PH1 = data.frame( t(PH1), stringsAsFactors = FALSE )
rownames(PH1) = paste('PAADHussain1:',rownames(PH1),sep = '')
PH1$Type = NA
PH1[which(PH1$TISSUE == 'Tumor'),'Type'] = 'Tumor'
PH1[which(PH1$GROUP_ID == 'Normal_donor'),'Type'] = 'Donor'
PH1[which(PH1$GROUP_ID == 'Nontumor'),'Type'] = 'Normal'
comptypes('PAADHussain1',alltypes,PH1)

########################################################################################
#6. PAADHussain2
PH2wb = loadWorkbook('../data/studies/PAADHussain2/NCIA-20-11VW CDT 14Feb2012.xlsx')
PH2 = readWorksheet(PH2wb,sheet = 'ScaledImpData',header = TRUE,check.names = FALSE,
                    startRow = 1, endRow = 9,startCol = 12, rownames = 1)
PH2 = data.frame( t(PH2), stringsAsFactors = FALSE)
rownames(PH2) = paste('PAADHussain2:',rownames(PH2),sep = '')
PH2$Type = PH2$GROUP_ID
comptypes('PAADHussain2',alltypes,PH2)

########################################################################################
#7. OV
OVwb = loadWorkbook('../data/studies/OV/LOUI-01-10VW CDT d1-22July2010.xlsx')
OV = readWorksheet(OVwb,sheet = 'ScaledImpData',header = TRUE,check.names = FALSE,
                   startRow = 2,endRow = 9,startCol = 11,rownames = 1)
OV = data.frame( t(OV), stringsAsFactors = FALSE)
rownames(OV) = paste('OV:',rownames(OV),sep = '')
OV$Type = NA
OV[which(OV$GROUP_DESC == 'Normal'),'Type'] = 'Normal'
OV[which(OV$GROUP_DESC == 'EOC'),'Type'] = 'Tumor'
OV[which(OV$GROUP_DESC == 'Metastatic EOC'),'Type'] = 'Metastasis'
comptypes('OV',alltypes,OV)

########################################################################################
# 8. PAAD
# This one is again easy, just check whether tumor or normal is in sample name
PAADwb = loadWorkbook('../data/studies/PAAD/KamphorstNofal_Metabolomics_Data.xlsx')
PAAD = readWorksheet(PAADwb,sheet = 'Processed_Data',header = TRUE,check.names = FALSE,
                     startRow = 1,endRow = 1,startCol = 1,rownames = 1)
PAAD = data.frame(t(PAAD))
rownames(PAAD) = paste('PAAD:',rownames(PAAD),sep = '')
PAAD$Type = NA
PAAD[grep('Tumor',rownames(PAAD)),'Type'] = 'Tumor'
PAAD[grep('Normal',rownames(PAAD)),'Type'] = 'Normal'
comptypes('PAAD',alltypes,PAAD)

########################################################################################
#9. PRADLODA
PLwb = loadWorkbook('../data/studies/PRADLODA/131853_1_supp_2658512_nbpsn6.xlsx')
PL = readWorksheet(PLwb,sheet = 'Raw data Human tumors',header = TRUE,check.names = FALSE,
                   startRow = 1,startCol = 7,rownames = 1,endRow = 2)
PL = data.frame( t(PL), stringsAsFactors = FALSE)
rownames(PL) = paste('PRADLODA:',rownames(PL),sep = '')
PL$Type = NA
PL[which(PL$TUMOR_NORMAL == 'T'),'Type'] = 'Tumor'
PL[which(PL$TUMOR_NORMAL == 'N'),'Type'] = 'Normal'
comptypes('PRADLODA',alltypes,PL)

########################################################################################
#10. LGG. Everything should be a tumor, easy!
LGG = alltypes[which(alltypes$Study == 'LGG'),]
if (length(which(LGG$Type == 'Tumor'))!= dim(LGG)[1]){
  stop()
}else{
  write('LGG is OK, all tumors!',logfile,append = TRUE)
}

########################################################################################
#11. PRAD

# From the Readme file, we know:
#N1 - N16	Benign prostate sample data
#T1 - T12	Localized tumor sample data
#M1 - M14	Metastatic tumor sample data

# Make sure this matches
PRAD = alltypes[which(alltypes$Study == 'PRAD'),]
PRAD$FirstCharacter = sapply(PRAD$JustID,function(x){substr(x,1,1)})
if (length(which(PRAD$FirstCharacter == 'T' & PRAD$Type != 'Tumor'))!=0){
  stop()
}
if (length(which(PRAD$FirstCharacter == 'N' & PRAD$Type != 'Normal'))!=0){
  stop()
}
if (length(which(PRAD$FirstCharacter == 'M' & PRAD$Type != 'Metastasis'))!=0){
  stop()
}
write('PRAD is OK!',logfile,append = TRUE)