# A script which checks that we correctly munged the metabolomics data

rm(list = ls())
options( java.parameters = "-Xmx6g" )
setwd('/Users/ereznik/pancanmet_analysis/checks/')
library(XLConnect)
library(data.table)

# Set a logfile
logfile = 'Excel2CSVcheck.txt'
write('Beginning to check importing of data...\n',logfile,append = FALSE)

# Issues:
# 1. Is there a reason the PRADLODA rownames are re-arranged? 
# 1-Resolved: August 5,2016: Not sure why it was re-arranged, but everything seems to check out.

readmunged = function( study ){
  # Function to read in the munged data that we analyze, accounting for our weird formatting
  f = paste('../data/metabolomics/',study,'/',study,'_metdata.csv',sep = '')
  
  # Read the data we use for analysis, note that first column will be the rownames
  temp = fread(f,sep = ',',header = FALSE,data.table = FALSE,skip = 2)
  tempanalyzed = data.frame( temp, row.names = 1 )
  
  # Get column names and tissue types, note that first element is empty (i.e. "")
  tempcolnames = strsplit( readLines(f, n=1), ',' )[[1]]
  colnames(tempanalyzed) = tempcolnames[2:length(tempcolnames)]
  
  # Get tissue types as well
  temptissue = strsplit( readLines(f, n=2), ',' )[[2]]
  temptissue = temptissue[2:length(temptissue)]
  
  reslist = list('analyzed' = tempanalyzed,'tissue' = temptissue)
  
  
}

compareDF = function(study,df1, df2, logfile){
  write(paste('\nChecking',study,'...'),logfile,append = TRUE)
  
  # Function which checks that rownames match, colnames match, and values match for two dataframes
  sortrnames1 = sort(rownames(df1))
  sortrnames2 = sort(rownames(df2))
  diffrnames = which( sortrnames1 != sortrnames2 )
  
  sortcnames1 = sort(colnames(df1))
  sortcnames2 = sort(colnames(df2))
    
  diffcnames = which( sortcnames1 != sortcnames2)
  if (length(diffrnames)!=0){
    write( paste(sortrnames1[diffrnames],sortrnames2[diffrnames]),logfile,append = TRUE )
    write( paste(study,'This rowname in munged data does not match.'),logfile,append = TRUE )
  }
  
  if (length(diffcnames)!=0){
    write( paste(sortcnames1[diffcnames], sortcnames2[diffcnames]),logfile,append = TRUE )
    write( paste(study,'This colname in munged data does not match'),logfile,append = TRUE )
  }
  diffvals = df1 - df2[rownames(df1),colnames(df1)]
  diffvals[which(is.na(diffvals),arr.ind= TRUE)] = 0
  maxdiff = max(diffvals)
  if (max( diffvals ) > 1e-6){stop(paste(study,'Values in munged and re-imported data do not match.'))}
  write(paste('Maximum difference in values between XLS and CSV in this study:',maxdiff),logfile,append = TRUE)
  
  return(maxdiff)
}

# 1. Bladder: email from Ali Shojaie:"The data is log-transformed and normalized (centered) according 
# to the procedure explained in our paper, which is a bit more involved than simple median centering.
# Also, Ali said: Yes, it is log transformed. But please keep in mind that the data has been normalized 
# AFTER log transformation.
BLCAwb = loadWorkbook('../data/studies/BLCA/PutlurietalData_profiling (1).xlsx')
BLCA = readWorksheet(BLCAwb,'Profiling_allcpds_processed',header = FALSE, check.names = FALSE,rownames = 1)

# The second row contains sample names, top row is tumor/normal. We only use rows up to 85, 
# other metabolites are not sufficiently identified
BLCAmunge = as.matrix(BLCA[3:85,] )
class(BLCAmunge) = 'numeric'
rownames(BLCAmunge) = rownames(BLCA)[3:85]
colnames(BLCAmunge) = BLCA[2,]

# Untransform from log2
BLCAmunge = 2^BLCAmunge
BLCAmunge = as.data.frame(BLCAmunge)

readBLCA = readmunged('BLCA')
BLCAanalyzed = readBLCA$analyzed
BLCAtissue = readBLCA$tissue

# Compare data frames
resBLCA = compareDF( 'BLCA', BLCAmunge, BLCAanalyzed, logfile )

# Compare tissues
mungedtissue = BLCA[1,]
mungedtissue[ which(mungedtissue == 'Benign') ]  = 'Normal'
mungedtissue[ which(mungedtissue == 'Cancer') ]  = 'Tumor'
BLCAtissuecheck = which(mungedtissue != BLCAtissue)

########################################################################################

# 2. BRCA
#BRCAwb = loadWorkbook('../data/studies/BRCA/BRCA_alldata.xlsx')
BRCAwb = loadWorkbook('../data/studies/BRCA/JCI71180sd2.xlsx')
BRCA = readWorksheet(BRCAwb,'ScaledImpData',header = FALSE, check.names = FALSE,startRow = 2)
BRCAmunge = as.matrix( BRCA[5:337,11:dim(BRCA)[2]] )
class(BRCAmunge) = 'numeric'
BRCAmunge = data.frame(BRCAmunge,check.names = F)
rownames(BRCAmunge) = BRCA[5:337,1]
colnames(BRCAmunge) = BRCA[1,11:dim(BRCA)[2]]

# Read in the BRCA data we analyzed
readBRCA = readmunged('BRCA')
BRCAanalyzed = readBRCA$analyzed
BRCAtissue = readBRCA$tissue

# Compare
resBRCA = compareDF( 'BRCA',BRCAmunge,BRCAanalyzed, logfile) 


########################################################################################
#3. BRCATang
BRCATangwb = loadWorkbook('../data/studies/BRCATang/s13058-014-0415-9-s2.xlsx')
BRCATang = readWorksheet(BRCATangwb,'Sheet1',header = FALSE, check.names = FALSE,startRow = 1,
                         startCol = 2, rownames = 1)
BRCATangmunge = as.matrix( BRCATang[20:dim(BRCATang)[1],2:dim(BRCATang)[2]] )
class( BRCATangmunge )= 'numeric'
BRCATang_colnames = BRCATang[1,2:dim(BRCATang)[2]]

# We need to rename the normal samples the same way we did before
BT_normalix = which(BRCATang_colnames %in% c('Normal Breast','Nomal Breast'))
BT_numnormal = length(BT_normalix)
BT_normalnames = sapply(seq(1,BT_numnormal),function(x){paste('N',x,sep='')})
BRCATang_colnames[ BT_normalix ] = BT_normalnames
colnames(BRCATangmunge) = BRCATang_colnames

# Read in BRCATang data we analyzed
readBRCATang = readmunged('BRCATang')
BRCATanganalyzed = readBRCATang$analyzed
BRCATangtissue = readBRCATang$tissue

# Compare 
resBRCATang = compareDF( 'BRCATang',BRCATangmunge, BRCATanganalyzed, logfile)

########################################################################################
#4. KIRC
KIRCwb = loadWorkbook('../data/studies/KIRC/Supp Table 2_MedianNormalized.xlsx')
KIRC = readWorksheet(KIRCwb,sheet = 'Median Normalized',header = TRUE,check.names = FALSE,
                     startRow = 2, startCol = 2, rownames = 1)
KIRCmunge = KIRC[12:588,11:dim(KIRC)[2]]

# Drop mannitol
mannitolix = which(rownames(KIRCmunge) == 'mannitol')
KIRCmunge = KIRCmunge[-mannitolix,]

# Convert to numeric
KIRCmunge = as.matrix(KIRCmunge)
class(KIRCmunge) = 'numeric'

# Read in KIRC data we analyzed
readKIRC = readmunged('KIRC')
KIRCanalyzed = readKIRC$analyzed
KIRCtissue = readKIRC$tissue

resKIRC = compareDF( 'KIRC',KIRCmunge,KIRCanalyzed, logfile)

########################################################################################
#5. PAADHussain1
PH1wb = loadWorkbook('../data/studies/PAADHussain1/NCIA-08-10VW CDT (final).xlsx')
PH1 = readWorksheet(PH1wb,sheet = 'ScaledImpData',header = TRUE,check.names = FALSE,
                     startRow = 2, startCol = 1, rownames = 1)

PH1munge = as.matrix( PH1[9:dim(PH1)[1],11:dim(PH1)[2]] )
class(PH1munge) = 'numeric'

# Read in PH1 data we analyzed
readPH1 = readmunged('PAADHussain1')
PH1analyzed = readPH1$analyzed
PH1tissue = readPH1$tissue

resPH1 = compareDF( 'PAADHussain1',PH1munge,PH1analyzed, logfile)

########################################################################################
#6. PAADHussain2
PH2wb = loadWorkbook('../data/studies/PAADHussain2/NCIA-20-11VW CDT 14Feb2012.xlsx')
PH2 = readWorksheet(PH2wb,sheet = 'ScaledImpData',header = TRUE,check.names = FALSE,
                    startRow = 1, startCol = 2, rownames = 1)

PH2munge = as.matrix( PH2[9:dim(PH2)[1],11:dim(PH2)[2]] )
class(PH2munge) = 'numeric'

# Read in PH2 data we analyzed
readPH2 = readmunged('PAADHussain2')
PH2analyzed = readPH2$analyzed
PH2tissue = readPH2$tissue

resPH2 = compareDF( 'PAADHussain1',PH2munge,PH2analyzed, logfile)

########################################################################################
#7. OV
OVwb = loadWorkbook('../data/studies/OV/LOUI-01-10VW CDT d1-22July2010.xlsx')
OV = readWorksheet(OVwb,sheet = 'ScaledImpData',header = TRUE,check.names = FALSE,
                   startRow = 2,startCol = 1,rownames = 1)

OVmunge = as.matrix( OV[10:dim(OV)[1], 11:dim(OV)[2]] )
class(OVmunge) = 'numeric'

# Read in OV data we analyzed
readOV = readmunged('OV')
OVanalyzed = readOV$analyzed
OVtissue = readOV$tissue

resOV = compareDF( 'OV',OVmunge,OVanalyzed,logfile)

########################################################################################
# 8. PAAD: generated by Rabinowitz group, data is normalized in natural units of ion counts, and then log2-transformed. Simply untransform the data
PAADwb = loadWorkbook('../data/studies/PAAD/KamphorstNofal_Metabolomics_Data.xlsx')
PAAD = readWorksheet(PAADwb,sheet = 'Processed_Data',header = TRUE,check.names = FALSE,
                     startRow = 1,startCol = 1,rownames = 1)
PAADmunge = PAAD
PAADmunge[which(PAADmunge == 'NA',arr.ind = TRUE)] = NA
PAADmunge = as.matrix(PAADmunge)
class(PAADmunge) = 'numeric'

# Convert to natural units
PAADmunge = 2^PAADmunge

# Calculate the minimum, non-NA value for each row
PAAD_rowmins = apply(PAADmunge,1,function(x){min(x,na.rm = TRUE)})

# For each row, replace NA with min value
for (rname in rownames(PAADmunge)){
  naix = which(is.na(PAADmunge[rname,]))
  PAADmunge[rname,naix] = PAAD_rowmins[rname]
}

# Compare to what we analyzed
readPAAD = readmunged('PAAD')
PAADanalyzed = readPAAD$analyzed
PAADtissue = readPAAD$tissue

resPAAD = compareDF( 'PAAD',PAADmunge,PAADanalyzed,logfile)

########################################################################################
#9. PRADLODA
PLwb = loadWorkbook('../data/studies/PRADLODA/131853_1_supp_2658512_nbpsn6.xlsx')
PL = readWorksheet(PLwb,sheet = 'Post-normalization Human tumors',header = TRUE,check.names = FALSE,
                   startRow = 2,startCol = 1,rownames = 1)

# Drop the columns which are all NA
c2drop = which(colnames(PL) %in% c('Col27','Col28','Col89'))
PLmunge = PL[,-c2drop]

# There are 2 docasapentanoates, we modified the rownames manually
rownames(PLmunge)[130] = 'docosapentaenoate (n3 DPA; 22:5n3)'
rownames(PLmunge)[131] = 'docosapentaenoate (n6 DPA; 22:5n6)'

# Drop the last column, where we don't know if its a tumor or not
PLmunge = PLmunge[,-dim(PLmunge)[2]]

# Read in PRADLODA data we analyzed
readPL = readmunged('PRADLODA')
PLanalyzed = readPL$analyzed
PLtissue = readPL$tissue

resPL = compareDF( 'PRADLODA',PLmunge,PLanalyzed,logfile)

########################################################################################
#10. LGG
LGGwb = loadWorkbook('../data/studies/LGG/glioma_metabolomics_2001_normalized.xlsx')
LGGvalues = readWorksheet(LGGwb,sheet = 'ScaledImpData',header = FALSE,check.names = FALSE,
                   startRow = 14,startCol = 2,rownames = 1)

# The line below gets the header, and it seems to skip the empty leading entries, so we don't have to 
# trim the first 10 values off
LGGheader = readWorksheet(LGGwb,sheet = 'ScaledImpData',header = FALSE,check.names = FALSE,
                          startRow = 1,startCol = 2,endRow = 1,rownames = 1)
LGGmunge = LGGvalues[,11:dim(LGGvalues)[2]]
colnames(LGGmunge) = LGGheader

readLGG = readmunged('LGG')
LGGanalyzed = readLGG$analyzed
LGGtissue = readLGG$tissue

resLGG = compareDF( 'LGG',LGGmunge,LGGanalyzed,logfile)

########################################################################################
# Print results
print(resBLCA)
print(resBRCA)
print(resBRCATang)
print(resKIRC)
print(resPH1)
print(resPH2)
print(resOV)
print(resPAAD)
print(resPL)
print(resLGG)
