# Script to manually prepare some metabolomics data which was not available as appropriately scaled/normalized

# Here are the studies this script addresses:
# 1. PAAD: the data is reported as log2-normalized after normalization. Missing values are reported as NA. 
# We impute the data to minimum measured value, and then re-exponentiate.

rm(list = ls())
analysisDir = normalizePath(".")
setwd(analysisDir)

library(data.table)
library(XLConnect)




##############################################################################################################
# Import PAAD
PAADwb = loadWorkbook('../data/studies/PAAD/KamphorstNofal_Metabolomics_Data.xlsx')
PAAD = readWorksheet(PAADwb,sheet = 'Processed_Data',header = TRUE,check.names = FALSE,
                     startRow = 1,startCol = 1,rownames = 1)

PAADfinal = matrix(NA,dim(PAAD)[1],dim(PAAD)[2])
rownames(PAADfinal) = rownames(PAAD)
colnames(PAADfinal) = colnames(PAAD)

# For every row in PAAD, find the minimum value, set all NAs to it, and re-exponentiate
for (rname in rownames(PAAD)){
  rdata = PAAD[rname,]
  naix = which(rdata == 'NA')
  rdata[naix] = NA

  if (length(which(is.na(rdata))) == length(rdata)){print(rname)}
  
  class(rdata) = 'numeric'
  
  impval = min(rdata,na.rm = TRUE)
  
  # Assign imputed data and re-exponentiate
  rdata[naix] = impval
  rdata = 2^rdata
  
  PAADfinal[rname,] = rdata
  
}

# Add row for tissue type
PAADsamples = colnames(PAADfinal)
PAADtissue = rep('Normal',dim(PAADfinal)[2])
tumorix = grep('Tumor',PAADsamples)
PAADtissue[tumorix] = 'Tumor'
names(PAADtissue) = colnames(PAADfinal)

# We are having issues writing the column names with extra quotes being added,
# so to do so we just added an extra row with column names, and omit the actual column names when writing the csv file.

# Bind
PAADfinal = rbind(PAADsamples,PAADtissue,PAADfinal)

# Add metabolite names
PAADfinal = cbind(rownames(PAADfinal),PAADfinal)

# Correct column names
colnames(PAADfinal) = c('Metabolite',colnames(PAAD))
write.table(PAADfinal,'../data/metabolomics/PAAD/PAAD_metdata.csv',row.names = FALSE,col.names = FALSE,quote = c(1),sep = ',')


##############################################################################################################
# 2. Import PRAD
PRADwb = loadWorkbook('../data/studies/PRAD/Sreekumar_etal_2009_tissue_metabolites (1).xls')
PRAD = readWorksheet(PRADwb, sheet = 'MetaboliteData',header = TRUE,check.names = FALSE,
                     startRow = 1,startCol = 1,rownames = 1)

# Drop first 2 columns
PRAD = PRAD[,3:dim(PRAD)[2]]

# Only used the named metabolites
PRAD = PRAD[1:626,]

# Drop all rows that contain 'X-', these are unknowns
rows2drop = grep('X\\-',rownames(PRAD))
print('Dropping these rows from PRAD:')
print(rownames(PRAD)[rows2drop])
if (length(which(rows2drop != 194:623))!=0){stop('We did not drop PRAD unknown metabolites correctly.')}
PRAD = PRAD[-rows2drop,]

# Check that the last metabolite is Xylitol
if (rownames(PRAD)[dim(PRAD)[1]] !='Xylitol'){stop('PRAD data not imported correctly.')}

# Remove control columns from analysis (should be 1:42)
notcontrolcols = grep('C',colnames(PRAD),invert = TRUE)

# For each row, find the median measured value and the imputed value, then normalize
PRADfinal = data.frame( matrix(NA,0,length(notcontrolcols)) )
colnames(PRADfinal) = colnames(PRAD)[notcontrolcols]

# There are two entries for 2-hydroxybutyrate, GC and LC, and we will only keep the second one, which is LC
hydroxyflag = TRUE
for (rname in rownames(PRAD)){
 
  if (rname == "2-Hydroxybutyrate (AHB)" & hydroxyflag){
    # Skip this entry, it's the GC measurement of 2-hydroxybutyrate
    hydroxyflag = FALSE
    next
  }
  
  rdata = PRAD[rname,notcontrolcols]
  class(rdata) = 'numeric'
  
  naix = which(is.na(rdata))
  medval = median(rdata,na.rm = TRUE)
  impval = min(rdata,na.rm = TRUE)
  
  # Assign imputed data and re-exponentiate
  rdata[naix] = impval
  
  # Rescale by median
  rdata = rdata/medval
  
  if (rname == "2-Hydroxybutyrate (AHB).1" & !(hydroxyflag)){
    # Save this entry with the correct metabolite name
    rname2save = "2-Hydroxybutyrate (AHB)"
  }else{
    rname2save = rname
  }
  
  
  # Save
  PRADfinal[rname2save,] = rdata
  
}

# Make a tissue row
PRADtissue = rep('Normal',dim(PRADfinal)[2])
tumorix = grep('T',colnames(PRADfinal))
metix = grep('M',colnames(PRADfinal))
PRADtissue[tumorix] = 'Tumor'
PRADtissue[metix] = 'Metastasis'

# We are having issues writing the column names with extra quotes being added,
# so to do so we just added an extra row with column names, and omit the actual column names when writing the csv file.

# Bind
PRAD2write = rbind(colnames(PRADfinal),PRADtissue,PRADfinal)

# Set the first 2 rownames
rownames(PRAD2write)[1:2] = c('SampleName','SampleType')

# Correct column names
write.table(PRAD2write,'../data/metabolomics/PRAD/PRAD_metdata.csv',row.names = TRUE,col.names = FALSE,quote = c(1),sep = ',')

# Compare to deprecated version, imported using Python
PRAD_deprecated = read.csv('../data/metabolomics/PRAD/PRAD_metdata_DEPRECATED_April2016.csv',header = TRUE,row.names = 1)
PRAD_deprecated_values = as.matrix( PRAD_deprecated[-1,] )
class(PRAD_deprecated_values) = 'numeric'

# Check column names, rownames, and tissue types
PRAD2writemets = rownames(PRAD2write)[3:dim(PRAD2write)[1]]
PRAD_deprecated_mets = rownames(PRAD_deprecated)[2:dim(PRAD_deprecated)[1]]
if ( length( which(  PRAD2writemets != PRAD_deprecated_mets ) ) != 0){
  warning('PRAD rownames not equal to deprecated file.')
  print(setdiff(PRAD2writemets,PRAD_deprecated_mets))
}
if ( !(all.equal( colnames(PRAD2write),colnames(PRAD_deprecated) ) ) ){
  stop('PRAD colnames not equal to deprecated file.')
}

PRAD_deprecated_tissues = sapply(PRAD_deprecated[1,],as.character)
if ( length( which( as.character(PRAD2write[2,]) != PRAD_deprecated_tissues ) ) != 0 ){
  stop('PRAD rownames not equal to deprecated file.')
}

PRADdiffvals = PRADfinal - PRAD_deprecated_values
if (max(abs(PRADdiffvals)) > 1e-6){stop('The values in the newly imported PRAD file do not match.')}
