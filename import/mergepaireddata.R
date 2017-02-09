# Function to merge metabolomics data for paired tumor/normal samples

source('useful_metimport.R')

mergepaireddata = function(pairspath,logfile){
  
  # Read in big data file
  alldata = getbigdata('tempdir/merged_metabolomics.csv')
  
  # Read in the pairs data
  pairs = read.csv(pairspath,header = TRUE)
  
  # Modify the pairs file to make rownames identical to alldata file
  normalIDs = sapply(1:dim(pairs)[1],function(x){paste(pairs[x,3],pairs[x,1],'Normal',sep=':')})
  tumorIDs = sapply(1:dim(pairs)[1],function(x){paste(pairs[x,3],pairs[x,2],'Tumor',sep=':')})
  
  # Check that these column names are in alldata
  normcheck = length(which(normalIDs %in% colnames(alldata)))
  tumcheck = length(which(tumorIDs %in% colnames(alldata)))
  
  if (normcheck!=dim(pairs)[1] | tumcheck!=dim(pairs)[1]){
    stop('Some of the samples in the paired annotation file are not in the merged metabolomics file.')
  }
  
  # Make metabolomics files
  normal = alldata[,normalIDs]
  tumor = alldata[,tumorIDs]
  
  paired = tumor/normal
  write.csv(paired,'tempdir/tnpaired_fc.csv')
}