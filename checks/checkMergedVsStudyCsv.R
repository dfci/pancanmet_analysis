verbose <- TRUE

mergedDat <- read.table(paste0("import/tempdir/merged_metabolomics.csv"), sep=",", header=TRUE, row.names=1, check.names=FALSE, stringsAsFactors = FALSE)

mergedMapping <- read.table(paste0("import/tempdir/merged_mapping.csv"), sep=",", header=TRUE, row.names=1, check.names=FALSE, stringsAsFactors = FALSE)
rownames(mergedMapping) <- rownames(mergedDat)

for(study in colnames(mergedMapping)) {
  #study <- "BLCA"
  
  # READ DATE
  fname <- paste('import/tempdir/',study,'_met.csv',sep='')
  tempmet <- read.csv(fname,sep = ',',header = FALSE,skip = 2,row.names = 1)
  
  #colnames(dat) <- paste0(study, ":", colnames(dat), ":", dat[1,])
  #dat <- dat[2:nrow(dat),]
  
  sampname <- scan(fname, nlines = 1, what = character(), sep=',')
  samptype <- scan(fname, skip = 1, nlines = 1, what = character(), sep=',')
  
  tmpType <- sapply(samptype, function(x) {
    paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
  })
  
  # Make sure below that the third field, indicating tumor/normal/metastasis status, is in sentence case
  metcnames <- sapply(2:length(sampname), function(x){ paste(study,sampname[x],tmpType[x],sep=':') })
  
  dat <- data.frame(tempmet)
  colnames(dat) <- metcnames
  
  tmp <- sapply(dat, as.numeric)
  rownames(tmp) <- rownames(dat)
  colnames(tmp) <- colnames(dat)
  dat <- tmp
  
  for(i in 1:nrow(mergedMapping)) {
    datRow <- mergedMapping[i, study]
    mergedRow <- rownames(mergedMapping)[i]
    
    for(j in 1:ncol(dat)) {
      curCol <- colnames(dat)[j]
      
      if(verbose) {
        cat("I: ", i, " J: ", j, "\n")
      }
      
      if(!is.na(datRow)) {
        #if(dat[datRow, curCol] != mergedDat[mergedRow, curCol]) {
        # Put check that looks for "bigger" errors because write.table writes only 15 decimal points
        if(abs(mergedDat[mergedRow, curCol] - dat[datRow, curCol]) > 1e-10) {  
          stop("Data: ", datRow, ": ", dat[datRow, curCol], " not equal to Merged Data: ", mergedRow, ": ", mergedDat[mergedRow, curCol])
        } else {
          if(verbose) {
            cat("Study: ", study, " Row: ", datRow, " and ", mergedRow, "\n")
          }
        }
      }
    }
  }
  
  #colnames(dat)[i]
  #rownames(dat)[j]
}





