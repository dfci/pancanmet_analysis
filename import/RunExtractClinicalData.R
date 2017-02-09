# Extract clinical data from xlsx files 
## Attempt to grab STAGE and GRADE where possible 

library(xlsx)

results <- list()

# PROCESS STUDIES ----

## BLCA (PDF MUST BE MANUAL, IGNORED, MOST GRADES HIGH (ALL BUT 2), HIGHLY UNLIKELY TO PRODUCE ANYTHING INTERESTING)

## BRCA (DONE)
cancerType <- "BRCA"
file <- file.path("data", "studies", cancerType, "BRCA_alldata.xlsx")
sheetName <- "OrigData"
rowIndex <- c(2:18)
colIndex <- c(10:142)
header <- TRUE

dat <- read.xlsx(file, sheetName=sheetName, rowIndex=rowIndex, colIndex=colIndex, 
                 colClasses=NA, header=header, as.data.frame=TRUE, stringsAsFactors=FALSE, check.names = FALSE)

dat2 <- dat[,2:ncol(dat)]
rownames(dat2) <- dat[,1]

tmpType <- sapply(substr(dat2["TISSUE TYPE", ], 5, 10), function(x) {
  paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
})

tmp <- data.frame(sample=colnames(dat2),
                  study=rep(cancerType, ncol(dat2)),
                  type=tmpType,
                  stage=as.numeric(dat2["STAGE_NUM", ]), 
                  grade=as.numeric(dat2["TUMOR_GRADE", ]), 
                  notes=rep(NA, ncol(dat2)),
                  stringsAsFactors = FALSE)

rownames(tmp) <- colnames(dat2)
brcaResults <- tmp

## BRCATang (DONE)
cancerType <- "BRCATang"
file <- file.path("data", "studies", cancerType, "s13058-014-0415-9-s2.xlsx")
sheetName <- "Sheet1"
rowIndex <- c(1:18)
colIndex <- c(3:33)
header <- TRUE

dat <- read.xlsx(file, sheetName=sheetName, rowIndex=rowIndex, colIndex=colIndex, 
                 colClasses=NA, header=header, as.data.frame=TRUE, stringsAsFactors=FALSE, check.names=FALSE)

dat2 <- dat[,2:ncol(dat)]
rownames(dat2) <- dat[,1]

# Grab stage from TCGA data
clinFile <- file.path("data", "studies", cancerType, "gdac.broadinstitute.org_BRCA.Clinical_Pick_Tier1.Level_4.2016012800.0.0", "BRCA.clin.merged.picked.txt")

clinDat <- read.table(clinFile, sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)

# Get stages 
stageRow <- which(clinDat[,1] == "pathologic_stage")

stageResults <- NULL 

for(i in 1:ncol(dat2)) {
  # NOTE: Some names have trailing whitespace
  curSample1 <- trimws(colnames(dat2)[i])
  curSample2 <- paste0("tcga-", tolower(curSample1))
  
  stageIdx <- which(colnames(clinDat) == curSample2)
  
  if(length(stageIdx) == 1) {
    stageResults <- c(stageResults, clinDat[stageRow, stageIdx])
    cat("Sample: ", curSample2, " Stage: ", clinDat[stageRow, stageIdx], "\n")
  } else {
    stageResults <- c(stageResults, NA)
  }
}

## Fix stages
### Replace a/b in stage 
tmpStageResults <- stageResults
tmpStageResults <- gsub("ia", "i", tmpStageResults)
tmpStageResults <- gsub("ib", "i", tmpStageResults)

### Replace roman numerals and remove "stage
tmpStageResults <- gsub("stage ", "", tmpStageResults, ignore.case = TRUE)
tmpStageResults[tolower(tmpStageResults) == "i"] <- 1
tmpStageResults[tolower(tmpStageResults) == "ii"] <- 2
tmpStageResults[tolower(tmpStageResults) == "iii"] <- 3
tmpStageResults[tolower(tmpStageResults) == "iv"] <- 4
tmpStageResults[tolower(tmpStageResults) == "x"] <- NA

tmp <- data.frame(sample=colnames(dat2),
                  study=rep(cancerType, ncol(dat2)),
                  type=rep("Tumor", ncol(dat2)),
                  stage=tmpStageResults, 
                  grade=rep(NA, ncol(dat2)), 
                  notes=stageResults,
                  stringsAsFactors = FALSE)

# Remove normal samples 
tmp <- tmp[grepl("^B", tmp$sample),]
brcaTangResults <- tmp 

## COAD (IGNORED)

## KIRC (DONE)
cancerType <- "KIRC"
file <- file.path("data", "studies", cancerType, "KIRC_clin.csv")
dat <- read.table(file, sep=",", header=TRUE, stringsAsFactors = FALSE)

results <- NULL

for(i in 1:nrow(dat)) {
  tmp1 <- c(dat[i, "Metabolon.ID.Tumor"], cancerType, "Tumor", dat[i, "AJCC"], dat[i, "PathGrade1"])
  tmp2 <- c(dat[i, "ID.Paired.Normal"], cancerType, "Normal", dat[i, "AJCC"], dat[i, "PathGrade1"])
  results <- rbind(results, tmp1)
  results <- rbind(results, tmp2)
}

colnames(results) <- c("sample", "study", "type", "stage", "grade")
results <- cbind(results, notes=rep(NA, nrow(results)))

kircResults <- results

## LGG (DONE)
cancerType <- "LGG"
file <- file.path("data", "studies", cancerType, "glioma_metabolomics_2001_normalized.xlsx")
sheetName <- "OrigScale"
rowIndex <- c(1:11)
colIndex <- c(12:81)
header <-  FALSE

dat <- read.xlsx(file, sheetName=sheetName, rowIndex=rowIndex, colIndex=colIndex, 
                 colClasses=NA, header=header, as.data.frame=TRUE, stringsAsFactors=FALSE, check.names = FALSE)

dat2 <- dat[,2:ncol(dat)]
rownames(dat2) <- dat[,1]

# Fix grades
grades <- tolower(as.character(dat2["GRADE", ]))

tmpGrades <- grades
tmpGrades[tolower(tmpGrades) == "i"] <- 1
tmpGrades[tolower(tmpGrades) == "ii"] <- 2
tmpGrades[tolower(tmpGrades) == "iii"] <- 3
tmpGrades[tolower(tmpGrades) == "iv"] <- 4
tmpGrades[tolower(tmpGrades) == "x"] <- NA

tmp <- data.frame(sample=as.character(dat2["SAMPLE NAME", ]), 
                  study=rep(cancerType, ncol(dat2)),
                  type=rep("Tumor", ncol(dat2)), 
                  stage=rep(NA, ncol(dat2)),
                  grade=tmpGrades,
                  notes=grades,
                  stringsAsFactors = FALSE)

lggResults <- tmp

## OV (DONE)
cancerType <- "OV"
file <- file.path("data", "studies", cancerType, "LOUI-01-10VW CDT d1-22July2010.xlsx")
sheetName <- "OrigData"
rowIndex <- c(1:10)
colIndex <- c(11:41)
header <- FALSE

dat <- read.xlsx(file, sheetName=sheetName, rowIndex=rowIndex, colIndex=colIndex, 
                 colClasses=NA, header=header, as.data.frame=TRUE, stringsAsFactors=FALSE)

dat2 <- dat[,2:ncol(dat)]
rownames(dat2) <- dat[,1]

### Fix types
tmpType <- as.character(dat2["TISSUE", ])
tmpType[tmpType == "Ovary"] <- "Tumor"
tmpType[tmpType == "Met. Ovary"] <- "Metastasis"

### Fix stages
stages <- as.character(dat2["STAGE", ])
tmpStageResults <- stages
tmpStageResults <- gsub("IA", "I", tmpStageResults)
tmpStageResults <- gsub("IC", "I", tmpStageResults)

tmpStageResults <- gsub("stage ", "", tmpStageResults, ignore.case = TRUE)
tmpStageResults[tolower(tmpStageResults) == "i"] <- 1
tmpStageResults[tolower(tmpStageResults) == "ii"] <- 2
tmpStageResults[tolower(tmpStageResults) == "iii"] <- 3
tmpStageResults[tolower(tmpStageResults) == "iv"] <- 4
tmpStageResults[tolower(tmpStageResults) == "x"] <- NA

### Get stages
results <- data.frame(sample=as.character(dat2["SAMPLE_NAME", ]), 
                      study=rep(cancerType, ncol(dat2)),
                      type=tmpType,
                      stage=tmpStageResults,
                      grade=rep(NA, ncol(dat2)), 
                      notes=stages, 
                      stringsAsFactors = FALSE)

### Remove normals because they are not paired and have no stage
results <- results[results$stage != "Normal",]

ovResults <- results

## PAAD (DONE, NO STAGE/GRADE, IGNORE)
# file <- file.path("data", "studies", "PAAD", "KamphorstNofal_Metabolomics_Data.xlsx")
# sheetName <- "Sample_List"
# rowIndex <- c(1:184)
# colIndex <- c(1:7)
# header <- TRUE
# 
# dat <- read.xlsx(file, sheetName=sheetName, rowIndex=rowIndex, colIndex=colIndex, 
#                  colClasses=NA, header=header, as.data.frame=TRUE, stringsAsFactors=FALSE)
# 
# dat2 <- dat[,2:ncol(dat)]
# rownames(dat2) <- dat[,1]

## PAADHussain1 (DONE, NO STAGE/GRADE, IGNORE)
# file <- file.path("data", "studies", "PAADHussain1", "NCIA-08-10VW CDT (final).xlsx")
# sheetName <- "OrigScale"
# rowIndex <- c(1:9)
# colIndex <- c(11:70)
# header <- TRUE
# 
# dat <- read.xlsx(file, sheetName=sheetName, rowIndex=rowIndex, colIndex=colIndex, 
#                  colClasses=NA, header=header, as.data.frame=TRUE, stringsAsFactors=FALSE)
# 
# dat2 <- dat[,2:ncol(dat)]
# rownames(dat2) <- dat[,1]

## PAADHussain2 (DONE, NO STAGE/GRADE, IGNORE)
# file <- file.path("data", "studies", "PAADHussain2", "NCIA-20-11VW CDT 14Feb2012.xlsx")
# sheetName <- "OrigScale"
# rowIndex <- c(1:8)
# colIndex <- c(12:62)
# header <- TRUE
# 
# dat <- read.xlsx(file, sheetName=sheetName, rowIndex=rowIndex, colIndex=colIndex, 
#                  colClasses=NA, header=header, as.data.frame=TRUE, stringsAsFactors=FALSE)
# 
# dat2 <- dat[,2:ncol(dat)]
# rownames(dat2) <- dat[,1]

## PRAD ()
cancerType <- "PRAD"
file <- file.path("data", "studies", cancerType, "Sreekumar_etal_2009_tissue_metabolites (1).xls")
sheetName <- "SampleInfo"
rowIndex <- c(1:53)
colIndex <- c(1:13)
header <- TRUE

dat <- read.xlsx(file, sheetName=sheetName, rowIndex=rowIndex, colIndex=colIndex, 
                 colClasses=NA, header=header, as.data.frame=TRUE, stringsAsFactors=FALSE)

dat2 <- dat[,2:ncol(dat)]
rownames(dat2) <- dat[,1]

### Keep only tumors and normals
tumorNormalIdx <- which(grepl("^[N|T]", rownames(dat2)))
dat2 <- dat2[tumorNormalIdx, ]

### Get type
tmpType <- character(nrow(dat2))
tmpType[which(grepl("^N", rownames(dat2)))] <- "Normal"
tmpType[which(grepl("^T", rownames(dat2)))] <- "Tumor"

gsMajor <- dat2[, "GSMajor"]
gsMinor <- dat2[, "GSMinor"]

gs <- NULL

for(i in 1:length(gsMajor)) {
  major <- gsMajor[i]
  minor <- gsMinor[i]
  
  if(!is.na(major) && !is.na(minor)) {
    gs <- c(gs, paste0(major, "+", minor))
  } else {
    gs <- c(gs, NA)
  }
}  

tmpGrades <- sapply(gs, function(x) { 
  tmp <- eval(parse(text=x)) 
  
  if(!is.na(tmp)) {
    if(tmp <= 6) {
      return(1)
    } 
    
    if(tmp == 7) {
      return(2)
    }
    
    if(tmp > 7) {
      return(3)
    }    
  } else {
    return(NA)
  }
})

results <- data.frame(sample=rownames(dat2), 
                      study=rep(cancerType, nrow(dat2)),
                      type=tmpType,
                      stage=rep(NA, nrow(dat2)),
                      grade=tmpGrades, 
                      notes=gs,
                      stringsAsFactors = FALSE)

pradResults <- results

## PRADLODA (DONE, ???)
cancerType <- "PRADLODA"

pradLodaPairs <- read.csv(file.path("data", "merged_metabolomics", "tumor_normal_pairs.csv"), stringsAsFactors = FALSE)
pradLodaPairs <- pradLodaPairs[pradLodaPairs$Study == cancerType, ]

file <- file.path("data", "studies", cancerType, "131853_1_supp_2658512_nbpsn6.xlsx")
sheetName <- "Post-normalization Human tumors"
rowIndex <- c(2)
colIndex <- c(2:26)
header <- FALSE

### Get Normal IDs
dat <- read.xlsx(file, sheetName=sheetName, rowIndex=rowIndex, colIndex=colIndex, 
                 colClasses=NA, header=header, as.data.frame=TRUE, stringsAsFactors=FALSE)
normalIds <- as.character(dat)

### Get Tumor IDs
rowIndex <- c(2)
colIndex <- c(28:88)
header <- FALSE

dat <- read.xlsx(file, sheetName=sheetName, rowIndex=rowIndex, colIndex=colIndex, 
                 colClasses=NA, header=header, as.data.frame=TRUE, stringsAsFactors=FALSE)
tumorIds <- as.character(dat)

# Get Gleason scores 
file <- file.path("data", "studies", cancerType, "Sample_information_PrioloC.xlsx")
sheetName <- "Sheet1"
rowIndex <- c(1:63)
colIndex <- c(1:16)
header <- TRUE

dat <- read.xlsx(file, sheetName=sheetName, rowIndex=rowIndex, colIndex=colIndex, 
                   colClasses=NA, header=header, as.data.frame=TRUE, stringsAsFactors=FALSE)

dat2 <- dat[2:nrow(dat),]

# Get Gleason scores
tmpGs <- dat2[,"specimen.Gleason"]
tmpGs[tmpGs == "N/A"] <- NA

# Get IDs
tmpId <- dat2[,"Unique.ID"]

# Keep on complete cases
tmpDat <- cbind(tmpId, tmpGs)
tmpDat <- tmpDat[complete.cases(tmpDat), ]

tmpGrades <- sapply(tmpDat[,2], function(x) { 
  tmp <- eval(parse(text=x)) 
  
  if(!is.na(tmp)) {
    if(tmp <= 6) {
      return(1)
    } 
    
    if(tmp == 7) {
      return(2)
    }
    
    if(tmp > 7) {
      return(3)
    }    
  } else {
    return(NA)
  }
})
  
results <- data.frame(sample=tmpDat[,"tmpId"], 
                      study=rep(cancerType, nrow(tmpDat)),
                      type=rep("Tumor", nrow(tmpDat)),
                      stage=rep(NA, nrow(tmpDat)),
                      grade=tmpGrades, 
                      notes=tmpDat[,2],
                      stringsAsFactors = FALSE)

### Add any missing normal/tumor IDs 

#### Tumors
tmpIds <- tumorIds[!(tumorIds %in% results$sample)]
tmpResultsTumors <- data.frame(sample=tmpIds, 
                      study=rep(cancerType, length(tmpIds)),
                      type=rep("Tumor", length(tmpIds)),
                      stage=rep(NA, length(tmpIds)),
                      grade=rep(NA, length(tmpIds)), 
                      notes=rep(NA, length(tmpIds)),
                      stringsAsFactors = FALSE)

#### Normals 
tmpIds <- normalIds[!(normalIds %in% results$sample)]
tmpResultsNormals <- data.frame() 
naResults <- data.frame(sample=NA, 
                         study=cancerType,
                         type="Normal",
                         stage=NA,
                         grade=NA, 
                         notes=NA,
                         stringsAsFactors = FALSE)

for(i in 1:length(tmpIds)) {
  #i <- 4
  normalId <- tmpIds[i]
  idx <- which(pradLodaPairs$Normal.Sample == normalId)
  
  if(length(idx) > 0) {
    tumorId <- pradLodaPairs$Tumor.Sample[idx]
    
    if(length(tumorId) > 0) {
      curResults <- results[results$sample == tumorId,]
      
      if(nrow(curResults) > 0) {
        curResults$notes <- paste0("Tumor Pair ID: ", tumorId)
      } else {
        curResults <- naResults
      }
    } else {
      curResults <- naResults
    }
  } else {
    curResults <- naResults
  }
  
  curResults$sample <- normalId
  curResults$type <- "Normal"
  tmpResultsNormals <- rbind(tmpResultsNormals, curResults)
}

results <- rbind(results, tmpResultsTumors, tmpResultsNormals)
pradLodaResults <- results

## STAD (IGNORED)

# AGGREGATE DATA ----

allResults <- rbind(brcaResults, brcaTangResults, kircResults, lggResults, ovResults, pradResults, pradLodaResults)

write.table(allResults, file=file.path("data", "merged_metabolomics", "clinFeatures.csv"), sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)
