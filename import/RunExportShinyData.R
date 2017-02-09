# Update pancanmet Package Data

# NOTE: Assumption the root directories of pancanmet_analysis and pancanmet are in the same folder
metDataFile <- file.path("data", "merged_metabolomics", "merged_metabolomics.csv")
fcDataFile <- file.path("data", "merged_metabolomics", "tnpaired_fc.csv")
clinFeaFile <- file.path("data", "merged_metabolomics", "clinFeatures.csv")

outputFile <- file.path("..", "pancanmet", "data", "pancanmet_shinydata.Rdata")
shinyOutputFile <- file.path("..", "pancanmet", "inst", "shinyApp", "www", "db", "pancanmet_shinydata.Rdata")
  
# METDATA
metDat <- read.table(metDataFile, sep=",", header=TRUE, row.names=1, stringsAsFactors = FALSE, check.names = FALSE)

tmp <- strsplit(colnames(metDat), ":")
study <- unlist(lapply(tmp, "[[", 1))
type <- unlist(lapply(tmp, "[[", 3))

metdata <- t(metDat)
metdata <- as.data.frame(metdata)
metdata <- cbind(metdata, "Type"=type, "Study"=study, stringsAsFactors = FALSE)

# FC
fcDat <- read.table(fcDataFile, sep=",", header=TRUE, row.names=1, stringsAsFactors = FALSE, check.names = FALSE)

tmp <- strsplit(colnames(fcDat), ":")
study <- unlist(lapply(tmp, "[[", 1))

fc <- t(fcDat)
fc <- as.data.frame(fc)
fc <- cbind(fc, "Study"=study, stringsAsFactors = FALSE)

# CLINICAL FEATURES
clinFeaDat <- read.table(clinFeaFile, sep=",", header=TRUE, stringsAsFactors = FALSE, check.names = FALSE)

save(metdata, fc, clinFeaDat, file=outputFile)
save(metdata, fc, clinFeaDat, file=shinyOutputFile)
