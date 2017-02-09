setwd('/Users/ereznik/Documents/pancanmet/analysis/convert_metabolic_ids')
source("convertChemIds.R")

# Name of input study
study = 'PAADHussein2'

# Load Data
metabolonData <- read.csv(paste(study,'startingIDs.txt',sep='_'), sep = '\t',as.is=TRUE)
metabolonData[is.na(metabolonData)] <- c("")

#toIdTypes <- c("Chemical Name", "PubChem CID", "KEGG", "Human Metabolome Database", "ChEBI")
#fromIdTypes <- c("Chemical Name", "PubChem CID", "KEGG", "Human Metabolome Database")
toIdTypes <- c("PubChem CID", "KEGG", "Human Metabolome Database", "ChEBI")
fromIdTypes <- c("PubChem CID", "KEGG", "Human Metabolome Database")

outputMat <- matrix(NA, nrow=nrow(metabolonData), ncol=length(toIdTypes))
colnames(outputMat) <- toIdTypes
rownames(outputMat) <- metabolonData[,1]

# Use with name search; offset helps match the fromIdTypes to the toIdTypes
#outputMat[,1:4] <- as.matrix(metabolonData[,c(1,4:6)])
offset <- 2

# Use without name search
outputMat[,1:3] <- as.matrix(metabolonData[,c(2:4)])
offset <- 1

rows <- 1:nrow(metabolonData)
#rows <- 1:10

for(i in rows) {
	cat("I: ", i, "\n")

	for(j in 1:length(fromIdTypes)) {
		
		# Don't try to strsplit the names
		if(j == 1) {
			ids <- metabolonData[i,j]
		} else {
			ids <- strsplit(metabolonData[i,(j+offset)], ",")[[1]]		
		}
		
		outputIds <- NULL
		
		for(id in ids) {
			if(id != "NA" | is.na(id)) {
				fromIdType <- fromIdTypes[j]
				
				# Make sure not to grab the same ID type as the current fromIdType
				tmp <- sapply(toIdTypes, function(x) convertChemIds(id, fromIdType, x, debug=FALSE))
				
				outputIds <- unlist(lapply(tmp, paste, collapse=","))				
			} else {
				outputIds <- rep("NA", length(toIdTypes))
			}			
			
			# Put the new IDs in the correct columns of the outputMat
			for(k in 1:length(toIdTypes)) {
					
				tmp1 <- strsplit(outputMat[i,k], ",")[[1]]
				tmp2 <- strsplit(outputIds[k], ",")[[1]]
				
				tmp <- c(tmp1, tmp2)
			
				outputMat[i,k] <- paste(tmp, collapse=",")
			}
		}
	}
}

# Make IDs unique and fill in NA for any missing IDs
for(i in rows) {
	for(j in 1:length(toIdTypes)) {
		tmp <- strsplit(outputMat[i,j], ",")[[1]]
		tmp <- tmp[which(tmp != "NA")]
		
		if(length(tmp) > 0) {
			tmp <- unique(tmp) 
		} else {
			tmp <- "NA"
		}
	
		outputMat[i,j] <- paste(tmp, collapse="|")	
	}
}

write.table(outputMat, file=paste(study,'convertedIDs.txt',sep='_'), sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE) 

# Counts of missing values for the original and output data
orgDataCnts <- sapply(2:4, function(x) length(which(metabolonData[,x] == "")))
names(orgDataCnts) <- names(metabolonData[,2:4])

outputCnts <- sapply(1:4, function(x) length(which(outputMat[,x] == "NA")))
names(outputCnts) <- colnames(outputMat[,1:4])

#tmp <- read.table("tmp_name_check5.txt", sep="\t", quote="", header=TRUE, comment="")
#tmpCnts <- sapply(1:4, function(x) length(which(tmp[,x] == "NA")))
#names(tmpCnts) <- colnames(tmp[,1:4])


