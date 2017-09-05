rm(list = ls())

library("RColorBrewer");
library("glmnet");
library("beeswarm");

setwd('/Users/ereznik/Documents/pancanmet/analysis/')

# This is normalized data (t vs n)
allData <- as.matrix(read.csv("../data/merged_metabolomics/tnpaired_fc.csv", header=T, row.names=1));
allSamples <- colnames(allData);
cancers <- as.factor(sapply( strsplit(allSamples, "\\."), FUN=function(x) { return(x[1]); }));
cancerTypes <- levels(cancers);

# and this is all unnormalized data (t and n)
hugeData <- as.matrix(read.csv("../data/merged_metabolomics/alldata.csv", header=T, row.names=1));
hugeSamples <- colnames(hugeData);
sampleTypes <- as.factor(sapply( strsplit(hugeSamples, "\\."), FUN=function(x) { return(x[length(x)]); }));
hugeSampleIds <- sapply( strsplit(hugeSamples, "\\."), FUN=function(x) { return(x[2]); });
hugeCancers <- as.factor(sapply( strsplit(hugeSamples, "\\."), FUN=function(x) { return(x[1]); }));
tumorCounts <- table(hugeData[1, sampleTypes == "Tumor"]);
normalCounts <- table(hugeData[1, sampleTypes == "Normal"]);
hugeCancers <- as.factor(sapply( strsplit(hugeSamples, "\\."), FUN=function(x) { return(x[1]); }));
hugeCancerTypes <- levels(hugeCancers);
# Make the data numeric (this might introduce NAs for empty strings, but it is expected behaviour)
hugeData <- hugeData[3:nrow(hugeData), ];
class(hugeData) <- "numeric";

# this is all the clinical data we have
clinical <- read.csv("../data/merged_metabolomics/clinfeatures.csv", header=T, sep=",");
sampleIds <- sapply( strsplit(allSamples, "\\."), FUN=function(x) { return(x[2]); });

metabolites <- rownames(hugeData);
noOfMetabolites <- length(metabolites);
metaboliteCounts <- rep(0, length(metabolites));
names(metaboliteCounts) <- metabolites;
studyMetabolites <- list();

for(metabolite in metabolites) {
  if (nchar(metabolite)==0){
    next
  }
  for(study in hugeCancerTypes) {
    studySamples <- which(hugeCancers == study);
    studyData <- hugeData[metabolite, studySamples];
    allNa <- sum(is.na(studyData)) == length(studySamples);
    
    if(!allNa) {
      metaboliteCounts[metabolite] <- metaboliteCounts[metabolite] + 1;
      studyMetabolites[[study]] <- c(metabolite, studyMetabolites[[study]]);
    }
  }
}


pdf("../results/summary.pdf", width=11.7, height=8.27);

# Pie chart
oldPar <- par(mfrow=c(1,3));
metaboliteBins <- table(metaboliteCounts);
#metaboliteBinColors <- brewer.pal(length(metaboliteBins), "Dark2");
metaboliteBinColors <- colorRampPalette(brewer.pal(8, "Spectral"))(length(metaboliteBins)+1)
metabolitePercentage <- round(metaboliteBins / sum(metaboliteBins) * 100);
pieLabels <- paste(metaboliteBins, " metabolites\nin ", names(metaboliteBins), " studies (", metabolitePercentage, "%)", sep="");
pie(metaboliteBins, labels=pieLabels, col=metaboliteBinColors, main="Measured Metabolites");

# Study Distributions
stackedCounts <- matrix(0, nrow=length(metaboliteBins), ncol=length(hugeCancerTypes));
colnames(stackedCounts) <- hugeCancerTypes;
rownames(stackedCounts) <- names(metaboliteBins);

for(study in hugeCancerTypes) {
  theseMetabolites <- studyMetabolites[[study]];
  thisDistrib <- table( metaboliteCounts[theseMetabolites] );
  stackedCounts[names(thisDistrib), study] <- thisDistrib;
}

labelCoordinates <- barplot(
  stackedCounts, 
  main="Metabolite distribution by study",
  col=metaboliteBinColors, legend=paste(names(metaboliteBins), "studies"), 
  horiz=TRUE,
  xlim=c(0, noOfMetabolites * 1.2), xlab="Number of metabolites",
  las=2
);
abline(v=noOfMetabolites, col="darkgray", lty=2);
text(x=noOfMetabolites, y=max(labelCoordinates)/2, labels="Total number of\nmeasured metabolites", srt=90, col="darkgray");
perStudyCounts <- colSums(stackedCounts);
text(x=perStudyCounts + 10, y=labelCoordinates, labels=paste(perStudyCounts, "metabolites"), pos=4);


# Study Patient Counts
allSampleCounts <- matrix(0, nrow=2, ncol=length(hugeCancerTypes));
colnames(allSampleCounts) <- hugeCancerTypes;
rownames(allSampleCounts) <- c("Tumor", "Normal");
sampleCounts <- table(hugeCancers);
allSampleCounts["Tumor", names(sampleCounts)] <- tumorCounts[names(sampleCounts)]
allSampleCounts["Normal", names(sampleCounts)] <- normalCounts[names(sampleCounts)];
allSampleCounts[which(is.na(allSampleCounts),arr.ind = T)] <- 0

tnColors <- rev(brewer.pal(3, "Pastel2")[1:2]);

labelCoordinates <- barplot(
  allSampleCounts,
  main="Sample distribution by study",
  col=tnColors,
  horiz=TRUE,
  xlab="Number of samples",
  las=2,
  legend=rownames(allSampleCounts),
  beside=TRUE,
  xlim=c(0, max(allSampleCounts) * 1.2)
);

text(y=labelCoordinates[1, ], x=allSampleCounts["Tumor", ] + 1, pos=4, labels=paste(allSampleCounts["Tumor", ], " tumors"));
text(y=labelCoordinates[2, ], x=allSampleCounts["Normal", ] + 1, pos=4, labels=paste(allSampleCounts["Normal", ], " normals"));
par(oldPar);

tmp <- dev.off();