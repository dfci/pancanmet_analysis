library("RColorBrewer");
library("glmnet");
library("beeswarm");

setwd('/Users/ereznik/Documents/pancanmet/analysis/')

# These are functions
checkMetabolites <- function(study, studyData, clinicalData, clinicalVariable, customName) {
  ssamples <- intersect(colnames(studyData), rownames(clinicalData));
  y <- clinicalData[ssamples, clinicalVariable];
  names(y) <- ssamples;
  if(sum(is.na(y)) == length(y)) {
    cat("Study", study, "has no clinical attributes for: ", clinicalVariable, ". Skipping...\n");
    return();
  }
  
  knowns <- which(!is.na(y));
  ssamples <- ssamples[knowns];
  y <- y[ssamples];
  x <- t(studyData[, ssamples]);
  
  results <- rbind(NULL);
  
  # Now do it the other way around
  for(aMet in colnames(x)) {
    my <- x[, aMet];
    cx <- y;
    tmpCoefs <- summary(lm(cx ~ my))$coefficients;
    if(nrow(tmpCoefs) > 1) {
      results <- rbind(results, tmpCoefs["my", ]);
      rownames(results)[nrow(results)] <- aMet;
    }
  }
  
  p.raws <- results[, 4];
  p.fdrs <- p.adjust(p.raws, method="fdr");
  results <- cbind(results, p.fdrs);
  colnames(results)[5] <- "FDR";
  results <- results[sort(p.raws, decreasing=F, index.return=T)$ix, ];
  write.table(file=sprintf("../results/fitplots/lm_results_%s_%s_%s.tsv", study, clinicalVariable, customName), results, sep="\t", quote=F, col.names=NA);
  
  # Now, if we have significant hits; let's plot them
  sigMets <- names(which(sort(p.fdrs) < .05));
  if(length(sigMets) > 0) {
    pdf(file=sprintf("fitplots/lm_significant_plots_%s_%s_%s.pdf", study, clinicalVariable, customName));
    for(sigMet in sigMets) {
      metstudydf <- data.frame(row.names=row.names(x), Metabolite=x[, sigMet], Clinical=y);
      beeswarm(Metabolite ~ Clinical, data=metstudydf, pch=16, corral="wrap", ylab=sprintf("%s", sigMet), xlab=sprintf("%s", clinicalVariable), main=sprintf("Study: %s \n %s vs clinical attribute \n p-val: %f | FDR: %f", study, customName, results[sigMet, 4], p.fdrs[sigMet]));
      boxplot(Metabolite ~ Clinical, data=metstudydf, outline=F, add=T);
    }
    tmp <- dev.off();
  }
}; # end of function definition


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
normalCounts <- table(hugeData[1, sampleTypes == "Normal"]);
hugeCancers <- as.factor(sapply( strsplit(hugeSamples, "\\."), FUN=function(x) { return(x[1]); }));
hugeCancerTypes <- levels(hugeCancers);
# Make the data numeric (this might introduce NAs for empty strings, but it is expected behaviour)
hugeData <- hugeData[3:nrow(hugeData), ];
class(hugeData) <- "numeric";

# this is all the clinical data we have
clinical <- read.csv("../data/merged_metabolomics/clinfeatures.csv", header=T, sep=",");
sampleIds <- sapply( strsplit(allSamples, "\\."), FUN=function(x) { return(x[2]); });

### Let's do the clinical association
for(study in hugeCancerTypes) {
  incSamples <- which(clinical[, "Study"] == study);
  if(length(incSamples) == 0) {
    cat("Study", study, "has no clinical attributes. Skipping...\n");
    next;
  }
  
  sclinical <- clinical[incSamples, ];
  rownames(sclinical) <- sclinical[, "Sample"];
  
  tsclinical <- sclinical[which(sclinical[, "Type"] == "Tumor"), ];
  nsclinical <- sclinical[which(sclinical[, "Type"] == "Normal"), ];
  
  studySamples <- which(cancers == study);
  studyData <- allData[, studySamples];
  colnames(studyData) <- sampleIds[studySamples];
  # Filter out all-NAs
  studyData <- log(studyData[which(rowSums(is.na(studyData)) < ncol(studyData)), ]);
  
  studyTumors <- which((hugeCancers == study) & (sampleTypes == "Tumor"));
  sTumorData <- hugeData[, studyTumors];
  tumorData <- log(sTumorData[which(rowSums(is.na(sTumorData)) < ncol(sTumorData)), ]);
  colnames(tumorData) <- hugeSampleIds[studyTumors];
  
  studyNormals <- which((hugeCancers == study) & (sampleTypes == "Normal"));
  sNormalData <- hugeData[, studyNormals];
  normalData <- log(sNormalData[which(rowSums(is.na(sNormalData)) < ncol(sNormalData)), ]);
  colnames(normalData) <- hugeSampleIds[studyNormals];
  
  for(clinicalVariable in c("Stage", "Grade")) {
    checkMetabolites(study, studyData, tsclinical, clinicalVariable, "FC");
    checkMetabolites(study, tumorData, tsclinical, clinicalVariable, "Tumor");
    checkMetabolites(study, normalData, nsclinical, clinicalVariable, "Normal");
  }
  
}


metabolites <- rownames(allData);
noOfMetabolites <- length(metabolites);
metaboliteCounts <- rep(0, length(metabolites));
names(metaboliteCounts) <- metabolites;
studyMetabolites <- list();

for(metabolite in metabolites) {
  for(study in cancerTypes) {
    studySamples <- which(cancers == study);
    studyData <- allData[metabolite, studySamples];
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
stackedCounts <- matrix(0, nrow=length(metaboliteBins), ncol=length(cancerTypes));
colnames(stackedCounts) <- cancerTypes;
rownames(stackedCounts) <- names(metaboliteBins);

for(study in cancerTypes) {
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
allSampleCounts <- matrix(0, nrow=2, ncol=length(cancerTypes));
colnames(allSampleCounts) <- cancerTypes;
rownames(allSampleCounts) <- c("Tumor", "Normal");
sampleCounts <- table(cancers);
allSampleCounts["Tumor", names(sampleCounts)] <- sampleCounts;
allSampleCounts["Normal", names(sampleCounts)] <- normalCounts[names(sampleCounts)];

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