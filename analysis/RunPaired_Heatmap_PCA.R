# Function to visualize paired meabolomics data

rm(list = ls())
setwd('/Users/ereznik/pancanmet_analysis/analysis/')
library(ggplot2)
library(pheatmap)
source('plottingconventions.R')

# Modify names
names2plot = c('BLCA' = 'Bladder','BRCA' = 'Breast Terunuma','BRCATang' = 'Breast Tang',
               'KIRC' = 'Clear Cell Kidney','LGG' = 'Glioma','OV' = 'Ovarian','PAAD' = 'Pancreatic Kamphorst',
               'PAADHussain1' = 'Pancreatic Zhang1','PAADHussain2' = 'Pancreatic Zhang2',
               'PRAD' = 'Prostate Sreekumar','PRADLODA' = 'Prostate Priolo')

# Set color palette
redblue = colorRampPalette(c("navy", "white", "firebrick3"))(50)

# Read in paired data
d = read.csv('../data/merged_metabolomics/tnpaired_fc.csv',
             check.names = FALSE,stringsAsFactors = FALSE,header = TRUE,row.names = 1)

# Take log2
d = log2(d)

# Get only metabolites measured across all studies
d = d[which(complete.cases(d)),]

# Do principal components analysis
pcares =  prcomp(d,center = TRUE,scale = TRUE)
presults = data.frame(pcares$rotation )
presults$Study = sapply(rownames(presults),function(x){strsplit(x,'\\:')[[1]][1]})
presults$Name = names2plot[presults$Study]
presults$Tissue = sourcetissue[presults$Study]
ggplot(presults,aes(PC1,PC2,color = Tissue,shape = Tissue,fill = Tissue)) + geom_point() + stat_ellipse(geom = "polygon",alpha = 0.3,level = 0.5,type = 'norm')

# Make a heatmap
maxval = 4
d[which(d > maxval,arr.ind = TRUE)] = maxval
d[which(d < -maxval, arr.ind = TRUE)] = -maxval
annotation_col = data.frame(colnames(d))
rownames(annotation_col) = annotation_col[,1]
annotation_col$Study = sapply(rownames(annotation_col),function(x){strsplit(x,'\\:')[[1]][1]})
annotation_col$Name = names2plot[annotation_col$Study]
annotation_col$Tissue = sourcetissue[annotation_col$Study]

# Remove irrelevant first column
annotation_col = annotation_col[,-c(1,2,3),drop = FALSE]

# Specify colors
ann_colors = list(
  Time = c("white", "firebrick"),
  CellType = c(CT1 = "#1B9E77", CT2 = "#D95F02"),
  GeneClass = c(Path1 = "#7570B3", Path2 = "#E7298A", Path3 = "#66A61E"))

pheatmap(d,color = redblue,
         clustering_method = "complete",
         show_colnames = FALSE,treeheight_row = 0,treeheight_col = 0,
         annotation_col = annotation_col,clustering_distance_rows = 'correlation',
         clustering_distance_cols = 'correlation')
