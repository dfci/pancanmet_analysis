# Make sure all pairs have the same clinical data

rm(list = ls())
setwd('/Users/ereznik/pancanmet_analysis/checks/')

clin = read.csv('/Users/ereznik/pancanmet_analysis/data/merged_metabolomics/clinfeatures.csv',
                header = TRUE)

pairs = read.csv('/Users/ereznik/pancanmet_analysis/data/merged_metabolomics/tumor_normal_pairs.csv',
                 header = TRUE)

pairs$NormalID = paste(pairs$Study,pairs$Normal.Sample,'Normal',sep = ':')
pairs$TumorID = paste(pairs$Study,pairs$Tumor.Sample,'Tumor',sep = ':')
rownames(clin) = paste(clin$study,clin$sample,clin$type,sep = ':')

# for now, drop the one non-unique normal/tumor pair
pairs = pairs[-which(duplicated(pairs$NormalID)),]
rownames(pairs) = pairs$NormalID

normal = clin[which(clin$type == 'Normal'),]
tumor = clin[which(clin$type == 'Tumor'),]

normal$TumorGrade = tumor[pairs[rownames(normal),'TumorID'],'grade']
normal$TumorStage = tumor[pairs[rownames(normal),'TumorID'],'stage']

table(normal$stage,normal$TumorStage)
table(normal$grade,normal$TumorGrade)
