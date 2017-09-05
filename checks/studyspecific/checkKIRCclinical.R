# Script to confirm that we have correct stage and grade for KIRC. Can only be run from Ed's machine.
setwd('/Users/ereznik/pancanmet_analysis/checks/studyspecific/')
library(XLConnect)

d1 = readWorksheetFromFile('/Users/ereznik/Documents/ccrc/checks/randomdata/metabolomics clinical without ccPap_June2016_fromAri.xlsx',
                                   sheet = 'Sheet1')
rownames(d1) =d1$Metabolon.ID.Tumor
d2 = read.csv('../../data/studies/KIRC/KIRC_clin.csv',header = TRUE,row.names = 1)
colnames(d2)[3] = 'AJCC.'

ixcols = intersect(colnames(d1),colnames(d2))
for (col in c('PathGrade1','AJCC.')){
  print(col)
  print(table(d1[,col],d2[rownames(d1),col]))
}

print(length(which(d1$ID.Paired.Normal != d2[rownames(d1),'ID.Paired.Normal'])))
print(length(which(d1$ID.Paired.Normal == d2[rownames(d1),'ID.Paired.Normal'])))
