# Script to compare old and new clinical calculations

rm(list = ls())
setwd('/Users/ereznik/pancanmet_analysis/checks/')
library(ggplot2)
library(cowplot)
library(XLConnect)
library(metap)

# Determine whether we want to use stage or grade
f2use = 'grade'

# Read in old and new results
new = readWorksheetFromFile('../results/clinical/updated results stage and grade Jan 2017 max1.xlsx',sheet = f2use,
                              header = TRUE,rownames = 1,startRow = 2)
old = readWorksheetFromFile('../results/clinical/updated results stage and grade v2_OLDRESULTS.xlsx',sheet = f2use,
                            header = TRUE,rownames = 1,startRow = 2)

for (d2use in c('new','old')){
  d = get(d2use)
  cnames = colnames(d)
  rnames = rownames(d)
  d =apply(d,2,as.numeric)
  rownames(d) = rnames
  colnames(d) = cnames
  assign(d2use,d)
  
}

# Intersect
ixr = intersect(rownames(new),rownames(old))
ixc = intersect(colnames(new),colnames(old))

# for each intersecting column, compare p-values
for (col in ixc){
  d2plot = data.frame(new[ixr,col],old[ixr,col])
  colnames(d2plot) = c('new','old')
  d2plot = apply(d2plot,2,log10)
  d2plot = as.data.frame(d2plot)
  
  # Plot
  print(qplot(d2plot$old,d2plot$new) + theme_bw() + ggtitle(col) + geom_abline())
  
}