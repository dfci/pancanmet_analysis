# Script to import data for shiny app

rm(list=ls())
setwd('/Users/ereznik/Documents/pancanmet/shinyapp/')

selectfirst <- function(s) strsplit(s, ":")[[1]][1]
selectlast <- function(s) strsplit(s, ":")[[1]][ length( strsplit(s, ":")[[1]] )]

# Import the fold change data
fc = read.csv('/Users/ereznik/Documents/pancanmet/data/merged_metabolomics/tnpaired_fc.csv', row.names = 1,check.names = FALSE)
splittypes = strsplit(colnames(fc),'\\:')
tissuetype = sapply(splittypes, "[", 1)

fc = data.frame( log2( t(fc) ) )
fcstudy = sapply(rownames(fc),selectfirst)
fc$Study = fcstudy

# Import all of the data, regardless of matched existence
alldata = read.csv('/Users/ereznik/Documents/pancanmet/data/merged_metabolomics/alldata.csv', row.names = 1,,check.names = FALSE)
metdata = data.frame( sapply(alldata[3:dim(alldata)[1],1:(dim(alldata)[2]-2)], function(x) { as.numeric(as.character(x)) } ), check.names = F )
rownames(metdata) = rownames(alldata)[3:dim(alldata)[1]]

metdata = data.frame( log2( t(metdata) ) )
ttype =  sapply(rownames(metdata),selectlast)
study = sapply(rownames(metdata),selectfirst)
metdata$Type = ttype
metdata$Study = study

save(metdata,fc,file = 'pancanmet_shinydata.RData')
