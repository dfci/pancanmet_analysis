# A script to confirm that the only changes between the deprectated KIRC metabolomics and the new one is one row.

rm(list = ls())
setwd('/Users/ereznik/pancanmet_analysis/checks/studyspecific/')

df1 = read.csv( '../../data/metabolomics/KIRC/KIRC_metdata.csv', header = TRUE,row.names = 1 )
df2 = read.csv( '../../data/metabolomics/KIRC/KIRC_metdata_DEPRECATEDinAPRIL2016.csv', header = TRUE,row.names = 1 )

# Convert to numeric
df1 = as.matrix( df1[-1,] )# drop tumor/normal row
df2 = as.matrix( df2[-1,] )

class(df1) = 'numeric'
class(df2) = 'numeric'

# Remove the row corresponding to methylglutaroylcarnitine
methix = which(rownames(df1) == 'methylglutaroylcarnitine')
df1 = df1[-methix,]

# Function which checks that rownames match, colnames match, and values match for two dataframes
diffrnames = which( sort(rownames(df1)) != sort(rownames(df2)))
diffcnames = which( sort(colnames(df1)) != sort(colnames(df2)))
if (length(diffrnames)!=0){print(diffrnames);warning(paste(study,'This rowname in munged data does not match.'))}
if (length(diffcnames)!=0){print(diffcnames);warning(paste(study,'This colname in munged data does not match'))}
diffvals = df1 - df2[rownames(df1),colnames(df1)]
diffvals[which(is.na(diffvals),arr.ind= TRUE)] = 0
maxdiff = max(diffvals)
if (max( diffvals ) > 1e-6){stop(paste(study,'Values in munged and re-imported data do not match.'))}