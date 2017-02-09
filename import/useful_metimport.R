# Library of useful functions for merging metabolomics data
library(data.table)

getmode <- function(x) {
  ux <- unique(x)
  tab <- tabulate(match(x, ux)); 
  ux[tab == max(tab)]
}

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), tolower(substring(s, 2)),
        sep="", collapse=" ")
}

getbigdata <- function( pathtodata ){
  
  # A function using data tables to read in a big metabolomics data file in csv format and 
  #   maintains the format of the header (colons don't get deleted). 
  
  temp = fread(pathtodata,sep = ',',header = FALSE,skip = 1)
  header = read.csv(pathtodata, header = TRUE, nrow = 1,check.names = FALSE)
  alldata = data.frame( temp, row.names = 1 )
  colnames(alldata) = colnames(header)[2:dim(header)[2]]
  
  return(alldata)
  
}
