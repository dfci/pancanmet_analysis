#' Function to calculate co-variation of a dataframe
#'
#' @param df: dataframe for which to calculate a covariation matrix. Calculated across all pairs of rows using spearman correlations
#' @param pcrit: optional, if provided, set any correlations with p-value greater than p-crit to zero
#' @return r: correlation matrix
#' @return p: matrix of p-values
#' @return rlist: flattened list of correlations

library(Hmisc)
covariation = function( df, pcrit = NA ){
  
  tempdf = t( as.matrix(df) )
  res = rcorr( tempdf, type = 'spearman' )
  
  r = res$r
  p = res$P
  
  # In addition to returning correlation and p values, make a list of the correlations
  # Note that each pair appears twice, so use upper triangle
  upper = upper.tri(r)
  upperix = which(upper,arr.ind = TRUE)
  uppernames = data.frame(rownames(r)[upperix[,1]],colnames(r)[upperix[,2]])
  uppernames[,3] = paste(uppernames[,1],uppernames[,2],sep = ':')
  rownames(upperix) = uppernames[,3]
  
  if (!is.na(pcrit)){
    insig = which(p > pcrit,arr.ind = TRUE)
    r[insig] = 0
  }
  
  rlist = r[ upperix ]
  names(rlist) = rownames(upperix)
  
  returnlist = list(r,p,rlist)
  names(returnlist) = c('r','p','rlist')
  
  return( returnlist )
  
}