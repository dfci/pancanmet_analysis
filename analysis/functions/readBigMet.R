#' Function to read big merged metabolomics data file
#' 
#' @param study2drop: list of studies to remove from the metabolomics file
#' @return met: a matrix of metabolite values for each sample
#' @return studytype: a list of the originating study for each sample/column in met
#' @return tissuetype: a list of the kind of sample (tumor/normal/metastasis) for each sample/column in met
#' 
readBigMet = function( study2drop ){
  met = read.csv('../data/merged_metabolomics/merged_metabolomics.csv', row.names = 1,check.names = FALSE,header = TRUE)
  splittypes = strsplit(colnames(met),'\\:')
  studytype = sapply(splittypes, "[", 1)
  names(studytype) = colnames(met)
  tissuetype = sapply(splittypes, "[", 3)
  names(tissuetype) = colnames(met)
  
  if (length(study2drop)!= 0){
    
    # Drop data corresponding to these studies
    dropix = which(studytype %in% study2drop)
    if (length(dropix)>0){
      met = met[,-dropix]
      studytype = studytype[-dropix]
      tissuetype = tissuetype[-dropix]
      
      # Check that we haven't messed up names
      if (length(which(colnames(met)!=names(studytype))) >0 | length(which(colnames(met)!=names(tissuetype))) >0){
        stop()
      }
    }

  }
  
  res = list(met,studytype,tissuetype)
  names(res) = c('met','studytype','tissuetype')
  return(res)
}
