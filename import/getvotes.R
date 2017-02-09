# A function which gets the votes for the metabolite mapping

getvotes = function( metIDs,dicts ){
  
  # Initialize some variables for this metabolite
  votes = c() # A list of which rows each dictionary indicates we should map to
  alternates = list() # A list of alternate identifiers for each type (e.g. alternate Chebi IDs)
  
  for (IDtype in names(metIDs)){
    
    # Get the actual string of IDs
    idstring = metIDs[[IDtype]]
    
    # Separate out each alternative identifier
    if (grepl('\\|',metIDs[IDtype])){
      alternates[[IDtype]] = unlist(strsplit(idstring,'\\|'))
    }else{
      alternates[[IDtype]] = idstring
    }
    
    # Get the full list of alternate names
    separatenames = alternates[[IDtype]]
    
    # To further prevent weird behavior in R-dataframes rownames, manually search for which names
    trimnames = separatenames[separatenames %in% rownames(dicts[[IDtype]])]
    
    newvotes = dicts[[IDtype]][trimnames,1]
    newvotes = na.omit( newvotes )
    votes = c(votes,newvotes)
    
  }
  
  return(list(votes,alternates))
  
}