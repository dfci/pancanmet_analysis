# A function which adds a metabolite to the merged mapping dataframe and adds all of its IDs to the dictionaries.

# Issues:
#1. When we add to dictionaries here, we could possibly overwrite prior entries. 
#   I'm not really sure how to handle this conflict. Perhaps actually removing the entry and blacklisting, 
#   rather than overwriting?

newmet4map = function(dicts,mapping,nextmet,alternates,metname,study,rowsused){
  
  # Add this metabolite to the mapping data frame
  mapping[nextmet,study] = metname
  
  # For each ID type, add new data regarding where the identifier maps to.
  for (IDtype in names(alternates)){
    
    # Get the dictionary you want to use
    d2use = dicts[[IDtype]]
    
    # Remove NA and remove ''
    alternateIDs = na.omit( alternates[[IDtype]] )
    alternateIDs = setdiff(alternateIDs,'')
    if (length(alternateIDs) == 0){next}
    
    # Use only unique alternateIDs
    uq_alternateIDs = unique(alternateIDs)
    
    d2use[uq_alternateIDs,1] = nextmet
    
    # Now re-assign dictionary in dicts list
    dicts[[IDtype]] = d2use
    
  }
  
  # Update the rowsused
  rowsused = c(rowsused,nextmet)
  
  # Update the next metabolite
  metcounter = strsplit(nextmet,' ')[[1]][2]
  newmetcounter = as.integer( metcounter ) + 1
  nextmet = paste('Metabolite',newmetcounter)
  
  # Return
  return( list(dicts,mapping,nextmet,rowsused) )
  
}