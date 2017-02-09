# A function which adds a metabolite to the merged mapping dataframe and adds all of its IDs to the dictionaries.

# Issues
#1. When we add to dictionaries here, we could possibly overwrite prior entries.
#   I'm not really sure how to handle this conflict. Perhaps actually removing the entry and blacklisting,
#   rather than overwriting? CAN YOU FLAG (PRINT AN OUTPUT) OF WHEN THIS HAPPENS?

append2map = function(dicts,mapping,row2use,alternates,metname,study,rowsused,nextmet,nextFlag){

  if (row2use == ''){stop('We are using an empty rowname in append2map!')}

  # Add this metabolite to the mapping data frame
  mapping[row2use,study] = metname

  # For each ID type, add new data regarding where the identifier maps to.
  for (IDtype in names(alternates)){

    # Get the dictionary you want to use
    d2use = dicts[[IDtype]]

    # Remove NA and remove ''
    alternateIDs = na.omit( alternates[[IDtype]] )
    alternateIDs = setdiff(alternateIDs,'')
    if (length(alternateIDs) == 0){next}

    # Use only unique alternate IDs
    uq_alternateIDs = unique(alternateIDs)

    for (eachID in uq_alternateIDs){
      if (eachID %in% rownames(d2use)){

        # We've already seen this ID before, append the sources
        d2use[eachID,1] = row2use # PRINT SOMETHING HERE WHEN OVERRIDING HAPPENS!!!
        d2use[eachID,2] = paste(d2use[eachID,2],paste(study,metname,sep = ':'),sep=',')
      }else{

        # We've never seen this ID before
        d2use[eachID,1] = row2use
        d2use[eachID,2] = paste(study,metname,':')
      }

    }


    # Now reassign the dictionary in the dicts list
    dicts[[IDtype]] = d2use

  }

  # Update the rowsused
  rowsused = c(rowsused,row2use)

  if (nextFlag){
    # Update the next metabolite
    metcounter = strsplit(nextmet,' ')[[1]][2]
    newmetcounter = as.integer( metcounter ) + 1
    nextmet = paste('Metabolite',newmetcounter)
  }


  # Return
  return( list(dicts,mapping,nextmet,rowsused) )

}
