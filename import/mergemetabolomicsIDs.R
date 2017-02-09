# Function to merge metabolomics ID files
library(data.table)
library(stringr)
source('useful_metimport.R')
source('append2map.R')
source('getvotes.R')

# Issues:
# 1. The case (upper/lower) of the metabolite name in config ID is not controlled, could cause things to be missed.  CURRENTLY USING tolower(). 

mergemetabolomicsIDs = function( studies, logfile ){

  # Read in the config file
  configIDs = read.csv('configIDs.csv',header = TRUE)
  rownames(configIDs) = sapply(1:dim(configIDs)[1],function(x){paste(configIDs[x,1],configIDs[x,2],sep=':')})

  # Initialize a bunch of dataframes to use as dictionaries. Dataframes preferred so we can't have duplicate keys.
  nulldf = data.frame('ID' = character(),'Source' = character(),stringsAsFactors=FALSE)
  dicts = list( nulldf,nulldf,nulldf,nulldf,nulldf )
  names(dicts) = c('Name','PubChem.CID','KEGG','Human.Metabolome.Database','ChEBI' )

  # Initialize a dataframe for storing the mappings of every metabolite
  mapping = data.frame( matrix(NA,0,length(studies)) )
  colnames(mapping) = studies

  # Set the current metabolite counter
  nextmet = 'Metabolite 1'

  for (study in studies){

    # Keep a running tally of the rows which have already been mapped to for a study.
    # Don't allow overwriting, force a new row to be created.
    rowsused = c()

    # Read in the ID file
    ids = read.csv( paste('tempdir/',study,'_IDs.csv',sep=''),header = TRUE,row.names = 1,stringsAsFactors=FALSE )

    # For each field in the file plus the name, check if it's in a dictionary
    for (metname in rownames(ids)){

      # Check if this metabolite is in configIDs file, and if so, add a new metabolite and move onto next one.
      mergename = paste(metname,study,sep=':')

      if (tolower(mergename) %in% tolower(rownames(configIDs))){
        configix = which(rownames(configIDs) == mergename) # Get the index
        alternates = c()
        returnlist = append2map(dicts,mapping,nextmet,alternates,metname,study,rowsused,nextmet,TRUE)
        if (configIDs[configix,'AddtoDict']){
          dicts = returnlist[[1]]; mapping = returnlist[[2]]; nextmet = returnlist[[3]]; rowsused = returnlist[[4]]
        }else{
          # Don't update the dictionaries, update everything else though
          mapping = returnlist[[2]]; nextmet = returnlist[[3]]; rowsused = returnlist[[4]]
        }

        # Move on to next metabolite
        next
      }

      # Get the IDs for this metabolite
      metIDs = c(metname,ids[metname,])
      names(metIDs) = c('Name',colnames(ids))

      # Get votes based on these ids
      returnvotes = getvotes( metIDs,dicts )
      votes = returnvotes[[1]]
      alternates = returnvotes[[2]]

      # If we have never seen this metabolite before, then add a new entry
      if (length(votes) == 0){

        # We are going to map this metabolite to nextmet
        returnlist = append2map(dicts,mapping,nextmet,alternates,metname,study,rowsused,nextmet,TRUE)
        dicts = returnlist[[1]]; mapping = returnlist[[2]]; nextmet = returnlist[[3]]; rowsused = returnlist[[4]]

        next
      }

      # Check if name already in name dictionary, if so, map there by default (names overrule other IDs)
      if (metname %in% rownames(dicts$Name)){

        # Get the row to map to
        row2use = dicts$Name[metname,1]

        # Check if we have already mapped to this row, if so, just make a new row
        if (row2use %in% rowsused){
          # We have already mapped to this metabolite once in this study, make a new row
          returnlist = append2map(dicts,mapping,row2use,alternates,metname,study,rowsused,nextmet,TRUE)
        }else{
          # Map to this row
          returnlist = append2map(dicts,mapping,row2use,alternates,metname,study,rowsused,nextmet,FALSE)
        }

        # Unpack the returned values
        dicts = returnlist[[1]]; mapping = returnlist[[2]]; nextmet = returnlist[[3]]; rowsused = returnlist[[4]]

        # Move on to next row
        next
      } # end if metname %in% rownames(dicts$Name)


      # If we have seen this metabolite before, check the consensus
      if (length(unique(votes)) == 1){

        # We have consensus
        row2use = unique(votes)

        # Check if we have already mapped to this row, if so, just make a new row
        if (row2use %in% rowsused){

          # We have already mapped to this metabolite once in this study, make a new row
          returnlist = append2map(dicts,mapping,nextmet,alternates,metname,study,rowsused,nextmet,TRUE)

        }else{

          # Map to this row
          returnlist = append2map(dicts,mapping,row2use,alternates,metname,study,rowsused,nextmet,FALSE)

        }

        # Unpack the returned values
        dicts = returnlist[[1]]; mapping = returnlist[[2]]; nextmet = returnlist[[3]]; rowsused = returnlist[[4]]

      }else{

        # Use the plurality, and if there is a tie, use first entry, write to the logfile
        trimvotes = votes[ which(!(votes %in% rowsused)) ]
        topvote = getmode( trimvotes )
        writeline = paste('Using plurality to map metabolite',metname,'in',study,'with votes:')
        writeline2 = votes
        write(writeline,logfile,append = TRUE)
        write(writeline2,logfile,append = TRUE)

        # If there is a tie, note this and choose first among topvotes
        if (length(topvote) != 1){
          writeline = paste('There was a tie mapping metabolite',metname,'in',study,'!')
          write(writeline,logfile,append = TRUE)
          topvote = topvote[1]

        }

        # Add to mapping
        if (is.na(topvote)){
          # We didn't have any places to map that weren't already mapped for this study, make a new row
          returnlist = append2map(dicts,mapping,nextmet,alternates,metname,study,rowsused,nextmet,TRUE)
        }else{
          returnlist = append2map(dicts,mapping,topvote,alternates,metname,study,rowsused,nextmet,FALSE)
        }

        dicts = returnlist[[1]]; mapping = returnlist[[2]]; nextmet = returnlist[[3]]; rowsused = returnlist[[4]]

      } #end if/else length(votes) == 1

    } #end for (metname in rownames(ids))

  } #end for (study in studies)

  # Write results to files
  write.csv(mapping,'tempdir/merged_mapping.csv')

  for (i in 1:length(dicts)){
    write.csv(dicts[[i]],paste('tempdir/',names(dicts)[[i]],'_dictionary.csv',sep=''))
  }
}
