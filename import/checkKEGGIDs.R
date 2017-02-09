# A function which validates the KEGG dictionary that we have imported from our merged data

# NB for Ed: remove metabolites which have no pathways they map to

checkKEGGIDs = function( logfile ){
  
  library(stringdist)
  library(KEGGREST)
  
  # The default location for the KEGG data is in a dictionary in tempdir/
  keggdict = read.csv('tempdir/KEGG_dictionary.csv',header = TRUE,row.names = 1,stringsAsFactors = FALSE)
  colnames(keggdict) = c('CompoundID','MetaboliteRow','Source','FinalName')
  
  # Because of conflicts, there is sometimes a KEGG ID which maps to two different metabolites, get rid of this
  dupix = which(duplicated(keggdict$CompoundID))
  dupnames = keggdict[dupix,'CompoundID']
  if (length(dupix)!=0){
    
    # Remove these duplicates
    pastednames = paste(dupnames,collapse = ',')
    writeline = paste('Removing the following metabolites from the KEGG dictionary:',pastednames,sep='')
    write(writeline,logfile,append = TRUE)
    
    ix2remove = which(keggdict$CompoundID %in% dupnames)
    if (length(ix2remove) == 0){stop('Error in removing duplicate IDs from KEGG dictionary.')}
    keggdict = keggdict[-ix2remove,]
  }
  
  # Set rownames of KEGG dictionary as compound IDs
  rownames(keggdict) = keggdict$CompoundID
  
  # Initialize a dataframe for storing distances, separately tracking compound IDs and drug IDs
  cdist = data.frame( 'Distance' = numeric(0), 'KEGGReference' = character(0),stringsAsFactors = FALSE)
  ddist = data.frame( 'Distance' = numeric(0), 'KEGGReference' = character(0),stringsAsFactors = FALSE)
  
  # For each row in keggdict, compute the string distance to the actual KEGG entry. Go 10 by 10 to reduce lag time in calling keggGet
  for (ctr in seq(1,dim(keggdict)[1],10)){
    cIDs2use = rownames(keggdict)[ctr:(ctr+9)]
    
    # Get KEGG data
    keggmet_list = keggGet(cIDs2use)
    
    # Name the keggmet_list
    names(keggmet_list) = sapply(keggmet_list,function(x){x$ENTRY})
    
    # Find the names which are not in the list, these couldn't be mapped, and write them to log file
    notinlist = setdiff(cIDs2use,names(keggmet_list))
    if (length(notinlist)!=0){
      notinlist = paste(notinlist,collapse = ',')
      write(paste('These metabolites could not be checked for distance to KEGG:',notinlist),logfile,append = TRUE)
    }
    
    # Now go through this list
    for (cID in names(keggmet_list)){
      
      metname = keggdict[ cID,'FinalName' ]
      print(metname)
      
      # If this metabolite does not map to a pathway, continue
      if ( !('PATHWAY' %in% names(keggmet_list[[cID]]) ) ){next}
      
      keggmet = keggmet_list[[cID]]$NAME
      
      # Make sure everything is lower case and strip semicolons
      keggmet = sapply(keggmet,tolower)
      keggmet = sapply(keggmet,function(x){gsub('[;]','',x)})
      
      # Strip d- and l- prefixes
      keggmet = sapply(keggmet,function(x){gsub('^[dl]-','',x)})
      metname = sapply(metname,function(x){gsub('^[dl]-','',x)})
      
      metname = tolower(metname)
      
      # Remove parenthesis content
      metname = sapply(metname,function(x){gsub(' \\(.*\\)','',x)})
      
      # Replaces dashes/spaces, remove asterisks
      keggmet = sapply(keggmet,function(x){gsub(' ','',x)})
      keggmet = sapply(keggmet,function(x){gsub('-','',x)})
      
      metname = sapply(metname,function(x){gsub(' ','',x)})
      metname = sapply(metname,function(x){gsub('-','',x)})
      metname = sapply(metname,function(x){gsub('\\*','',x)})
      
      # Calculate distances
      tempdist = stringdist(metname,keggmet)
      
      # Determine if drug or compound
      if (substr(cID,1,1)=='D'){
        ddist[paste(cID,metname,sep=':'),'Distance'] = min( tempdist )
        ddist[paste(cID,metname,sep=':'),'KEGGReference'] = paste(keggmet,collapse = ':')
      }else{
        cdist[paste(cID,metname,sep=':'),'Distance'] = min( tempdist )
        cdist[paste(cID,metname,sep=':'),'KEGGReference'] = paste(keggmet,collapse = ':')
      }
      
    }
    
    
  }
  
  # Write the files to tempdir
  write.csv(cdist,'tempdir/KEGGID_Check_Compounds.csv')
  write.csv(ddist,'tempdir/KEGGID_Check_Drugs.csv')
}