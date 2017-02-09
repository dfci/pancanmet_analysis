#' Function to map metabolites to KEGG pathways
#' @param logfile: filename of the logfile we are using to keep track of what we do
#' 
#' 

map2KEGG = function( logfile ){
  
  # The default location for the KEGG data is in a dictionary in tempdir/ and a list of lists in ../data/
  keggdict = read.csv('tempdir/KEGG_dictionary.csv',header = TRUE,row.names = 1,stringsAsFactors = FALSE)
  colnames(keggdict) = c('CompoundID','MetaboliteRow','Source','FinalName')
  keggpathways = load('../data/KEGGpathays.RData')
  
  # Because of conflicts, there is sometimes a KEGG ID which maps to two different metabolites, get rid of this
  dupix = which(duplicated(keggdict$CompoundID))
  dupnames = keggdict[dupix,'CompoundID']
  if (length(dupix)!=0){
    
    # Remove these duplicates
    pastednames = paste(dupnames,collapse = ',')
    writeline = paste('Removing the following metabolites from the KEGG dictionary:',pastednames,sep='')
    write(writeline,logfile,append = TRUE)
    
    ix2remove = which(keggdict$CompoundID %in% dupnames)
    keggdict = keggdict[-ix2remove,]
  }
  
  rownames(keggdict) = keggdict$CompoundID
  
  # Make a big matrix of pathways x metabolites
  uniquemets = unique(keggdict$FinalName)
  pathwaymap = matrix(0,length(cpdmap),length(uniquemets))
  rownames(pathwaymap) = names(cpdmap)
  colnames(pathwaymap) = uniquemets
  
  # For each pathway, add a 1 for each metabolite in the pathway to pathwaymap
  for (pathway in names(cpdmap)){
    
    compounds = cpdmap[[pathway]]
    
    # Find compounds in dictionary
    ixcpds = intersect(names(compounds),rownames(keggdict))
    
    # Find the corresponding name of the metabolite in merged metabolomics data using kegg dictionary
    cpdsInMerged = keggdict[ixcpds,'FinalName']
    
    # Add a 1 to these compounds
    pathwaymap[pathway, cpdsInMerged] = 1
    
    # Write to logfile
    cpdString = paste(cpdsInMerged,collapse = ',')
    writeline = paste('The following metabolites are in',pathway,':',cpdString,'\n')
    write(writeline,logfile,append = TRUE)
  }
  
  # Rename rownames of pathwaymap
  rownames(pathwaymap) = sapply(rownames(pathwaymap),function(x){y = strsplit(x,'Homo')[[1]][1]; z = substr(y,1,nchar(y)-3)})
  
  write.csv(pathwaymap,'tempdir/KEGG_pathwaymap.csv')
}

