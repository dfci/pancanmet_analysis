# Function to merge metabolomics data files.
library(data.table)
library(stringr)
source('useful_metimport.R')
source('append2map.R')

# Issues:
# 1. When making the merged data frame, the way I am subsetting the data (fullmapping,mapping) scares me, re-write? (PROBABLY OK WITH DIAGNOSTICS)

mergemetabolomicsdata = function( logfile ){

  # Read the merged IDs file
  ids = read.csv('tempdir/merged_mapping.csv',header = TRUE,
                    row.names = 1,stringsAsFactors=FALSE)

  # Set the rownames of each metabolite ID to the mode of the row
  metnames = data.frame(rowname = character(0),stringsAsFactors = FALSE)
  for (rname in rownames(ids)){
    altnames = ids[rname,]
    altnames = altnames[!is.na(altnames)]
    topname = getmode( altnames )[[1]] # Just choose the first name
    metnames[rname,1] = topname
  }

  # Check that metnames are all unique and that there are no NAs
  isdup = which(duplicated(metnames))
  NAcheck = which( !complete.cases(metnames) )
  if (length(isdup) !=0 | length(NAcheck) != 0){
    stop('There is an error assigning metabolite names to the merged data file in mergemetabolomics.R')
  }else{
    rownames(ids) = metnames[,1]

    # Write this to a file
    write.csv(metnames,'tempdir/FinalMetaboliteNames.csv')

    # Update each dictionary
    dictnames = list.files('tempdir/',pattern = 'dictionary')
    for (dictname in dictnames){
      write(paste('Updating ',dictname,'with final metabolite names...'),logfile,append = TRUE)
      dict = read.csv(paste('tempdir/',dictname,sep=''),header = TRUE,stringsAsFactors = FALSE)
      dict$FinalName = metnames[ dict$ID, 1]
      write.csv(dict, paste('tempdir/',dictname,sep=''))
    }
  }

  # Get unique studies
  studies = colnames(ids)

  # Set flag to initialize data frame
  firstflag = TRUE

  # Read in studies one by one
  for (study in studies){
    fname = paste('tempdir/',study,'_met.csv',sep='')
    tempmet = read.csv(fname,sep = ',',header = FALSE,skip = 2,row.names = 1)

    sampname = scan(fname, nlines = 1, what = character(),sep=',')
    samptype = scan(fname, skip = 1, nlines = 1, what = character(),sep=',')

    # Make sure below that the third field, indicating tumor/normal/metastasis status, is in sentence case
    metcnames = sapply(2:length(sampname),function(x){paste(study,sampname[x],simpleCap(samptype[x]),sep=':')})

    met = data.frame( tempmet )
    colnames(met) = metcnames

    if (firstflag){

      merged = data.frame( matrix(NA,dim(ids)[1],dim(met)[2]) )
      rownames(merged) = rownames(ids)
      colnames(merged) = colnames(met)
    }

    # Get the correct mapping of IDs and drop NAs
    fullmapping = data.frame( ids[,study],stringsAsFactors = FALSE )
    rownames(fullmapping) = rownames(ids)
    ## Drop NA columns
    notNAix = which(complete.cases(fullmapping))
    mapping = fullmapping[notNAix,]
    names(mapping) = rownames(fullmapping)[notNAix]

    # Add
    if (firstflag){
      merged[names(mapping),colnames(met)] = met[mapping,colnames(met)]
      firstflag = FALSE
    }else{
      d2merge = data.frame( matrix(NA,dim(ids)[1],dim(met)[2]) )
      rownames(d2merge) = rownames(ids)
      colnames(d2merge) = colnames(met)

      d2merge[names(mapping),colnames(met)] = met[mapping,colnames(met)]

      # Merge with merged
      merged = cbind(merged,d2merge)
    }

  }

  write.csv(merged,'tempdir/merged_metabolomics.csv')

}
