# Function to identify metabolomics studies to import and confirm that they have necessary files
library(data.table)
library(stringr)

# Issues:
# 1. Need to add code for handling instance when no ID file is present (AL)
# 2. Check that column names of the IDs file are in the appropriate order (AL)

# Description
# IDdir: path to directory containing subdirectories for each study.
  # In those subdirectories, there should be an metabolite ID file.
  # Naming convention: *_FinalMetIDs.csv, where * corresponds to a study name
# metdir: path to directory containing subdirectories for each study.
  # In those subdirectories, there should be a metabolomics data matrix.
  # Naming convention: *_metdata.csv, where * corresponds to study name
# logfile: a time-stamped text file containing a log of all steps in the import process.

getstudies = function( IDdir,metdir,logfile ) {

  studies = list.dirs( IDdir, recursive = FALSE, full.names = FALSE )

  # Remove the directory 'oldstudies_EXCLUDED', this is for storing old data that might be included later
  studies = setdiff(studies, c('COAD','STAD','oldstudies_EXCLUDED') )

  for (study in studies){

    # 1. Check that metabolomics file exists, get rownames for metabolites
    metfile = paste(metdir,'/',study,'/',study, '_metdata.csv',sep='')
    if (file.exists( metfile )){
      tempmet = fread(metfile,sep = ',',header = TRUE)
      met = data.frame( tempmet, row.names = 1, check.names = FALSE )

      # Drop the first row of the metabolomics data, as this is an indicator of tumor/normal status
      metnames = sapply( rownames(met)[2:dim(met)[1]], tolower )

      # Trim whitespace too
      metnames = sapply( metnames, str_trim)
      
      ##Added July 2016: Change the tumor/normal status to lower case
      tnstatus = sapply(met[1,], tolower)
      met[1,] = tnstatus

    }else{
        stop(paste('Metabolomics file does not exist for study:',study))
    }

    # 2. Check that final metabolite IDs file exists, and that rownames match those for metabolomics data
    IDfile = paste(IDdir,'/',study,'/',study,'_FinalMetIDs.csv',sep='')
    if (file.exists(IDfile)){

      # Check that metabolite names match, change to lowercase and trim whitepsace
      IDs = read.csv(IDfile,header = TRUE,row.names = 1,na.strings = c('NA','nan'))
      IDnames = sapply(rownames(IDs),tolower)
      IDnames = sapply(IDnames,str_trim)

      diffnames = which(IDnames != metnames)

      if (length( diffnames )!=0){
        write(paste('The following differences are present in ID file rownames and metabolomics file rownames for',study,':'),
              logfile,append = TRUE)
        for (ii in diffnames){
          write( paste(IDnames[ii],metnames[ii] ),logfile,append = TRUE )
        }
        write('\n',logfile,append = TRUE)
        #stop(paste('Rownames in metabolomics and metabolite ID files do not match for study:',study))

        # Change rownames of met to IDnames for consistency and write to tempdir
        write(paste('Changing rownames in metabolomics file to ID names in',study,'.\n'),logfile,append = TRUE)
      }

    }else{
      # We need to put something here to autmoatically run ID conversion
      stop(paste('No ID file present for',study))
      # CREATE ID FILE HERE
      # ORDER MATTERS ChemicalName	PubChem CID	KEGG	Human Metabolome Database	ChEBI, TRY TO USE NAME
    }

    # 3. Check that column names are correct, and if not, change them
    id_cnames = colnames(IDs)
    standard_cnames = c('PubChem.CID','KEGG','Human.Metabolome.Database','ChEBI' )
    diff_cnames = which(id_cnames != standard_cnames)
    if (length(diff_cnames)!=0){
      write(paste('WARNING: Changing column names in metabolite IDs file to match nomenclature for',study,'.\n'),logfile,append = TRUE)
      colnames(IDs) = standard_cnames
    }

    # 4. Write ID file and metabolomics file to temporary directory
    # Use pretty version of IDnames
    rownames(IDs) = IDnames
    rownames(met)[2:dim(met)[1]] = IDnames
    write.csv(IDs,paste('tempdir/',study,'_IDs.csv',sep=''))
    write.csv(met,paste('tempdir/',study,'_met.csv',sep=''))

    # 5. Write to logfile
    writeline = paste('Done checking for presence of ID and metabolomics files for ',study,'.\n')
    write(writeline,logfile,append = TRUE)

  }

  return(studies)
}
