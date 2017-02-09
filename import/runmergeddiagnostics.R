# A function to evaluate whether the output of 

# Issues:
# 1. Minor issue, but since we never renamed the rownames of the mergedIDs file, 
#   we do it here by using the ronwames of the merged metabolomics file. 
#   This should be fine, but would be better addressed earlier in the import scripts.

library(data.table)
source('useful_metimport.R')

runmergeddiagnostics = function( studies,logfile ){
  
  # 1. Print a file indicating the different names aligned to each row of a metabolite, and their counts
  ids = read.csv('tempdir/merged_mapping.csv',header = TRUE,row.names = 1,stringsAsFactors = FALSE)
 
  # Create a file to write to
  diagfile1 = 'tempdir/uniquemets.csv'
  write('NULL HEADER,NULL HEADER,NULL HEADER,NULL HEADER,NULL HEADER ',diagfile1,append = TRUE)
  
  for (rname in rownames(ids)){
    idtable = table(t(ids[rname,]))
    
    # Make a string
    writelist = sapply(1:length(idtable),function(x){paste(names(idtable)[x],idtable[x],sep=':')})
    writelist = c(rname,length(idtable),writelist)
    writeline = paste(writelist,collapse = ',')
    write(writeline,diagfile1,append = TRUE)
  }
  
  # 2. Confirm that we can "unmap" the data
  alldata = getbigdata('tempdir/merged_metabolomics.csv')
  
  for (study in studies){
    # Read in the original data file
    metfile = paste('tempdir/',study,'_met.csv',sep='')
    met = read.csv(metfile,header = FALSE,skip = 2,row.names = 1)
    
    sampname = scan(metfile, nlines = 1, what = character(),sep=',')
    samptype = scan(metfile, skip = 1, nlines = 1, what = character(),sep=',')
    metcnames = sapply(2:length(sampname),function(x){paste(study,sampname[x],simpleCap(samptype[x]),sep=':')}) # Skip first element, which is blank
    colnames(met) = metcnames 
    
    # Get the relevant column in ids
    studyids = ids[,study]
    names(studyids) = rownames(alldata) # This is the cause of issue 1 in header
    
    # Drop NAs
    studyids = studyids[which(!is.na(studyids))]
    
    # "Undo" the merge mapping
    undomet = alldata[names(studyids),colnames(met)]
    
    # Rename rows of undomet to original names
    rownames(undomet) = studyids
    
    # Check that columnnames and rownames match
    diffrow = which(!(rownames(undomet)%in%rownames(met)))
    diffcol = which(!(colnames(undomet)%in%colnames(met)))
    
    if (length(diffrow) !=0 | length(diffcol) != 0){
      stop(paste('Error in post-merge diagnostics when "unmerging" data:',study))
    }
    diffmet = met - undomet[rownames(met),]
    diffval = max(abs(diffmet))
   
    if (diffval > 1e-10){
      print(paste(study,diffval,sep=':'))
      stop(paste('Error in post-merge diagnostics when "unmerging" data:',study))   
    }
    
  }
  
  
}