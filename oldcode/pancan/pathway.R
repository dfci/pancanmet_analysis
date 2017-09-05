# A script to simplify the KEGG pathway analysis for the Pancan Project. Remember to choose rings2 or rings3 (with/without annotations)

rm(list=ls())
setwd('/Users/ereznik/Documents/pancanmet/analysis/pancan/')

library(reshape)
library(ggplot2)
library(metabologram)
library(KEGGREST)
printpath = '/Users/ereznik/Documents/pancanmet/results/pathway_v2/'
datapath = '/Users/ereznik/Documents/pancanmet/data/pathview/'
source('../pathview/rings2.R')
source('../pathview/rings3.R')
source('../pathview/makemgram.R')

# Decide whether to make JPEGs or not
dojpeg = 0

# Do we save the data or not?
savedataparam = 1

# Decide whether to do all pathways (92) or not
maxpath = 92

################
# Get KEGG pathway names
keggdbpaths = keggList("pathway", "hsa") 
keggdbpaths = keggdbpaths[1:maxpath]
paths = sapply(strsplit(names(keggdbpaths[1:length(keggdbpaths)]),':'), "[", 2)
pathnames = sapply(strsplit(keggdbpaths[1:length(keggdbpaths)],' - Homo'), "[", 1)
keggnames = data.frame(pathnames); rownames(keggnames) = paths

# Read in hugo to entrez conversion
h2e = read.csv('/Users/ereznik/Documents/useful/hugo2entrez.txt',row.names = 1,header = T,sep = '\t',
               colClasses = 'character')

################
# Read the table of KEGG pathway annotations
keggdata = read.csv('../../data/pathway_v2/kegg_details_munged.csv',header = F,row.names = 1,colClasses='character')

studies = c('BLCA','BRCA','BRCATang','COAD','KICH','KIRC','OV','PAAD','PAADHussain1','PAADHussain2','PRAD','PRADLODA','STAD')
#studies = c('KIRC')

################
# Initialize some data arrays
pathmetx = matrix(0,dim(keggnames)[1],length(studies)); rownames(pathmetx) = keggnames[paths,1]; colnames(pathmetx) = studies
pathgenex = matrix(0,dim(keggnames)[1],length(studies)); rownames(pathgenex) = keggnames[paths,1]; colnames(pathgenex) = studies
pathmety = matrix(0,dim(keggnames)[1],length(studies)); rownames(pathmety) = keggnames[paths,1]; colnames(pathmety) = studies
pathgeney = matrix(0,dim(keggnames)[1],length(studies)); rownames(pathgeney) = keggnames[paths,1]; colnames(pathgeney) = studies
sizemet = matrix(0,dim(keggnames)[1],length(studies)); rownames(sizemet) = keggnames[paths,1]; colnames(sizemet) = studies
sizegene = matrix(0,dim(keggnames)[1],length(studies)); rownames(sizegene) = keggnames[paths,1]; colnames(sizegene) = studies

for (study in studies){
  allmgrams = list()
  mgramnames = c()
  
  if ( file.exists( paste(datapath,study,'_gene.csv',sep='') ) ){
    isg = TRUE
    g = read.csv(paste(datapath,study,'_gene.csv',sep=''),row.names = 1, header = F)
  } else{isg = FALSE}
  
  # Create the directories we are printing to
  outdirs = c(paste(printpath,study,sep=''),paste(printpath,study,'/metabolograms',sep=''))
  for (outdir in outdirs){
    do.call(file.remove,list(list.files(outdir,full.names=TRUE)))
    dir.create(file.path( outdir ), showWarnings = FALSE)
  }
  
  m = read.csv(paste(datapath,study, '_met.csv',sep=''),row.names = 1, header = F)
  
  for (path in paths){
    print(paste(study,path))
    
    setwd(paste(printpath,study,sep=''))
    if (isg){genes = g[,1]; names(genes) = rownames(g)} 
    mets = m[,1]; names(mets) = rownames(m)
    pathname = pathnames[paste('path',path,sep=':')]
     
    # Additionally, write the study and pathway specific genes and metabolites to a table
    metsinpath = unlist( strsplit(keggdata[path,2],'\\|') )
    intm =  intersect(metsinpath,rownames(m))
    if (dim(m)[2]==1){
      mdata = data.frame(m[intm,1]); rownames(mdata) = intm
    } else{mdata = m[intm,]}
    
    ################################
    # Make the variable "ringdata" 
    ################################
    if (isg){
      ginpath = unlist( strsplit(keggdata[path,1],'\\|') )
      entrezidx = which( h2e[,1] %in% ginpath )
      hugonames = h2e[entrezidx,1]
      entreznames = h2e[entrezidx,3]
      
      # Set the names of hugonames to entrez names to make the renaming easier
      names(hugonames) = entreznames
      
      intg = intersect(entreznames,rownames(g))
      if (dim(g)[2] == 1){
        gdata = data.frame( g[intg,1] ); rownames(gdata) = intg;
      } else{gdata = g[intg,]}
      rownames(gdata) = hugonames[intg]
      
      #write.csv(gdata,paste(printpath, study,paste(study,path,'gene.csv',sep='_'), sep='/' ) )
      
      numgene = length( which( abs( gdata[,1] ) >0 ) )/ length(which(!is.na(g[intg,1])))
      signgene = ( length( which( gdata[,1] >0 ) ) - length( which( gdata[,1] < 0 ) ) ) /
        length(which(!is.na(g[intg,1])))
      pathgenex[as.character(keggnames[path,1]),study] = numgene
      pathgeney[as.character(keggnames[path,1]),study] = signgene
      sizegene[as.character(keggnames[path,1]),study] = length(which(!is.na(g[intg,1])))
      
      # If we have gene data make a ring plot
      if (length(intm)>0){
        ringdata = list()
        ringdata[[1]] = data.frame(gdata[,1]); rownames(ringdata[[1]]) = rownames(gdata)
        ringdata[[2]] = data.frame(m[intm,1]); rownames(ringdata[[2]]) = mdata[,2]
      }else{
        ringdata = list()
        ringdata[[1]] = data.frame(gdata[,1]); rownames(ringdata[[1]]) = rownames(gdata)
        ringdata[[2]] = data.frame()
      }
      
    }
    
    if (!isg & length(intm)>0){
      # Make a ring plot using only metabolite data
      ringdata = list()
      ringdata[[1]] = data.frame()
      ringdata[[2]] = data.frame(m[intm,1]); rownames(ringdata[[2]]) = mdata[,2]
    }
    
    # Sort each variable in ringdata in alphabetical order
    rnames_sort = sort(rownames(ringdata[[1]]))
    temp = data.frame( ringdata[[1]][rnames_sort,] )
    rownames(temp) = rnames_sort
    ringdata[[1]] = temp
    
    rnames_sort2 = sort(rownames(ringdata[[2]]))
    temp = data.frame( ringdata[[2]][rnames_sort2,] )
    rownames(temp) = rnames_sort2
    ringdata[[2]] = temp
    
    ################################
    # END: Make the variable "ringdata" 
    ################################
    
    # Now that we have ringdata formatted correctly, plot away (rings3 adds annotation, rings2 does not)
  
    rings2(ringdata, pathname, paste(printpath,study,
      paste('metabolograms/',study,'_',gsub("[[:punct:][:space:]]", "", pathname),'_metabologram.pdf',sep=''), sep='/' ) )
    if (dojpeg == 1){# Repeat as JPEG
      rings2(ringdata, pathname, paste(printpath,study,
        paste('metabolograms/',study,'_',gsub("[[:punct:][:space:]]", "", pathname),'_metabologram.jpg',sep=''), sep='/' ) ) }
    
    # Save the data
    if (dim(mdata)[1] == 0){next}
    nummet =  length(which(abs(mdata[,1])>0)) / length(which(!is.na(mdata[,1])))
    signmet = ( length(which(mdata[,1]>0)) - length(which(mdata[,1]<0))  ) /
      length(which(!is.na(mdata[,1])))
    pathmetx[as.character(keggnames[path,1]),study] = nummet
    pathmety[as.character(keggnames[path,1]),study] = signmet
    sizemet[as.character(keggnames[path,1]),study] = length(which(!is.na(mdata[,1])))
   
    # Save data for shinyapp
    temp = makemgram( ringdata )
    if (length(allmgrams)==0){{allmgrams = list(temp)}}else{allmgrams[[length(allmgrams)+1]] = temp}
    mgramnames = c(mgramnames,pathname)
    
    #pdf(paste(printpath,study,paste('metabolograms/',study,
    #          '_',gsub("[[:punct:][:space:]]", "", pathname),'_metabologram_annotated.pdf',sep=''), sep='/' ))
    #print( metabologram(temp,width = 500,height = 500) )
    #dev.off()
  }

  names(allmgrams) = mgramnames
  if (study == 'KIRC'){
    # Save with this name for CCRCC project
    save(allmgrams,file = '/Users/ereznik/Documents/pancanmet/data/pathway_v2/allmgrams.RData')
  }
  mvarname = paste('allmgrams_',study,sep='')
  assign(mvarname,allmgrams)
  save(list = mvarname,file = paste('/Users/ereznik/Documents/pancanmet/data/pathway_v2/allmgrams_',study,'.RData',sep=''))
}


# Save the results, and add fail-safe to avoid overwriting accidentally
if (length(studies)>4 & savedataparam == 1){
  save(pathmetx,pathmety,pathgenex,pathgeney,sizemet,sizegene,file = '/Users/ereznik/Documents/pancanmet/data/pathway_v2/allpathwaydata.Rda')
}


