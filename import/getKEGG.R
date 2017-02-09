# A script to download and arrange data on KEGG pathways

getKEGG = function(){
  
  library(KEGGREST)
  
  # Get list of human pathways
  hsapathways = keggList('pathway','hsa')
  
  # Make a list, and get compound IDs and gene IDs for each path. Ultimately, we will have a list of lists
  cpdmap = list()
  genemap = list()
  
  for (pathway in names(hsapathways)){
    print(pathway)
    pathID = hsapathways[[pathway]]
    pathdata = keggGet(pathway)
    
    if (length(pathdata)!=1){
      stop( paste('Multiple values returned for pathway',pathway)) 
    }else{
        pathdata = pathdata[[1]]
    }
    
    if ('COMPOUND' %in% names(pathdata)){
      cpdmap[[pathID]] = pathdata$COMPOUND
    }else{
      cpdmap[[pathID]] = c()
    }
    
    
    # For genes, we have to take every other entry so that we get the actual gene names
    if ('GENE' %in% names(pathdata)){
      genemap[[pathID]] = pathdata$GENE[seq(2,length(pathdata$GENE),2)]
    }else{
      genemap[[pathID]] = c()
    }
    
    
  }
  
  # Save the data
  save(list = c('cpdmap','genemap','hsapathways'),file = '../data/KEGGpathays.RData')
  
}

