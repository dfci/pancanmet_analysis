# A function to make data compatible for the shiny app
library(RColorBrewer)
## color and breaks 
minval = -3
maxval = 3
br = seq( minval, maxval,.25)
cl = colorRampPalette(c("blue","white","red"))(length(br)-1)
save(file = '/Users/ereznik/Documents/pancanmet/data/pathway_v2/allmgrams_colorbar.Rda',br,cl)


makemgram = function( ringdata ){
  
  # Flip ringdata around
  tempring = list(ringdata[[2]],ringdata[[1]])
  ringdata = tempring
  
  mfc = mean(ringdata[[1]][,1])
  cidx = which(br>mfc)
  if (length(cidx) <= 1){cidx = length(br) - 1 
  }else{cidx = min(cidx)}
  if (mfc == 0){mcolour = '#FFFFFF'}else{mcolour = cl[cidx]}
  
  if (length(ringdata[[2]])>0){
    
    gfc = mean(ringdata[[2]][,1])
    cidx = which(br>gfc)
    if (length(cidx) <= 1){cidx = length(br) - 1 
    }else{cidx = min(cidx)}
    if (gfc == 0){gcolour = '#FFFFFF'}else{gcolour = cl[cidx]}
    
    mgram = list( list(name = 'Metabolites', 
                       NumberIncreased = length(which(ringdata[[1]][,1] >0)),
                       NumberDecreased = length(which(ringdata[[1]][,1] <0)),
                       colour = mcolour,
                       children = list()), 
                  list(name = 'Genes', 
                       NumberIncreased = length(which(ringdata[[2]][,1] >0)),
                       NumberDecreased = length(which(ringdata[[2]][,1] <0)),
                       colour = gcolour,
                       children = list() ))
  }else{
    mgram = list( list(name = 'Metabolites', 
                       NumberIncreased = length(which(ringdata[[1]][,1] >0)),
                       NumberDecreased = length(which(ringdata[[1]][,1] <0)),
                       colour = mcolour,
                       children = list()), 
                  list(name = 'Genes', children = list() ))
  }
  
 
  
  # Do metabolites first
  for (i in 1:dim( ringdata[[1]] )[1] ){
    if (dim(ringdata[[1]])[1] == 0){
      mgram[[1]]$children[[1]] = list(name = 'NA', Fold = 0, colour = '#6b6b6b')
      break
    }
    
    fc = ringdata[[1]][i,1]
    cidx = which(br>fc)
    if (length(cidx) <= 1){cidx = length(br) - 1 
    }else{cidx = min(cidx)}
    if(is.na(fc)){colour = '#6b6b6b'}else if (fc == 0){colour = '#FFFFFF'}else{colour = cl[cidx]}
    
    mgram[[1]]$children[[i]] = list(name = rownames(ringdata[[1]])[i], Fold = round(ringdata[[1]][i,1],3), colour = colour)
  }
  
  for (i in 1:dim( ringdata[[2]] )[1] ){
    if (dim(ringdata[[2]])[1] == 0){
      mgram[[2]]$children[[1]] = list(name = 'NA', Fold = 0, colour = '#6b6b6b')
      break
    }
    
    fc = ringdata[[2]][i,1]
    cidx = which(br>fc)
    if (length(cidx) <= 1){cidx = length(br) - 1 }else{cidx = min(cidx)}
    if(is.na(fc)){colour = '#6b6b6b'}else if (fc == 0){colour = '#FFFFFF'}else{colour = cl[cidx]}
    mgram[[2]]$children[[i]] = list(name = rownames(ringdata[[2]])[i], Fold = round( ringdata[[2]][i,1],3), colour = colour)
  }
  
  return( mgram )
}