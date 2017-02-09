# Script to permute differential abundance results and determine cross-cancer significance

rm(list = ls())
setwd('/Users/ereznik/pancanmet_analysis/analysis/')
library(ggplot2)

findrowname = function(name,df){
  ix = which(rownames(df) == name)
}

# Set number of permutations
numperm = 2000
d = read.csv('../results/DifferentialAbundance_SignedChanges.csv',header = TRUE,row.names = 1)
res = matrix(0,dim(d)[2]+1,numperm) # maximum value is the number of studies
rownames(res) = paste('MAX',0:dim(d)[2],sep = ':')

respos = res
resneg = res

for (i in 1:numperm){
  
  if (i%%100 == 0){print(i)}
  dperm = matrix(0,dim(d)[1],dim(d)[2])
  
  # For each study, permute the non-na values
  for (col in 1:dim(d)[2]){
    notna = which(!(is.na(d[,col])))
    shuffled = sample(notna,length(notna))
    
    # Replace values
    dperm[notna,col] = d[shuffled,col]
    
  }
  
  # Calculate sums
  rsums_abs = max( rowSums( abs( dperm ) ) )
  maxix = findrowname(paste('MAX:',rsums_abs,sep =''),res)
  res[1:maxix,i] = 1
  
  dpermpos = dperm
  dpermpos[dpermpos<1] = 0
  rsums_pos = max( rowSums( dpermpos ) )
  maxix = findrowname(paste('MAX:',rsums_pos,sep =''),respos)
  respos[1:maxix,i] = 1
  
  dpermneg = dperm
  dpermneg[dpermneg > -1] = 0
  rsums_neg  = abs( min( rowSums( dpermneg ) ) )
  maxix = findrowname(paste('MAX:',rsums_neg,sep =''),resneg)
  resneg[1:maxix,i] = 1
}

print(rowSums(res)/numperm)
print(rowSums(respos)/numperm)
print(rowSums(resneg)/numperm)
