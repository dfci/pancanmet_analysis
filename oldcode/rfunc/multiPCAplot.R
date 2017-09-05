# A function which makes multiplots of PCA using different principal components. Only works for 2-4 plots!
setwd('/Users/ereznik/Documents/pancanmet/analysis/rfunc/')
multiPCAplot = function(pca,dim2plot,dtype,fname){
  library(ggplot2)
  source('multiplot.R')
  
  # PCA is the object returned from the prcomp function
  # dim2plot is a list pairs of indices, indicating the PCA dimensions to plot
  # dtype is a list indicating the "type" of sample, e.g. tumor or normal or metastasis
  # fname is the filename to print the resulting plot to
  
  pctr = 1
  pnames = c()
  for (dims in dim2plot){
    dim1 = dims[1]
    dim2 = dims[2]
    
    pc1score = round(100*pca$sdev[dim1]^2/sum(pca$sdev^2),3)
    pc2score = round(100*pca$sdev[dim2]^2/sum(pca$sdev^2),3)
    
    pcaplot = data.frame(pca$rotation[,c(dim1,dim2)])
    colnames(pcaplot) = c('PC1','PC2')
    pcaplot$Type = dtype
    
    p = ggplot(pcaplot,aes(PC1,PC2,color = Type)) + geom_point() + theme_bw() + 
      xlab(paste('PC',dim1,': % Variance Explained:',pc1score)) + 
      ylab(paste('PC',dim2,': % Variance Explained:',pc2score))
    
    pname = paste('p',pctr,sep='')
    pctr = pctr + 1
    pnames = c(pnames,pname)
    assign(pname,p)
    print( p )
  }
  
  # Make a multiplot
  pnames = unlist(pnames)
  pdf(fname,height = 5,width = 10)
  if (length(dim2plot) == 2){
    multiplot(p1,p2,cols = 2)
  }
  if (length(dim2plot) == 3){
    multiplot(p1,p2,p3,cols = 2)
  }
  if (length(dim2plot) == 4){
    multiplot(p1,p2,p3,p4,cols = 2)
  }
  dev.off()
  
}