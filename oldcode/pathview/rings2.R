# Function to make donut ring plot in ggplot

rings2 = function( ringdata,pathway,path2save ){
  #print(ringdata)
  library("grid")
  
  # Set the color limits, and change all values outside of those limits to the max values
  maxval = 3
  minval = -3
  
  # Set the width of the lines for the ring
  linewidth = 0.5 # 0.1 for thick lines, linewidth for thin lines
  
  # Make half of data on left, half on right
  
  if (dim(ringdata[[1]])[1]>0){
    ringdata[[1]][ ringdata[[1]] > maxval] = maxval
    ringdata[[1]][ ringdata[[1]] < minval] = minval
    ringdata[[1]]$width =   0.5/dim(ringdata[[1]])[1]
    ringdata[[1]]$width = cumsum(ringdata[[1]]$width)
  }else{
    ringdata[[1]] = data.frame(V2 = NaN, width = 0.5)
  }
 
  # Make half of data on left, half on right
  if (dim(ringdata[[2]])[1]>0){
    ringdata[[2]][ ringdata[[2]] > maxval] = maxval
    ringdata[[2]][ ringdata[[2]] < minval] = minval
    ringdata[[2]]$width =   0.5/dim(ringdata[[2]])[1]
    ringdata[[2]]$width = 0.5 + cumsum(ringdata[[2]]$width)
  }else{
    ringdata[[2]] = data.frame(V2 = NaN, width = 1)
  }
  
  plotdata1 = ringdata[[1]]; colnames(plotdata1) = c('V2','width')
  plotdata2 = ringdata[[2]]; colnames(plotdata2) = c('V2','width')
  
  # Make some extra plotdata for the center of the ring
  plotdata3 = data.frame('V2' = mean(plotdata1$V2))
  plotdata4 = data.frame('V2' = mean(plotdata2$V2))
  
  # If the number of genes is larger than 25, cut them down
  if (dim( plotdata1 )[1] > 25){
    olddim = dim(plotdata1)[1]
    plotdata1 = plotdata1[1:25,]
    plotdata1[,2] = plotdata1[,2]*olddim/25
  }
  
  if (length(grep('jpg',path2save))>0){w = 2;h = 2}else{w=5;h=5}
  
  # Set names as well
  plotdata1$Name = rownames(plotdata1)
  plotdata2$Name = rownames(plotdata2)
  
  #pdf(path2save,width = 8,height= 8)
  print( ggplot(plotdata1, aes(fill=V2, ymax=width, ymin=c(0, head(width, n=-1)), xmax=5, xmin=3)) +
           geom_rect(colour="grey30",size = linewidth) + #optionally change geom_rect size = .1 OR linewidth for very thin lines
           #geom_rect(colour="#ffffff",size = 0) +
           coord_polar(theta="y") +
           xlim(c(0, 5)) +
           theme_bw() +
           theme(panel.background = element_rect(fill = NA),panel.border = element_blank()) +
           theme(panel.grid=element_blank()) +
           theme(axis.text=element_blank()) +
           theme(axis.ticks=element_blank()) +
           #labs(title=pathway) + 
           theme(plot.margin = unit(c(0,0,0,0),'cm')) + 
           scale_fill_gradient2(low = '#0000FF', mid = '#FFFFFF', high = '#FF0000', midpoint = 0,guide = FALSE,na.value = "#6B6B6B",limits = c(minval,maxval))  +
           geom_rect(data=plotdata2, size = linewidth, colour="grey30", xmax=5, xmin=3, aes(ymax=width, ymin=c(0.5, head(width, n=-1)))) +
           geom_rect(data = plotdata3, size = linewidth, colour="grey30", xmax=2.5, xmin=0, aes(ymax=0.5, ymin= 0)) + 
           geom_rect(data = plotdata4, size = linewidth, colour="grey30", xmax=2.5, xmin=0, aes(ymax=1, ymin= 0.5)) + 
           #ggtitle(pathway) + 
           geom_vline(xintercept=5,size = linewidth) + 
           geom_vline(xintercept=3,size = linewidth) +
           ggsave( file = path2save,width =w, height = h,dpi = 2000 )
         )
  dev.off()
}