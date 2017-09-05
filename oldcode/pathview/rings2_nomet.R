# Function to make donut ring plot in ggplot

rings2_nomet = function( ringdata,pathway,path2save ){
  #print(ringdata)
  
  # Make half of data on left, half on right
  ringdata[[1]]$width =   0.5/dim(ringdata[[1]])[1]
  ringdata[[1]]$width = cumsum(ringdata[[1]]$width)
  
  pdf(file = path2save)
  print(ringdata)
  plotdata1 = ringdata[[1]]; colnames(plotdata1) = c('V2','width')
  
  # Make some extra plotdata for the center of the ring
  plotdata3 = data.frame('V2' = mean(plotdata1$V2))
  
  print( ggplot(plotdata1, aes(fill=V2, ymax=width, ymin=c(0, head(width, n=-1)), xmax=5, xmin=3)) +
           geom_rect(colour="grey30") +
           coord_polar(theta="y") +
           xlim(c(0, 5)) +
           theme_bw() +
           theme(panel.background = element_rect(fill = NA)) +
           theme(panel.grid=element_blank()) +
           theme(axis.text=element_blank()) +
           theme(axis.ticks=element_blank()) +
           labs(title=pathway) + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0)  +
           geom_rect(data = plotdata3, colour="grey10", xmax=2.5, xmin=0, aes(ymax=0.5, ymin= 0))
  )
  dev.off()
  
}