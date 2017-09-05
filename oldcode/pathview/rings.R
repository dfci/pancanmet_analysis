# Function to make donut ring plot in ggplot

rings = function( ringdata,pathway,path2save ){
  print(ringdata)
  for (y in 1:2){
    ringdata[[y]]$width =   1/dim(ringdata[[y]])[1]
    ringdata[[y]]$width = cumsum(ringdata[[y]]$width)
  }
  
  pdf(file = path2save)
  print(ringdata)
  plotdata1 = ringdata[[1]]
  plotdata2 = ringdata[[2]]
  print( ggplot(plotdata1, aes(fill=V2, ymax=width, ymin=c(0, head(width, n=-1)), xmax=5, xmin=3)) +
    geom_rect(colour="grey30") +
    coord_polar(theta="y") +
    xlim(c(0, 10)) +
    theme_bw() +
    theme(panel.background = element_rect(fill = NA)) +
    theme(panel.grid=element_blank()) +
    theme(axis.text=element_blank()) +
    theme(axis.ticks=element_blank()) +
    labs(title=pathway) + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0)  +
    geom_rect(data=plotdata2, colour="grey30", xmax=7.5, xmin=5.5, aes(ymax=width, ymin=c(0, head(width, n=-1)))) )
  dev.off()

}