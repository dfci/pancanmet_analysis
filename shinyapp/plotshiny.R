# Plot Shiny data
rm(list=ls())
setwd('/Users/ereznik/Documents/reporter/shinyapp/')

library("devtools")
#install_github("ropensci/plotly")
#devtools::install_github("ropensci/plotly")

library(plotly)
#set_credentials_file("ereznik", "6r4hibg89c")


library(ggplot2)

# Load data
load('reporter_shinydata.RData')

#type = 'Fold Change'
type = 'All Data'
metabolite = 'Lactate'

if (type == 'Fold Change'){
    
    fcdata = data.frame(x = fc[,metabolite],y = fc$Study)
    colnames(fcdata) = c('X','Study')
    
    # Find studies which have notNA data
    notNA = unique( fcdata[ which(!is.na(fcdata[,1])),2 ])
    keepidx = which(fcdata[,2] %in% notNA)
    
    plotdata = fcdata[keepidx,]
    plotdata$Study = factor(plotdata$Study)
    
    # Plot
    p2show = ggplot(plotdata,aes(x=X,fill=Study) ) + geom_bar(aes(y = ..density..)) + theme_classic() + 
        geom_density(alpha = 0.5) + ylab(metabolite) + xlab('Log2 Ratio, Tumor:Normal') + 
        geom_vline(xintercept = 0,linetype = 'longdash') + 
        facet_grid(Study~.)
}

if (type == 'All Data'){
    
    alldata = data.frame(x = metdata[,metabolite],y = metdata$Study, z = metdata$Type)
    colnames(alldata) = c('X','Study','Type')
    
    # Find studies which have notNA data
    notNA = unique( alldata[ which(!is.na(alldata[,1])),2 ])
    keepidx = which(alldata[,2] %in% notNA)

    plotdata = alldata[keepidx,]
    plotdata$Study = factor(plotdata$Study)
    
    # Plot
    p2show = ggplot(plotdata,aes(x=X,fill=Type) ) + geom_bar(aes(y = ..density..)) + theme_classic() + 
        geom_density(alpha = 0.5) + ylab(metabolite) + xlab('Log2 Ratio, Tumor:Normal') + 
        facet_grid(Study~.,scales = 'free')
    
}

print(p2show)
# py = plotly(username="ereznik", key="6r4hibg89c")
# r = py$ggplotly(p2show)
# r$response$url
