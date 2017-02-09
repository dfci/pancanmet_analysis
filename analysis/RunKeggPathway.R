# Script to do simple KEGG pathway analysis

rm(list = ls())
setwd('/Users/ereznik/pancanmet_analysis/analysis/')
library(ggplot2)
library(reshape)

# User parameters
minmet = 5 # The minimum number of metabolites in a pathway required for it to be plotted in the summary figure for the study
minstudies = 6 # The minimum number of studies with at least minmet metabolites
allpath = FALSE # If true, use all KEGG pathways, if false, go up to 83, which seems to include all "detailed" metabolic pathways

# Read in the pathway map
map = read.csv('../data/KEGG_pathwaymap.csv',header = TRUE,row.names = 1,check.names = FALSE)
map = as.matrix(map)
if (!allpath){
  map = map[1:83,]
}

# Read in the metabolite differential abundance file
metchange = read.csv('../results/diffabundance/DifferentialAbundance_SignedChanges.csv',header = TRUE,row.names = 1)

# Restrict analysis only to those metabolites in the pathway map file
ixmets = intersect(rownames(metchange),colnames(map))
if (length(ixmets)!=dim(map)[2]){
  stop('Some metabolites in KEGG_pathwaymap.csv do not map to the differential abundance file!')
}
metchange = metchange[colnames(map),]

# Make two separate matrices, indicating which metabolites go up or down in each study
sigup = as.matrix( metchange )
sigup[ which(sigup <= 0, arr.ind = TRUE) ] = 0
sigup[ which(is.na(sigup),arr.ind = TRUE) ] = 0

sigdown = as.matrix( -metchange ) # Notice that we take a negative of everything here
sigdown[which(sigdown <= 0,arr.ind = TRUE)] = 0
sigdown[ which(is.na(sigdown),arr.ind = TRUE) ] = 0

# Calculate an up and down score for each pathway across each study
upscore = map %*% sigup
downscore = map %*% sigdown

# Calculate the number of measured metabolites
allmeasured = as.matrix( metchange )
allmeasured[which(!is.na(allmeasured),arr.ind = TRUE)] = 1
allmeasured[which(is.na(allmeasured),arr.ind = TRUE)] = 0
nummeasured = map %*% allmeasured

# Calculate differential abundance score for each pathway
diffabundance = (upscore - downscore)/nummeasured

# Save results
write.csv(diffabundance,'../results/pathway/AllPathway.csv')

# For each study, make a summary plot, including any pathway with at least minnum metabolites mapping to it
for (study in colnames(metchange)){
  pdata = data.frame('DA' = diffabundance[,study], 'Size' = nummeasured[,study])
  pdata = pdata[ -which(pdata$Size < minmet), ]
  pdata$Name = as.factor(rownames(pdata))
  
  ggplot(pdata,aes(DA,Name,color = DA)) + geom_point(aes(size = Size)) + 
    theme_bw() + scale_colour_gradient2(low = 'blue', mid = 'grey', high ='red',  midpoint = 0) +
    theme(panel.border = element_blank(), panel.grid.minor = element_blank(),  
          axis.line = element_line(colour = "black")) +
    geom_vline(xintercept = 0,linetype = "dotted") + 
    geom_segment(aes(x=0,y=Name,xend = DA,yend = Name)) + xlab('Differential Abundance Score') + 
    scale_x_continuous(limits = c(-1.1, 1.1)) + ylim(rev(levels(pdata$Name))) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle(study) +
    ggsave(paste('../results/pathway/studysummary/metpathway_',study,'.pdf',sep=''),width = 9, height = 9, useDingbats=FALSE)
  
}

# Get all pathways which measure at least minmet metabolites in at least minstudy studies
totalcount = map %*% allmeasured
totalcount[ which(totalcount < minmet,arr.ind = TRUE) ] = 0
totalcount[ which(totalcount >= minmet,arr.ind = TRUE) ] = 1
totalcount_studies = rowSums(totalcount)
path2keep = which(totalcount_studies > minstudies)
pathnames2keep = names(path2keep)

# Make a big plot
bigdata = diffabundance[path2keep,]
bigsize = nummeasured[path2keep,]

# Rank the pathways
pathrank = rowSums(bigdata,na.rm = TRUE)

bigdata_melt = melt(bigdata)
rownames(bigdata_melt) = paste(bigdata_melt[,1],bigdata_melt[,2],sep = ':')
bigsize_melt = melt(bigsize)
rownames(bigsize_melt) = paste(bigsize_melt[,1],bigsize_melt[,2],sep = ':')

# Bind the data
mergebig = bigdata_melt
colnames(mergebig) = c('Pathway','Study','DA')
mergebig$Size = bigsize_melt[rownames(mergebig),3]
mergebig$Pathway = factor(mergebig$Pathway,levels = rev(names(pathrank)[order(pathrank)]) )

# Plot
ggplot(mergebig, aes(Study,Pathway,color = DA, size = Size)) + geom_point(shape=19) +
  theme_bw() + geom_point(data = subset(mergebig,DA == 0),aes(size = Size),color = 'grey50',shape = 1)+ 
  geom_point(data = subset(mergebig,Size == 0),shape =4,color = 'grey50',size = 5) +
  theme(panel.border = element_blank(), 
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_size_area(name = 'Pathway\nSize') + xlab('') + ylab('') + 
  scale_colour_gradient2(low = 'blue', high ='red',  midpoint = 0,name = 'Differential\nAbundance\nScore') +
  ylim(rev(levels(mergebig$Pathway))) + theme(legend.position="bottom",legend.key = element_blank()) + 
  ggsave('../results/pathway/AllPathway.pdf',height = 10,width = 10,useDingbats = FALSE)

