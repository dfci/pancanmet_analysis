# Script to compare our differential abundance results with those from Goveia met-analysis
# We only considered metabolites measured in at least 10 studies in Goveia

rm(list = ls())
setwd('/Users/ereznik/pancanmet_analysis/analysis//')
library(ggplot2)
library(XLConnect)
library(ggrepel)

# Read in the two results
gvdata = readWorksheetFromFile('../checks/checkdata/GoveiaData.xlsx',sheet = 'Sheet1',rownames = 1)
rezdata = read.csv('../results/diffabundance/DifferentialAbundanceSummary.csv',header = TRUE,row.names = 1)

# Drop any metabolites where we could not unambiguously match the name to our dataset
gvdata = gvdata[-which(gvdata$OurName == 'NA'),]

# Check that our name matching is correct
if (length(which(gvdata$OurName %in% rownames(rezdata))) != dim(gvdata)[1]){stop('Error in name matching!')}

# Change rownames and compare
gvdata$gvnames = rownames(gvdata)
rownames(gvdata) = gvdata$OurName

# Grab data we want to compare
rezdata$DAScore = (rezdata$NumUp - rezdata$NumDown)/rezdata$notNA
gvdata$Reznik_DAScore = rezdata[ rownames(gvdata),'DAScore']
gvdata$Reznik_Size = rezdata[ rownames(gvdata),'notNA']
gvdata$Goveia_DAScore = gvdata$vote.count/gvdata$number.of.studies.that..report.on.the.metabolite
gvdata$LogP = log10(gvdata$multiple.testing.adjusted.pNAvalues..permutation..for.metabolites.reported..at.least.6.times)

ggplot(gvdata,aes(Reznik_DAScore,Goveia_DAScore, label = gvnames)) + 
   geom_point(aes(size = -LogP)) + theme_bw(base_size = 14) + 
   geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
   scale_color_gradient2() + geom_text_repel(color = 'black',size = 2,force = 2) + 
   scale_size_continuous(name = paste0('-Log10 P Value,\nGoveia')) +
   #geom_smooth(method='lm',se = FALSE,formula = y~x) + 
   xlim(-1,1) + ylim(-1,1) + 
   xlab('Reznik, Differential Abundance Score') + ylab('Goveia, Differential Abundance Score') + 
   ggsave('../results/diffabundance/CompareGoveia.pdf',height = 8,width = 8, useDingbats = FALSE)

# Calculate correlation
corval = cor.test(gvdata$Reznik_DAScore,gvdata$Goveia_DAScore,method = 'pearson')

# Make a table too
gvdata$ReznikSign = sign(gvdata$Reznik_DAScore)
gvdata$GoveiaSign = sign(gvdata$Goveia_DAScore)
