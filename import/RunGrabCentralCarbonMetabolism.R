# A script to grab the relevant glycolysis metabolites from the data file

rm(list = ls())
setwd('/Users/ereznik/pancanmet_analysis/import/')
library(data.table)

metnames = c('glucose','glucose-6-phosphate (g6p)','fructose-6-phosphate',
             'fructose-1,6-bisphosphate','dihydroxyacetone phosphate (dhap)',
            'dihydroxyacetonephosphate (dhap)','glyceraldehyde-3-phosphate',
            '3-phosphoglycerate','2-phosphoglycerate',
            'phosphoenolpyruvate (pep)','lactate','citrate','cis-aconitate',
            'succinate','fumarate','malate','alpha-ketoglutarate','glutamate','glutamine')

# Read in the differential abundance results and write the resulting table
res = read.csv('../results/diffabundance/DifferentialAbundanceSummary.csv',header = TRUE,row.names = 1)
res2 = res[metnames,]
write.csv(res2,'../results/diffabundance/DiffAbundance_CentralCarbonMetabolism.csv')