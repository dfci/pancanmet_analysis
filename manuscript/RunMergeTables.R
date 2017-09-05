options(java.parameters = "-Xmx8192m")

library(xlsx)

# Read tables 
siTable1 <- read.csv("tables/merged_metabolomics.csv")
siTable2 <- read.csv("tables/tumor_normal_pairs.csv")
siTable3 <- read.csv("tables/clinfeatures.csv")
siTable4 <- read.csv("tables/DifferentialAbundanceSummary.csv")
siTable6 <- read.csv("tables/hmdbClassification.csv")
siTable7 <- read.csv("tables/AllPathway_diffAbundanceBig.csv")
siTable8 <- read.csv("tables/kegg_drugbank_results.csv")
siTable9 <- read.csv("tables/merged_IDs.csv")
siTable10 <- read.csv("tables/pairscores.csv")
siTable11 <- read.csv("tables/Covariation_Results.csv")

file.copy("tables/updated results stage and grade Jan 2017 max1.xlsx", "finalfigures_tables/SI_Table2.xlsx", overwrite=TRUE)

# Write tables
write.xlsx2(siTable1, "finalfigures_tables/SI_Table1.xlsx", sheetName="merged_metabolomics")
write.xlsx2(siTable2, "finalfigures_tables/SI_Table1.xlsx", sheetName="tumor_normal_pairs", append=TRUE)
write.xlsx2(siTable3, "finalfigures_tables/SI_Table1.xlsx", sheetName="clinfeatures", append=TRUE)
write.xlsx2(siTable4, "finalfigures_tables/SI_Table1.xlsx", sheetName="DifferentialAbundanceSummary", append=TRUE)
write.xlsx2(siTable6, "finalfigures_tables/SI_Table1.xlsx", sheetName="hmdbClassification", append=TRUE)
write.xlsx2(siTable7, "finalfigures_tables/SI_Table1.xlsx", sheetName="AllPathway_diffAbundanceBig", append=TRUE)
write.xlsx2(siTable8, "finalfigures_tables/SI_Table1.xlsx", sheetName="kegg_drugbank_results", append=TRUE)
write.xlsx2(siTable9, "finalfigures_tables/SI_Table1.xlsx", sheetName="merged_IDs", append=TRUE)
write.xlsx2(siTable10, "finalfigures_tables/SI_Table1.xlsx", sheetName="pairscores", append=TRUE)
write.xlsx2(siTable11, "finalfigures_tables/SI_Table1.xlsx", sheetName="Covariation_Results", append=TRUE)


