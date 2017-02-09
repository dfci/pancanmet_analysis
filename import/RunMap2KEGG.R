#' Top-level script to map pathway data to KEGG
#' useKEGGREST: if TRUE, new KEGGREST run used to map compound IDs to pathways, otherwise; if FALSE, old data in ../data/KEGGPathways.RData used

# Last Modified: Ed Reznik, March 2, 2016
rm(list = ls())
analysisDir = normalizePath(".")
shinyDir = '../pancanmet/'
#setwd(file.path(analysisDir, "import"))
setwd('/Users/ereznik/pancanmet_analysis/import/')
writedir = '../data/'

# Load necessary scripts and libraries
source('getKEGG.R')
source('map2KEGG.R')
source('checkKEGGIDs.R')

# Parameters for users to set
useKEGGREST = FALSE

# Create a log file for storing run information
logfile = paste('log/PathwayImport_',Sys.time(),'.txt',sep = '')
write('****************************************',logfile,append = TRUE)
write('****************************************',logfile,append = TRUE)
write('Importing pathway data',logfile,append = TRUE)
write('****************************************',logfile,append = TRUE)

# 1. Decide whether to use KEGGREST
if (useKEGGREST){
  write('Starting new KEGGREST run...',logfile,append = TRUE)
  getKEGG()
}else{
  write('Loading old KEGGREST data...',logfile,append = TRUE)
  load('../data/KEGGpathays.RData')
}
write('****************************************',logfile,append = TRUE)

# 2. Use KEGG dictionary to map to pathways
write('Mapping merged metabolomics row names to KEGG pathways...',logfile,append = TRUE)
map2KEGG( logfile )
write('****************************************',logfile,append = TRUE)

# 3. Check that the mapping was reasonable
write('Checking the KEGG dictionary to make sure mappings are reasonable...',logfile,append = TRUE)
checkKEGGIDs( logfile )
write('****************************************',logfile,append = TRUE)

# 4. Copy files to final directory
write('Copying pathway data to ../data/ directory ...',logfile,append = TRUE)
write('****************************************',logfile,append = TRUE)
file.copy(from = 'tempdir/KEGG_pathwaymap.csv',to = paste(writedir,'KEGG_pathwaymap.csv',sep=''),overwrite = TRUE)

write('Done.',logfile,append = TRUE)

# Cleanup: Reset path to project directory
setwd(analysisDir)