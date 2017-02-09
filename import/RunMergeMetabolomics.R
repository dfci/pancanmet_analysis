# Top-level script to import metabolomics data

# Last Modified: Ed Reznik, March 2, 2016
rm(list = ls())
analysisDir = normalizePath(".")
shinyDir = '../pancanmet/'
setwd('/Users/ereznik/pancanmet_analysis/import/')
#setwd(file.path(analysisDir, "import"))

# Input parameters
seedstudy = 'KIRC'
writedir = '../data/merged_metabolomics/'

# Import functions
source('getstudies.R')
source('mergemetabolomicsIDs.R')
source('mergemetabolomicsdata.R')
source('mergepaireddata.R')
source('runmergeddiagnostics.R')
source('useful_metimport.R')
source('append2map.R')

# Create a log file for storing run information
logfile = paste('log/MetabolomicsImport_',Sys.time(),'.txt',sep = '')
write('****************************************',logfile,append = TRUE)
write('****************************************',logfile,append = TRUE)
write('Pan-Cancer Metabolomics Merge',logfile,append = TRUE)
write('****************************************',logfile,append = TRUE)

# 0. Make sure that tempdir is created and empty
dir.create(file.path('tempdir'), showWarnings = FALSE)
fs2del = list( list.files('tempdir',full.names = TRUE) )
if (length(fs2del)!=0){
  do.call( file.remove,fs2del )
}

# 1. Detect all study directories
IDdir = '../data/studies'
metdir = '../data/metabolomics'
write('****************************************',logfile,append = TRUE)
write('Checking for metabolomics data and metabolite ID files...\n',logfile,append = TRUE)
studies = getstudies(IDdir,metdir,logfile)

# 1.5 If you like, re-order the studies so that the "most-trusted" ID file is first, this acts as a seed for the rest of the merger
if (!is.na(seedstudy)){
  # Modify order of studies so that seedstudy is first
  studies = c(seedstudy,setdiff(studies,seedstudy))
}

# 2. Check that data is median normalized and imputation done correctly (Deprecated July 2016)
#write('****************************************',logfile,append = TRUE)
#write('Beginning to check data is appropriately normalized...\n',logfile,append = TRUE)

# 3. Match metabolite IDs
write('****************************************',logfile,append = TRUE)
write('Beginning to merge metabolomics ID files...\n',logfile,append = TRUE)
mergemetabolomicsIDs( studies,logfile )

# 4. Map actual metabolomics data
write('****************************************',logfile,append = TRUE)
write('Beginning to merge metabolomics data files...\n',logfile,append = TRUE)
mergemetabolomicsdata(logfile)

# 5. Generate metabolomics file for paired tumor/normal samples.
write('****************************************',logfile,append = TRUE)
write('Beginning to merge paired tumor/normal data...\n',logfile,append = TRUE)
path2pairs = '../data/merged_metabolomics/tumor_normal_pairs.csv'
mergepaireddata(path2pairs,logfile)

# 6. Run diagnostics
write('****************************************',logfile,append = TRUE)
write('Beginning to run diagnostics on merged files...\n',logfile,append = TRUE)
runmergeddiagnostics( studies,logfile)

write('****************************************',logfile,append = TRUE)
write('Done merging...\n',logfile,append = TRUE)

# 7. Copy files to final directory
file.copy(from = 'tempdir/merged_metabolomics.csv',to = paste(writedir,'merged_metabolomics.csv',sep=''),overwrite = TRUE)
file.copy(from = 'tempdir/tnpaired_fc.csv',to = paste(writedir,'tnpaired_fc.csv',sep=''),overwrite = TRUE)

# Cleanup: Reset path to project directory
setwd(analysisDir)
