# A script to import and normalize the prostate cancer metabolomics data from Sreekumar. Our strategy will be to impute in the same way as the remainder of the Metabolon data, which is to use the lowest recorded value as the imputed value for the remainder of the samples, then median normalize each row.

import os, sys, numpy as np, scipy as sp, csv, pdb

fin = open('../../data/metabolomics/LGG/LGG_metabolon_raw_unnormalized.csv','rU')
fout = open('../../data/metabolomics/LGG/LGG_metdata.csv','w')

r = csv.reader( fin )
w = csv.writer( fout )
rowctr = 0
numskipped = 0
for row in r:
	
	if rowctr == 0:
		w.writerow( ['Biochemical Name'] + range(len(row)-1) )
		rowctr = 1
		
		w.writerow([''] + ['Tumor']*(len(row)-1))
	
	else:
		
		# Note that we need to get rid of commas
		dataidx = [item for item in np.arange(1,len(row)) if row[item] != '']
		naidx = [item for item in np.arange(1,len(row)) if row[item] == '']
		med = np.median( [np.float( row[item].replace(',', '') ) for item in dataidx] )
		
		# Impute data
		imp = np.min( [float(row[item].replace(',', '')) for item in dataidx] )
		normdata = [float(row[item].replace(',', ''))/med if row[item] != '' else imp/med for item in np.arange(1, len(row))] 
		
		#normdata = [item/med for item in newdata]
		
		finaldata = [row[0]] + normdata
		
		w.writerow( finaldata )
		
fin.close()
fout.close()