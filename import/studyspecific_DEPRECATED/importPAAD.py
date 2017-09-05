# A script to import and the prostate cancer metabolomics data

import os, sys, numpy as np, scipy as sp, csv, pdb

fin = open('../../data/studies/PAAD/tempPAAD_metabolomics.csv','rU')
fout = open('../../data/metabolomics/PAAD/PAAD_metdata.csv','w')

r = csv.reader( fin )
w = csv.writer( fout )
rowctr = 0
numskipped = 0
skiprows = []
for row in r:
	
	if rowctr == 0:
		w.writerow( row )
		patids = row
	
		ttype = ['Tumor' if 'Tumor' in item else 'Normal' if 'Normal' in item else '' for item in patids]
		w.writerow( ttype )

	else:
		
		dataidx = [item for item in np.arange(1,len(row)) if 'NA' not in row[item]]
		naidx = [item for item in np.arange(1,len(row)) if 'NA' in row[item]]		
			
		# Impute data and re-expeonentiate
		imp = np.min( [float(row[item]) for item in dataidx] )
		newdata = [2**float(row[item]) if 'NA' not in row[item] else 2**imp for item in np.arange(1, len(row))]
		
		med = np.median( [float(newdata[item-1]) for item in dataidx] ) # Note that the indexing is off by 1
		
		normdata = [item/(med) for item in newdata]
	
		finaldata = [row[0]] + normdata
	
		w.writerow( finaldata )
		
	rowctr += 1
		
fin.close()
fout.close()

