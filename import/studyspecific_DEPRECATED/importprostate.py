# A script to import and normalize the prostate cancer metabolomics data from Sreekumar. Our strategy will be to impute in the same way as the remainder of the Metabolon data, which is to use the lowest recorded value as the imputed value for the remainder of the samples, then median normalize each row.

# NOTE June 1, 2015: Noticed that 2-Hydroxybutyrate (AHB) is profiled by both GC and LC. Decided to keep LC and toss GC. Note that this needs to be accounted for in the metabolite IDs file as well.

import os, sys, numpy as np, scipy as sp, csv, pdb

fin = open('../../data/studies/PRAD/prostate_raw_unnormalized.csv','rU')
fout = open('../../data/metabolomics/PRAD/PRAD_metdata.csv','w')

ahbcounter = 0

r = csv.reader( fin )
w = csv.writer( fout )
rowctr = 0
numskipped = 0
for row in r:
	
	if rowctr == 0:
		
		# Also write tumor type
		ttype = ['Tumor' if 'T' in row[item] else 'Normal' if 'N' in row[item] else 'Metastasis' if 'M' in row[item] else 'NA' for item in range(len(row))]
		
		col2write = [item for item in range(len(ttype)) if ttype[item] != 'NA']
		
		head2write = [row[0]] + [row[item] for item in col2write]
		ttype2write = [ttype[item] for item in col2write]
		
		w.writerow( head2write )
		w.writerow([''] + ttype2write)
		rowctr += 1
		
	else:
	
		if ahbcounter == 0 and row[0] == '2-Hydroxybutyrate (AHB)':
			ahbcounter = 1
			continue #i.e., skip this entry
		
		dataidx = [item for item in col2write if row[item] != '']
		naidx = [item for item in col2write if row[item] == '']
		
		# Impute data
		med = np.median( [float(row[item]) for item in dataidx] )
		imp = np.min( [float(row[item]) for item in dataidx] )
		newdata = [float(row[item]) if row[item] != '' else imp for item in col2write] 
		
		# Because the data is on the wrong scale, let's normalize
		normdata = [item/med for item in newdata]
		
		finaldata = [row[0]] + normdata
		w.writerow( finaldata )
		
fin.close()
fout.close()