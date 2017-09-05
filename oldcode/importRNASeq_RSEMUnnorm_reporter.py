# This is a script to import expression data for metabolic genes from RSEM UN-NORMALIZED (!!!) data into a usable format for R. It's purpose is to take a giant file and reduces down into something much smaller that R can process quickly. 

import os, sys, csv, pdb, linecache, pandas as pd, numpy as np, scipy as sp   

# Get only the unique genes
d = np.load('../data/uniquegenes.npz')
entrezIDs = d['uqgenes']
	
# Now, open the raw data files
currdir = '/Users/ereznik/Documents/MuClub/Data/RNASeq_Sep112013/RawUnNormData/'
outdir = '../data/rna_unnorm/'
files_in_dir = [f for f in os.listdir(currdir) if f.endswith('.txt')]
files_in_dir = [files_in_dir[0]] # 3 = COAD, KIRC = 8, BRCA = 1, 17 = PRAD, 15 = OV, 1= BLCA

for currentfile in files_in_dir:
	# Initialize an empty dictionary to store Entrez Gene IDs
	TCGAEntrezDict = {}
	
	# Open the file
	openfile = open(currdir+currentfile, "rU")
	
	# First, we need to read the data in
	reader = csv.reader(openfile, delimiter = '\t')
	
	# Do a quick initial run to determine the number of rows and columns
	rownum = 0
	for row in reader:
		# If you are in the first row, count the columns
		if rownum == 0:
			
			# The number of columns - 1 is important
			numcol = len(row) - 1
			
			# The column headers are unique TCGA identifiers which may be useful later
			patientstemp = row 
		
		elif rownum == 1:
			# Retain the column indices of the raw counts
			rawidx = [i for i in range(len(row)) if row[i] == 'raw_count'  ]
			patients = [ patientstemp[ item ] for item in rawidx ]
			
		elif rownum >1:
			# It's not the first row, save the gene name
			TCGAEntrezDict[ row[0].split('|')[1] ]= rownum
			
		rownum += 1

	# Now, for each line in the file that is relevant, get the required information
	parseddata = np.zeros( (len(entrezIDs),len(rawidx)) )
	finalentrezIDs = np.copy(entrezIDs)
	row2del = []
	for i in range(len(entrezIDs)):
		if entrezIDs[i] in TCGAEntrezDict.keys():
			line = linecache.getline(currdir+currentfile, TCGAEntrezDict[entrezIDs[i]]+1)
			splitline = line.split('\t')
			
			# Check that the gene ID of the line we got matches
			if splitline[0].split('|')[1] == entrezIDs[i]:
				tempdata = [np.float( splitline[item] ) for item in rawidx]
				parseddata[i,:] = np.round(tempdata)
			else:
				print 'Entrez ID does not match!'
				pdb.set_trace()
		else:
			# If we don't have data for this gene, get rid of it
			row2del.append(i)
	
	# Get rid of genes we don't have data for
	finalentrezIDs = sp.delete(finalentrezIDs, row2del, 0)
	parseddata = sp.delete(parseddata,row2del,0)
	
	# Save the data by converting to pandas dataframe for easy import to R
	p = pd.DataFrame( parseddata, index = finalentrezIDs, columns = patients )
	p.to_csv( outdir+currentfile.split('.')[0]+'.csv')
		
		