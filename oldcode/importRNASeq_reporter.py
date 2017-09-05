# A script to import TCGA RNA-Seq data for reporter metabolite analysis.

import os, sys, csv, pdb, numpy as np, linecache, pandas as pd

study = 'KIRC'

# Get only the unique genes
d = np.load('../data/uniquegenes.npz')
genes = d['uqgenes']

# Go to the TCGA RNASeq normalized file for KIRC
currfile = '/Users/ereznik/Documents/MuClub/Data/RNASeq_Sep112013/RawData/' + study + '.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt'
f = open(currfile,'r')

r = csv.reader(f, delimiter = '\t')

# Quickly scan the file, get the gene names
rowctr = 0
genedict = dict()
for row in r:
	
	if rowctr == 0:
		pats = row[1:]
	
	if rowctr > 1:
		genedict[ row[0].split('|')[1] ] = rowctr
	
	rowctr += 1

# Now for each gene in genes, go to the appropriate line and read in the data
intgenes = [item for item in genes if item in genedict.keys()]
parseddata = np.zeros((len(intgenes),len(pats)))
itemctr = 0
for item in intgenes:
	
	line = linecache.getline(currfile, genedict[item] + 1)
	splitline = line.split('\t')
	
	# Check that the gene ID of the line we got matches
	if splitline[0].split('|')[1] == item:
	
		parseddata[itemctr,:] = [np.float( item ) for item in splitline[1:]]
	
	else:
		print 'Entrez ID does not match!'
		pdb.set_trace()
	
	itemctr += 1

p = pd.DataFrame( parseddata, index = intgenes, columns = pats )
p.to_csv( '../data/rna/' + study + '_RNA_Norm.csv')
