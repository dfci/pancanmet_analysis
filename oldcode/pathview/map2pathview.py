# Script to find map differential expression/abundance into matrix for pathview

##########################################################################################
# Imports and paramaters
##########################################################################################

import os, sys, pdb, csv, numpy as np, scipy as sp, matplotlib.pyplot as plt, itertools as it, seaborn as sns, scipy.stats as st, pickle, networkx as nx

sys.path.append('..')
import reportertools as rt

h2e = pickle.load( open( '../../data/h2e.p', 'rb' ) )
e2h = {v:k for k, v in h2e.items()}

plt.ion()
plt.close('all')

mthresh = float(sys.argv[1]) # Threshold for calling significance
ethresh = float(sys.argv[2]) # Threshold for calling significance for expression data
study = str(sys.argv[3]) # The name of the study we are using

# Load Differential Abundance
d = np.load( '../../data/metabolomics/' + study + '/' + study + '_diffabund.npz')
ab = d['diffab'] # Values are fold change, Mann-Whitney p value, t-test pval

# Load Differential Expression
isexpr = False # Added this line of code because it may have been deleted! check!!!

if os.path.isdir('../../data/studies/'+study+'/diffexpr/'):
	
	if study+'nvsp.csv' in os.listdir('../../data/studies/'+study+'/diffexpr/'):
		isexpr = True
	
		f = open('../../data/studies/'+study+'/diffexpr/'+study+'nvsp.csv','r')
		r = csv.reader(f)
		numrows = -1 # First row is header, start at -1
		for row in r:
			numrows += 1
		expr = np.zeros((numrows,3)) # We can modify this if we don't know # of genes a priori
		rowctr = -1
		exprgenes = list()
		f.seek(0) # Restart at beginning of file
		for row in r:
			if rowctr > -1:

				# Store FC, adjusted p-value, and abundance
				expr[rowctr,:] = [float(row[1]),float(row[5]),float(row[2])]
				exprgenes.append( str(row[0]) )

			rowctr += 1

		# Change last column of expr data so that min is zero
		expr[:,2] = expr[:,2] - np.min(expr[:,2])

##########################################################################################
# Call significant differential expression and differential abundance
##########################################################################################
if isexpr:
	siggenes = np.where( (expr[:,1] < ethresh) & ( np.abs(expr[:,0])>0 ) )[0]

sigmets = np.where( (ab[:,1] < mthresh ) & ( np.abs(ab[:,0])>0) )[0]

##########################################################################################
# Also, read the metabolite annotation file and make dictionary
##########################################################################################
f = open('../../data/studies/' + study + '/' + study + '_FinalMetIDs.csv','rU')
r = csv.reader( f, delimiter = ',' )
rowctr = -1
keggdict = dict()
for row in r:
	if rowctr == -1:
		rowctr = 0
		continue
	else:
		if '|' in row[2]:
			if row[2].split('|')[0] not in [item[0] for item in keggdict.values()]:
				keggid = row[2].split('|')[0]
				keggdict[ rowctr ] = [keggid, row[0]]
		elif (row[2] != 'NA' and row[2]!='nan'): 
			if row[2] not in [item[0] for item in keggdict.values() ]:
				keggid = row[2]
				keggdict[ rowctr ] = [keggid, row[0]]
	rowctr += 1
f.close()
nummets = rowctr

# Write to files
f = open('../../data/pathview/' + study + '_met.csv', 'w')
w = csv.writer( f )
for item in range(nummets):

	if item in keggdict.keys():
		if item in sigmets:
			w.writerow( [keggdict[item][0],ab[item,0],keggdict[item][1]] )
		else:
			w.writerow( [keggdict[item][0],0,keggdict[item][1]] )
f.close()

if isexpr:
	f = open('../../data/pathview/' + study + '_gene.csv', 'w')
	w = csv.writer( f )
	for item in range(len(exprgenes)):
		genename = exprgenes[item]
		if genename in e2h.keys():
			hugoname = e2h[genename]
		else:
			hugoname = ''
		if item in siggenes:
			w.writerow( [genename, expr[item,0],hugoname] )
		else:
			w.writerow( [genename, 0,hugoname] )
	f.close()