# Script to munge the KEGG pathway data Arman handed to us

import os, sys, csv, pdb, pickle

# Read in the h2e dictionary and invert it
h2e = pickle.load( open( "../data/h2e.p", "rb" ) )
e2h = {v:k for k,v in h2e.items()}

f = open('/Users/ereznik/Documents/pancanmet/data/pathway_v2/kegg_details.tsv','rU')
r = csv.reader( f, delimiter = '\t' )
fout = open('/Users/ereznik/Documents/pancanmet/data/pathway_v2/kegg_details_munged.csv','w')
w = csv.writer( fout )
rowctr = 0
for row in r:
	if rowctr == 0:
		rowctr = 1
		continue
	
	pathway = row[0].split(':')[1]
	temp = row[1].replace(',','|')
	tempgenes = temp.split('|')
	tempgenes2 = [item.split(':')[1] for item in tempgenes]
	hugogenes = [e2h[item] for item in tempgenes2 if item in e2h.keys()]
	genes = '|'.join( hugogenes )
	compounds = row[4]
	
	w.writerow( [pathway,genes,compounds] )
	
	print pathway

fout.close()
f.close()
	
	