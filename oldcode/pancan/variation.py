# A script to study patterns in variation of metabolites in pancan metabolomics data

import os, sys, pandas as pd, numpy as np, scipy as sp, csv, scipy.stats as st, matplotlib.pyplot as plt, matplotlib.cm as cm, itertools as it, pdb, pickle, time, seaborn as sns

sns.set_style('ticks')

sys.path.append('..')
import reportertools as rt

from matplotlib.backends.backend_pdf import PdfPages

plt.ion()

#########################################################################################
# Input parameters
#########################################################################################
n = int(sys.argv[2]) # Number of permutations to do
numplots = 9
#########################################################################################
# Read in Data
#########################################################################################
# Read in tumor/normal pairs
tndict = pickle.load( open( "../../data/tndict.p", "rb" ) )

# Read in the pancan metabolite file
alldata = pd.io.parsers.read_csv('../../data/merged_metabolomics/alldata.csv', index_col = 0, header = 0)
# Make sure to drop first two rows
studylist = alldata.loc['Study']
tissuelist = alldata.loc['TissueType']

# Define the studies you are interested in, or leave as 'all' to do all
if sys.argv[1] != 'all':
	uqstudy = sys.argv[1].split()  # Input is study names separated by space
	
	colidx = [item for item in range(len(studylist)) if studylist[item] in uqstudy]
	alldata = alldata[colidx]
	studylist = [studylist[item] for item in colidx]
	tissuelist = [tissuelist[item] for item in colidx]
	
	alldata = alldata.ix[2:,:]
	
else:
	uqstudy = 'all'
	# Find the unique kinds of studies
	uqstudy = np.unique( studylist )
	uqstudy = [item for item in uqstudy if item !='nan']
	
	alldata = alldata.ix[2:,:]

# Make a dictionary of the column names
coldict = dict()
for item in range(alldata.shape[1]):
	coldict[ alldata.columns[item] ] = item
	
# Make a pandas array to store results
savedata = pd.DataFrame( np.zeros((len(alldata.index),len(uqstudy))), index = alldata.index, columns = uqstudy )

#########################################################################################
# Analyze Data
#########################################################################################	
studyctr = 0

for study in uqstudy:
	print study
	
	plt.close('all')
	plots = [None]*20
	pp = PdfPages('../../results/variation/' + study + '.pdf')
	
	tumdata, normdata = rt.pairsamples( studylist, tissuelist, study, alldata, tndict, coldict )
	
	if 'Empty' in normdata:
		continue
	
	# Format data
	normdata = normdata.dropna()
	tumdata = tumdata.dropna()
	
	normmat = normdata.as_matrix()
	tummat = tumdata.as_matrix()
	
	# Calculate coefficient of variation for each metabolite
	tumcv = np.std( tummat, axis = 1 ) / np.average( tummat, axis = 1 )
	normcv = np.std( normmat, axis = 1) / np.average( normmat, axis = 1 )
	diffcv = np.asarray([tumcv[item] - normcv[item] if (tumcv[item]!=0 and normcv[item]!=0) else 0.0 for item in range(len(tumcv)) ] )
	
	# Compare this CV to experimentally determined CV
	catdata = np.hstack((tumdata,normdata))
	numtum = tummat.shape[1]
	numnorm = normmat.shape[1]
	shufidx = range( numtum+numnorm )
	
	# Create storage array
	permcv = np.zeros(( diffcv.shape[0], n ))
	p_high = np.zeros(( diffcv.shape[0] ))
	p_low = np.zeros(( diffcv.shape[0] ))
	for perm in range(n):
		print perm
		
		# Shuffle
		np.random.shuffle( shufidx )
		temptum = catdata[:,shufidx[0:numtum]]
		tempnorm = catdata[:,shufidx[numtum:]]
		
		# Test
		temptumcv = np.std( temptum, axis = 1 ) / np.average( temptum, axis = 1 )
		tempnormcv = np.std( tempnorm, axis = 1 ) / np.average( tempnorm, axis = 1 )
		permcv[:,perm] = np.asarray( [temptumcv[item] - tempnormcv[item] if (temptumcv[item]!=0 and tempnormcv[item]!=0) else 0.0 for item in range(len(temptumcv)) ] )
		
	# Calculate p-values
	for met in range(diffcv.shape[0]):
		p_high[met] = float(len( np.where( permcv[met,:]>diffcv[met] )[0] ))/float(n)
		p_low[met] = float(len( np.where( permcv[met,:]<diffcv[met] )[0] ))/float(n)
	
	p_high[ np.where(diffcv< 2)[0] ] = 1
	p_low[ np.where(diffcv > 0.5 )[0] ] = 1
	
	zmets = np.where( (np.std( tummat, axis = 1 ) < 0.1) | (np.std( normmat, axis = 1 ) < 0.1) )[0]
	
	p_high[zmets] = 1
	p_low[zmets] = 1

	# Find significant mets
	sigmets_high = np.where(p_high*diffcv.shape[0]<0.05)[0]
	sigmets_low = np.where(p_low*diffcv.shape[0]<0.05)[0]
	
	# Plot results
	plt.plot( tumcv, normcv, 'o')
	plt.xlabel('Tumor CV')
	plt.ylabel('Normal CV')
	
	anndata = zip( [tumdata.index[item] for item in sigmets_high], tumcv[sigmets_high],normcv[sigmets_high] )
	for label, x, y in anndata:
		plt.annotate( label, xy = (x,y) )
		plt.plot(x,y,'go')
		
	anndata = zip( [tumdata.index[item] for item in sigmets_low], tumcv[sigmets_low],normcv[sigmets_low] )
	for label, x, y in anndata:
		plt.annotate( label, xy = (x,y) )
		plt.plot(x,y,'ro')
	
	plt.savefig('../../results/variation/'+study+'_allpts.pdf')
	plt.close()
	
	# For each sigmet, plot as a boxplot plot
	sigmets = list(sigmets_high) + list(sigmets_low)
	p = [p_high[item] for item in sigmets_high] + [p_low[item] for item in sigmets_low]
	
	if len(sigmets) == 0:
		continue
		
	for plotctr in range(len(sigmets)):
		
		plots[int(np.ceil(plotctr/numplots))] =  plt.figure(1 + np.ceil(plotctr/numplots),figsize = (10,10))
		plt.subplot(np.sqrt(numplots),np.sqrt(numplots),np.mod(plotctr+1,numplots))			
		normmean = np.average( normmat[ sigmets[plotctr],: ] )
		tummean = np.average( tummat[ sigmets[plotctr],: ] )
		
		plotdata = [normmat[ sigmets[plotctr],: ]/normmean,tummat[ sigmets[plotctr],: ]/tummean ]
		sns.violinplot( plotdata, color = 'pastel',names=["Normal", "Tumor"] )
		plt.title(normdata.index[sigmets[plotctr]])
		sns.despine(trim=True)
		plt.xlabel( ' p = ' + str(round(p[plotctr],4)) )
		plt.subplots_adjust(hspace=0.6)
		
	for i in range(int(np.ceil(float(plotctr+1)/float(numplots)))):
		pp.savefig(1+i)
	pp.close()
	
	for met in sigmets_high: 
		savedata.ix[ normdata.index[met],study ] = 1
	for met in sigmets_low: 
		savedata.ix[ normdata.index[met],study ] = -1
		
# Now make a figure showing the results
plotdata = savedata.dropna().as_matrix()
nzmets = np.where( np.sum( np.abs( plotdata ), axis = 1 ) != 0 )[0]
plotdata = plotdata[nzmets,:]

# Cluster plotdata
tmp1,tmp2,index = rt.cluster( plotdata, 0 )

plt.figure()
plt.imshow( plotdata[index,:], interpolation = 'nearest',cmap=cm.seismic )
plt.yticks( range(len(nzmets)), [alldata.index[nzmets[item]] for item in index], fontsize=5 )
plt.xticks( range(len(uqstudy)), uqstudy, fontsize=5, rotation  = 90 )
plt.colorbar()