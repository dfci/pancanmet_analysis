# A script to examine the co-variation of one metabolite in particular in the pancan metabolomics data.

import os, sys, pandas as pd, numpy as np, scipy as sp, csv, scipy.stats as st, matplotlib.pyplot as plt, matplotlib.cm as cm, itertools as it, pdb

sys.path.append('..')
import reportertools as rt

plt.close('all')
plt.ion()

#########################################################################################
# Input parameters
#########################################################################################
hotmet = '2-Hydroxyglutarate'

# Read in the pancan metabolite file
alldata = pd.io.parsers.read_csv('../../data/merged_metabolon/alldata.csv', index_col = 0, header = 0)
# Make sure to drop first two rows
studylist = alldata.loc['Study']
tissuelist = alldata.loc['TissueType']

# Define the studies you are interested in, or leave as 'all' to do all
if len(sys.argv) > 2:
	uqstudy = sys.argv[2].split()  # Input is study names separated by space
	colidx = [item for item in range(len(studylist)) if studylist[item] in uqstudy]
	alldata = alldata[colidx]
	studylist = [studylist[item] for item in colidx]
	tissuelist = [tissuelist[item] for item in colidx]
	
	alldata = alldata.ix[2:]
else:
	uqstudy = 'all'
	
	# Find the unique kinds of studies
	uqstudy = np.unique( studylist )
	uqstudy = [item for item in uqstudy if item !='nan']
	colidx = [item for item in range(len(studylist)) if studylist[item] in uqstudy]
	studylist = [studylist[item] for item in colidx[:-2]]
	tissuelist = [tissuelist[item] for item in colidx[:-2]]
	
	alldata = alldata.ix[2:]


# For each type of unique study, calculate the co-variation among metabolites
studyctr = 0
plotctr = 1
numsamples = np.zeros(len(uqstudy))

# Make one big dataframe
tempcols = [item + '_Tumor' for item in uqstudy] + [item + '_Normal' for item in uqstudy]
res = pd.DataFrame(index = alldata.index, columns = tempcols)

for study in uqstudy:

	print study
	
	# Find all tumors
	tumidx = [item for item in range(len(studylist)) if (studylist[item] == study and tissuelist[item] == 'Tumor') ]
	normidx = [item for item in range(len(studylist)) if (studylist[item] == study and tissuelist[item] == 'Normal') ]
	
	# Retain only non-NA data
	d = alldata.ix[:,tumidx + normidx]
	d = d.dropna()
	
	# Make sure this study has measured the metabolite in question
	if hotmet not in d.index:
		continue
	else:
		hotidx = np.where(d.index == hotmet)[0][0]
	
	tumdata = d.ix[:,0:len(tumidx)].as_matrix()
	normdata = d.ix[:,len(tumidx):].as_matrix()
	
	# Create a dictionary indicating how metabolites map
	mapdict = dict()
	for item in d.index:
		mapdict[ item ] = [idx for idx in range(alldata.shape[0]) if alldata.index[idx] == item][0]
	
	# Calculate covariation of metabolites in this study
	corrcalc = st.spearmanr( np.transpose( tumdata ) )
	corrcalc_normal = st.spearmanr( np.transpose( normdata ) )
	
	# Get rid of nans
	corrcalc[0][ np.isnan(corrcalc[0]) ] = 0
	corrcalc[1][ np.isnan(corrcalc[0]) ] = 1
	
	corrcalc_normal[0][ np.isnan(corrcalc_normal[0]) ] = 0
	corrcalc_normal[1][ np.isnan(corrcalc_normal[0]) ] = 1
	
	# Filter out insignificant correlations
	corrcalc[0][np.where(corrcalc[1]>0.05)] = 0
	corrcalc_normal[0][np.where(corrcalc_normal[1]>0.05)] = 0
	
	# Extract data of interest
	hotcorr = corrcalc[0][hotidx,:]
	hotcorr_normal = corrcalc_normal[0][hotidx,:]
	
	# Save it in the big data file
	nz = np.where(hotcorr!=0)[0]
	temp =[d.index[item] for item in nz]
	restumidx = [mapdict[item] for item in temp]
	res.ix[restumidx,study+'_Tumor'] = hotcorr[nz]
	
	nz_normal = np.where(hotcorr_normal!=0)[0]
	temp =[d.index[item] for item in nz_normal]
	resnormidx = [mapdict[item] for item in temp]
	res.ix[resnormidx,study+'_Normal'] = hotcorr_normal[nz_normal]
	
	print [d.index[item] for item in np.argsort(hotcorr)]
	print '\n\n'
	print [d.index[item] for item in np.argsort(hotcorr_normal)]
	print '\n\n'
	
	mergedata = np.vstack((hotcorr,hotcorr_normal)).transpose()
	#tempdf = pd.DataFrame(mergedata,columns = ['Tumor','Normal'],index = d.index)
	#tempdf.to_csv('/Users/ereznik/Documents/randomprojects/inklehofer/2hg_pancan/' + study+'.csv')
	
# Make a heatmap
resdrop = res.dropna(how = 'all')
resdrop = resdrop.fillna(0)
rsums = np.sum( np.sign( np.abs( resdrop ) ), axis =1 )
nzmets = np.where(rsums>4)[0]
plotdata = resdrop.ix[nzmets,:].as_matrix()
temp1,temp2,index = rt.cluster(plotdata,0)
plotdata = plotdata[index,:]
nzmets = [nzmets[item] for item in index]
plt.imshow(plotdata,interpolation = 'nearest',cmap=cm.seismic, aspect = 'auto' )
plt.xticks( range(len(tempcols)), tempcols, fontsize=6, rotation = 90 )
plt.yticks( range(len(nzmets)), resdrop.index[nzmets], fontsize = 6)
plt.colorbar()