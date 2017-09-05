# A script to calculate association of each metabolite in PanCan metabolomics data with available clinical features like grade and stage.

import os, sys, pandas as pd, numpy as np, scipy as sp, csv, scipy.stats as st, matplotlib.pyplot as plt, matplotlib.cm as cm, itertools as it, pdb, seaborn as sns

sns.set_style('ticks')

sys.path.append('..')
import reportertools as rt

norm = rt.MidpointNormalize(midpoint=0) # Normalize colormaps

plt.close('all')
plt.ion()

#########################################################################################
# Input parameters
#########################################################################################
clintype = str( sys.argv[1] )
#########################################################################################

# Read in the clinical data
clin = pd.io.parsers.read_csv('../../data/merged_metabolon/clinfeatures.csv', header = 0)
clin['PatID'] = np.nan
for item in range(clin.shape[0]):
	clin.ix[item,0] = str(clin.ix[item,0])
	clin.ix[item,'PatID'] = ':'.join( [clin.ix[item,1],clin.ix[item,0],clin.ix[item,2]] )

# Read in the pancan metabolite file
alldata = pd.io.parsers.read_csv('../../data/merged_metabolon/alldata.csv', index_col = 0, header = 0)
# Make sure to drop first two rows
studylist = alldata.loc['Study']
tissuelist = alldata.loc['TissueType']

# Define the studies you are interested in, or leave as 'all' to do all
uqstudy = ['KIRC','BRCA','PRAD']
#uqstudy = 'all'

if uqstudy == 'all':
	# Find the unique kinds of studies
	uqstudy = np.unique( studylist )
	uqstudy = [item for item in uqstudy if item !='nan']
	
	alldata = alldata.ix[2:,:]

else:
	colidx = [item for item in range(len(studylist)) if studylist[item] in uqstudy]
	alldata = alldata[colidx]
	studylist = [studylist[item] for item in colidx]
	tissuelist = [tissuelist[item] for item in colidx]
	
	alldata = alldata.ix[2:,:]

# Merge all the data
alldata = alldata.transpose()
datamerge = pd.merge( right = alldata, left = clin, right_index = True, left_on = 'PatID')
datamerge = datamerge.transpose()
datamerge.columns = datamerge.ix[5,:] # This renames the columns so that we have appropriate IDs

# Create a numpy array to store the results
rescol = [item + ':Tumor' for item in uqstudy] + [item + ':Normal' for item in uqstudy]
res_r = pd.DataFrame(index = datamerge.index, columns = rescol)
res_p = pd.DataFrame(index = datamerge.index, columns = rescol)

# For each study, trim the data and then test differential abundance
studyctr = 0

for study in uqstudy:
	
	# Grab the patient names for this study
	tumidx = np.where( ( datamerge.ix[2,:] == 'Tumor' ) & ( datamerge.ix[1,:] == study ))[0]
	normidx = np.where( ( datamerge.ix[2,:] == 'Normal' ) & ( datamerge.ix[1,:] == study ))[0]
	
	tumdata = datamerge.ix[:,tumidx].dropna(how='all')
	normdata = datamerge.ix[:,normidx].dropna(how='all')
	
	if clintype not in tumdata.index:
		print 'No ' + clintype  + ' in study ' + study
		continue
	
	# Do the association for tumor samples
	for met in np.arange(clin.shape[1],tumdata.shape[0]):
		
		r = st.spearmanr( tumdata.ix[met,:], tumdata.ix[clintype,:] )
		res_r.ix[tumdata.index[met],study+':Tumor'] = r[0]
		res_p.ix[tumdata.index[met],study+':Tumor'] = r[1]#*tumdata.shape[0]
		
	# Do the association for normal samples
	for met in np.arange(clin.shape[1],normdata.shape[0]):
		
		r = st.spearmanr( normdata.ix[met,:], normdata.ix[clintype,:] )
		res_r.ix[normdata.index[met],study+':Normal'] = r[0]
		res_p.ix[normdata.index[met],study+':Normal'] = r[1]#*normdata.shape[0]
	
	studyctr += 1

# Now just clean the results up a little bit
rdrop = res_r.dropna(how='all')
pdrop = res_p.dropna(how='all')
rdropmat = rdrop.as_matrix().astype(float)
rdropmat[np.where(np.isnan(rdropmat))] = 0

# Filter out insignificant correlations
rdropmat[ np.where(pdrop>(0.05)) ] = 0

# Get rid of metabolites with no association
nzmets = np.where( np.sum(rdropmat,axis=1) != 0 )[0]
nzrdropmat = rdropmat[nzmets,:]
nzmetnames = [rdrop.index[item] for item in nzmets]
	
# Plot results
plt.figure( 1 )
plt.imshow( nzrdropmat, norm = norm,interpolation = 'nearest',cmap=cm.seismic, aspect='auto' )
plt.xticks( range(len(res_r.columns)), res_r.columns, fontsize=10 )
plt.yticks( range(len(nzmetnames)), nzmetnames , fontsize = 6 )
plt.colorbar()