# A script to explore the partitioning of metabolites in pancan data

import os, sys, pandas as pd, numpy as np, scipy as sp, csv, scipy.stats as st, matplotlib.pyplot as plt, matplotlib.cm as cm, pdb, pickle, time, seaborn as sns
sns.set_style('ticks')
sys.path.append('..')
import reportertools as rt
plt.ion()
plt.close('all')

thresh = float(sys.argv[1]) # Statistical threshold

alldata = pd.io.parsers.read_csv('../../data/merged_metabolon/tnpaired_fc.csv', index_col = 0, header = 0)
tissue = np.asarray([item.split(':')[0] for item in alldata.columns])
uqtissue = np.unique(tissue)

datamat = alldata.as_matrix()
datamat[np.isnan(datamat)] = 0

met = np.where(alldata.index == 'Glycine')[0]
p = np.ones(( datamat.shape[0], len(uqtissue) ))
sign = np.zeros(( datamat.shape[0], len(uqtissue) ))

# Partition
up = np.where(datamat[met,:] > 1)[1]
down = np.where(datamat[met,:] < 1)[1]
nz = np.hstack(( up,down ))

# Make sure there are adequate samples
if len(down) < 5 or len(up) < 5:
	print 'Error not enough up or down samples!'
	
# Now check that we have some spread across tissues
if len(np.unique(tissue[down])) < 2 or len(np.unique(tissue[up])) < 2:
	print 'Error not enough tissues for this metabolite!'

# Test other metabolites
for tissuectr in range(len(uqtissue)):
		
	downtis = [item for item in down if tissue[item] == uqtissue[tissuectr]]
	uptis = [item for item in up if tissue[item] == uqtissue[tissuectr]]
	
	for met2 in range(datamat.shape[0]):
		
		if np.std( datamat[met2,uptis] ) < 1e-10 or np.std( datamat[met2,downtis] ) < 1e-10:
			continue
		
		if len(downtis) == 0 or len(uptis) == 0:
			continue
		
		updata = datamat[met2,uptis]
		downdata = datamat[met2,downtis]
		p[met2, tissuectr] = st.ttest_ind( updata,downdata )[1]
		
		if p[met2, tissuectr]<thresh:

			# Determine whether data is larger in up or down partition
			if np.average( updata ) > np.average(downdata):
				sign[met2,tissuectr] = 1
			else:
				sign[met2,tissuectr] = -1
		
# Find interesting metabolites

rsums = np.sum(np.abs(sign),axis = 1)
partmets = np.where(rsums > 2)[0]

# Make heatmap of results
tmp1,tmp2,index = rt.cluster(sign[partmets,:],0)
plotdata = sign[partmets[index],:]

plt.imshow( plotdata, interpolation = 'nearest',cmap=cm.seismic,aspect='auto' )
plt.colorbar()
plt.yticks( range(len(partmets)), [alldata.index[partmets[item]] for item in index], fontsize=8 )
plt.xticks( range(len(uqtissue)), uqtissue, fontsize=5, rotation  = 90 )
	
	