# A script to compare differential abundance patterns in metabolites across cancer types

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
dropna = int(sys.argv[1])

#########################################################################################

# Read in the pancan metabolite file
alldata = pd.io.parsers.read_csv('../../data/merged_metabolomics/alldata.csv', index_col = 0, header = 0)
# Make sure to drop first two rows
studylist = alldata.loc['Study']
tissuelist = alldata.loc['TissueType']

if dropna == 1:
	alldata = alldata.dropna()

# Define the studies you are interested in, or leave as 'all' to do all
#uqstudy = ['KIRC','ccpap']
uqstudy = 'all'

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

# Create a numpy array to store the results
res_t = np.zeros(( alldata.shape[0], len(uqstudy) ))
res_p = np.ones(( alldata.shape[0], len(uqstudy) ))
res_t[:] = np.nan

# For each study, trim the data and then test differential abundance
studyctr = 0
for study in uqstudy:
	idx = [item for item in range(len(studylist)) if studylist[item] == study]

	# Retain only non-NA data
	d = alldata.ix[:,idx]
	
	# Find tumor and normal indices
	temptissue = [tissuelist[item] for item in idx]
	tumoridx = [item for item in range(len(temptissue)) if (temptissue[item].upper() == 'TUMOR' or temptissue[item].upper() == 'MET')]
	normalidx = [item for item in range(len(temptissue)) if temptissue[item].upper() == 'NORMAL']
	
	# For each metabolite, evaluate if we have enough samples, if so, calculate z score and do Mann Whitney U Test
	for met in range(d.shape[0]):
		
		currdata = d.ix[met,:]
		numdata = np.where( currdata.isnull() )[0]
		tumornna = np.intersect1d( tumoridx, np.where( currdata.notnull() )[0] )
		normalnna = np.intersect1d( normalidx, np.where( currdata.notnull() )[0] )
		
		if (len(numdata) > 0.5*len(idx) or len(tumornna)<2 or len(normalnna)<2 ):
		
			res_t[met,studyctr] = 0
			res_p[met,studyctr] = 1
		
		else:
			
			temp1 = [float(item) for item in currdata[tumornna]]
			temp2 = [float(item) for item in currdata[normalnna]]
			
			stdtumor = np.std( temp1 )
			stdnormal = np.std( temp2 )
			
			if (stdtumor > 1e-10 and stdnormal > 1e-10):
				
				res_p[met,studyctr] = -np.log10( st.mannwhitneyu( temp1,temp2 )[1] )
				if res_p[met,studyctr] > -np.log10(0.05):
					res_t[met,studyctr] = st.ttest_ind( temp1,temp2 )[0]
				else:
					res_t[met,studyctr] = 0
				
				#if np.isnan( res_t[met,studyctr] ):
				#	pdb.set_trace()
			
			else: 
			
				res_t[met,studyctr] = 0
				res_p[met,studyctr] = 1
	
	studyctr += 1

# If we are not dropping NA's, at least filter out rows which don't have at least X studies measuring the metabolite
if dropna == 0:
	tsums = np.sum( np.abs( np.sign( res_t ) ), axis = 1 )
	nzt = np.where( tsums >= 5 )[0]
	res_p = res_p[nzt,:]
	res_t = res_t[nzt,:]
	metnames = [alldata.index[item] for item in nzt]
	
	# Write the results to files
	fcfile = pd.DataFrame(res_t,index = metnames, columns = uqstudy)
	pfile = pd.DataFrame(res_p,index = metnames, columns = uqstudy)
	#fcfile.to_csv('../../data/merged_metabolomics/fc.csv')
	#pfile.to_csv('../../data/merged_metabolomics/res_p.csv')
	
else:
	metnames = alldata.index
	
# Plot results
plt.figure( 1 )
plt.imshow( res_p, interpolation = 'nearest',cmap=cm.seismic, aspect='auto' )
plt.xticks( range(len(uqstudy)), uqstudy, fontsize=10 )
plt.yticks( range(len(metnames)), metnames , fontsize = 6 )
plt.colorbar()

# Plot results
order = np.argsort( np.average( np.sign(res_t), axis = 1 ) )
plt.figure( 2 )
plt.imshow( res_t[order], norm = norm, interpolation = 'nearest',cmap=cm.seismic, aspect='auto' )
plt.xticks( range(len(uqstudy)), uqstudy, fontsize=10 )
plt.yticks( range(len(metnames)), [metnames[item] for item in order] , fontsize = 6 )
plt.colorbar()