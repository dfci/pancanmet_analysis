# A script to correlate normal metabolite patterns to tumor metabolite patterns

# Make sure matplotlib does not print to window
import matplotlib
matplotlib.use('Agg')

import os, sys, pandas as pd, numpy as np, scipy as sp, csv, scipy.stats as st, matplotlib.pyplot as plt, matplotlib.cm as cm, itertools as it, pdb, pickle, time, seaborn as sns, sklearn.linear_model

sns.set_style('ticks')

sys.path.append('..')
import reportertools as rt

from sklearn import cross_validation as cv
from sklearn.linear_model import MultiTaskLassoCV as mtlassocv
from sklearn.linear_model import LassoCV as lassocv
from sklearn.linear_model import Lasso as lasso
from sklearn.linear_model import Ridge as ridge

norm = rt.MidpointNormalize(midpoint=0)

plt.close('all')
plt.ion()

#########################################################################################
# Input parameters
#########################################################################################
thresh = 1e-2 #float(sys.argv[1])
k = int(sys.argv[3])

regressor = 'lasso'
regdict = dict()
regdict[ 'lasso' ] = lasso
regdict[ 'ridge' ] = ridge
#########################################################################################
# Read in Data
#########################################################################################
# Read in tumor/normal pairs
tndict = pickle.load( open( "../../data/tndict.p", "rb" ) )

# Read in the pancan metabolite file
alldata = pd.io.parsers.read_csv('../../data/merged_metabolon/alldata.csv', index_col = 0, header = 0)
# Make sure to drop first two rows
studylist = alldata.loc['Study']
tissuelist = alldata.loc['TissueType']

# Make a data frame for storing the results
beta_pos = pd.DataFrame( np.zeros(( alldata.shape[0],alldata.shape[0] )), index = alldata.index, columns = alldata.index )
beta_neg = pd.DataFrame( np.zeros(( alldata.shape[0],alldata.shape[0] )), index = alldata.index, columns = alldata.index )
fc = []

# Define the studies you are interested in, or leave as 'all' to do all
if sys.argv[2] != 'all':
	uqstudy = sys.argv[2].split()  # Input is study names separated by space
else:
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

# Make a dictionary of the column names
coldict = dict()
for item in range(alldata.shape[1]):
	coldict[ alldata.columns[item] ] = item

#########################################################################################
# Analyze Data
#########################################################################################	
studyctr = 0

for study in uqstudy:

	print study
	
	# Find all tumors
	tumidx = [item for item in range(len(studylist)) if (studylist[item] == study and tissuelist[item] == 'Tumor') ]
	normidx = [item for item in range(len(studylist)) if (studylist[item] == study and tissuelist[item] == 'Normal') ]
	
	# Now arrange tumors so that everything is paired and in order
	pair_tum = []
	pair_norm = []
	for item in tumidx:
		if alldata.columns[item] in tndict.keys():
			pair_tum.append( coldict[ alldata.columns[item] ] )
			pair_norm.append( coldict[ tndict[ alldata.columns[item] ] ] )
	
	# Construct the paired data
	normdata = alldata.ix[:,pair_norm]
	normdata = normdata.astype('float')
	tumdata = alldata.ix[:,pair_tum]
	tumdata = tumdata.astype('float')
	
	if len(pair_norm) == 0:
		continue
	
	# Format data for regression
	x = normdata.dropna()
	y = tumdata.dropna()
	xt = x.as_matrix().transpose()
	
	beta = np.zeros(( x.shape[0], x.shape[0] ))
	alphas = np.logspace( -3, 3, 300 )
	k_fold = cv.KFold(n=x.shape[1], n_folds=k)
	r = np.zeros(x.shape[0])
	p = np.zeros(x.shape[0])
	
	# Create a dictionary indicating how metabolites map
	mapdict = dict()
	for item in x.index:
		mapdict[ item ] = [idx for idx in range(alldata.shape[0]) if alldata.index[idx] == item][0]
	
	for i in range(y.shape[0]):
		# Set alpha array to use initially to be the one above
		alphas_temp = np.copy(alphas)
		
		curmet = x.index[i]
		print 'Iteration for metabolite ' + curmet + ' number ' + str(i)
			
		ytemp = y.as_matrix().transpose()[:,i]
		ytemp = ytemp.reshape((len(ytemp),1))

		if np.std(ytemp) < 0.001:
			continue

		alphaflag = 3 # A flag to indicate if we caught best fit at end of range of alphas, we will only run every metabolite at most 3 times

		while alphaflag != 0:
			# Do CV
			cv_fit = np.zeros((x.shape[1],len(alphas_temp)))
			mse = np.zeros(( k,len(alphas_temp) ))

			for alphactr in range(len(alphas_temp)):
				alpha = alphas_temp[alphactr]
				for train, test in k_fold:
					cv_fit[test,alphactr] = np.squeeze( regdict[regressor](alpha=alpha).fit(xt[train,:], ytemp[train]).predict(xt[test]) )
	

			# Calculate mse
			mse = np.sum( (cv_fit - np.tile(ytemp,(1,len(alphas))))**2, axis = 0 )/np.sum( (ytemp-np.average(ytemp))**2)
	
			if np.argmin(mse) == 0:
				alphas_temp = alphas_temp/1000
				alphaflag = alphaflag -1
				print 'Redoing ' + curmet
			elif np.argmin(mse) == len(mse)-1:
				alphas_temp = alphas_temp*1000
				alphaflag = alphaflag - 1
				print 'Redoing ' + curmet
			else:
				alphaflag = 0
		
		# Calculate minimum alpha
		minalpha = alphas_temp[np.argmin(mse)]
		
		# Calculate correlation
		tmpr = st.pearsonr( cv_fit[:,np.argmin(mse)], np.squeeze(ytemp) )
		r[i] = tmpr[0]
		p[i] = tmpr[1]
		
		# Display results
		if p[i] < thresh and r[i] > 0 and np.min(mse) < 0.8:
		
			# Repeat everything using LOOCV
			print 'Doing LOOCV for metabolite ' + curmet 
			k_fold2 = cv.KFold(n=x.shape[1], n_folds=x.shape[1])
			tempfit = np.zeros( ytemp.shape[0] )
			for train, test in k_fold2:
					tempfit[test] = np.squeeze( regdict[regressor](alpha=alphas_temp[np.argmin(mse)]).fit(xt[train,:], ytemp[train]).predict(xt[test]) )
			
			ploocv = st.pearsonr( tempfit, np.squeeze(ytemp) )[1]
			mseloocv = np.sum( (np.squeeze(ytemp) - tempfit)**2, axis = 0 )/np.sum( (np.squeeze(ytemp)-np.average(np.squeeze(ytemp)))**2)
			
			if ploocv > thresh or mseloocv > 0.8:
				print 'Skipping ' + curmet + ', LOOCV failed.'
				continue # This means that our LOOCV predictions are not good
			
			plt.figure()
			plt.subplot(2,1,1)
			plt.semilogx(alphas_temp, mse, 'k', label='Average across the folds', linewidth=2)
			plt.axvline(alphas_temp[np.argmin(mse)], linestyle='--', color='k', label='alpha: CV estimate')

			plt.legend()

			plt.xlabel('-log(alpha)')
			plt.ylabel('Mean square error')
			plt.title('Mean square error over all folds')
			plt.axis('tight')
			
			plt.subplot(2,1,2)
			plt.plot( tempfit, np.squeeze(ytemp),'o' )
			plt.xlabel('Predicted Value')
			plt.ylabel('Measured Value')
			plt.title( curmet )
		
			plt.savefig('../../results/cornormtumor/' + study + '_' + curmet + '_' + regressor + '.pdf')
			plt.close()
			
			# Save good predictors
			beta[:,i] = regdict[regressor](alpha = minalpha).fit(xt,ytemp).coef_
	
	# Also, store the positive beta's in a big matrix
	nzidx_pos = np.where( beta > 0 )
	if len(nzidx_pos[0]) == 0:
		continue
	xidx_pos = np.array( [mapdict[item] for item in x.index[nzidx_pos[0]]] )
	yidx_pos = np.array( [mapdict[item] for item in x.index[nzidx_pos[1]]] )
	nzidx_mute = tuple(( xidx_pos, yidx_pos ))

	# Add to beta_neg
	beta_pos.ix[nzidx_mute] += 1

	# Also, store the negative beta's in a big matrix
	nzidx_neg = np.where( beta < 0 )
	if len(nzidx_neg[0]) == 0:
		continue
	xidx_neg = np.array( [mapdict[item] for item in x.index[nzidx_neg[0]]] )
	yidx_neg = np.array( [mapdict[item] for item in x.index[nzidx_neg[1]]] )
	nzidx_mute = tuple(( xidx_neg, yidx_neg ))

	# Add to beta_neg
	beta_neg.ix[nzidx_mute] += 1

	# Print the significant interactions
	nzpredictors = np.where( np.sum(beta, axis = 1)!=0 )[0]
	nzy = np.where( np.sum(beta, axis = 0)!= 0)[0]
	trimbeta = beta[nzpredictors,:]
	trimbeta = trimbeta[:,nzy]
	
	nzpredictors_names = [x.index[item] for item in nzpredictors]
	nzy_names = [x.index[item] for item in nzy]
	
	plt.figure()
	plt.imshow(trimbeta, norm = norm, interpolation = 'nearest',cmap=cm.seismic )
	plt.yticks( range(len(nzpredictors)), nzpredictors_names , fontsize=8 )
	plt.xticks( range(len(nzy)), nzy_names, fontsize=8, rotation  = 90 )	
	plt.colorbar()
	plt.savefig('../../results/cornormtumor/betamaps/' + study + '_' + regressor + '.pdf')
	
	np.savez('../../data/cor_norm_tum/beta_' + study + '.npz',trimbeta = trimbeta, nzpredictors_names = nzpredictors_names, nzy_names = nzy_names)
	
	# Also write to CSV
	tempbeta = pd.DataFrame( trimbeta, index = nzpredictors_names, columns = nzy_names )
	tempbeta.to_csv('../../data/cor_norm_tum/beta_' + study + '.csv')
	
	studyctr += 1

if uqstudy != 'all':
	sys.exit()

# After all is said and done, we can plot the recurrent predictors
beta_posmat = beta_pos.dropna().as_matrix()
nzidx = np.where(beta_posmat > 1 )
tmp = beta_posmat[np.unique(nzidx[0]),:]
tmp = tmp[:,np.unique(nzidx[1])]

plt.figure()
plt.imshow(tmp, interpolation = 'nearest',cmap=cm.seismic )
plt.yticks( range(len(np.unique(nzidx[0]))), [alldata.index[item] for item in np.unique(nzidx[0])], fontsize=5 )
plt.xticks( range(len(np.unique(nzidx[1]))), [alldata.index[item] for item in np.unique(nzidx[1])], fontsize=5, rotation  = 90 )
plt.colorbar()
plt.subplots_adjust(bottom=0.2)
plt.savefig('../../results/cornormtumor/betamaps/all_pos_' + regressor + '.pdf')

tempbeta = pd.DataFrame( betaposmat, index = [alldata.index[item] for item in np.unique(nzidx[0])], columns = [alldata.index[item] for item in np.unique(nzidx[1])] )
tempbeta.to_csv('../../data/cor_norm_tum/beta_all_pos.csv')

# Repeat for the negative data
beta_negmat = beta_neg.dropna().as_matrix()
nzidx = np.where(beta_negmat > 1 )
tmp = beta_negmat[np.unique(nzidx[0]),:]
tmp = tmp[:,np.unique(nzidx[1])]

plt.figure()
plt.imshow(tmp, interpolation = 'nearest',cmap=cm.seismic )
plt.yticks( range(len(np.unique(nzidx[0]))), [alldata.index[item] for item in np.unique(nzidx[0])], fontsize=5 )
plt.xticks( range(len(np.unique(nzidx[1]))), [alldata.index[item] for item in np.unique(nzidx[1])], fontsize=5, rotation  = 90 )	
plt.colorbar()
plt.subplots_adjust(bottom=0.2)
plt.savefig('../../results/cornormtumor/betamaps/all_neg_' + regressor + '.pdf')

tempbeta = pd.DataFrame( betanegmat, index = [alldata.index[item] for item in np.unique(nzidx[0])], columns = [alldata.index[item] for item in np.unique(nzidx[1])] )
tempbeta.to_csv('../../data/cor_norm_tum/beta_all_neg.csv')

# Save recurrent predictors
if uqstudy == 'all':
	np.savez('../../data/cor_norm_tum/beta.npz',beta_posmat = beta_posmat,beta_negmat = beta_negmat)