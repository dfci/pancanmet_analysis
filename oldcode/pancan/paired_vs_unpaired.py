# A script to analyze the difference between using paired tests to using unpaired tests in the pancanmet project.

import os, sys, numpy as np, scipy.stats as st, pandas as pd, pdb, matplotlib.pyplot as plt, seaborn as sns, statsmodels.sandbox.stats.multicomp as padjust

sns.set_style('ticks')

# Set p-value threshold
pthresh = 0.05

# Read in the two files
unpair = pd.read_csv('../../data/merged_metabolomics/alldata.csv',header = 0,index_col = 0, skiprows = [1,2] )

pair = pd.read_csv('../../data/merged_metabolomics/tnpaired_fc.csv',header = 0,index_col = 0)

# Get the unique studies
uqstudy = np.unique([item.split(':')[0] for item in pair.columns])

# For each study
for study in uqstudy:

	# Get the indices corresponding to these in both files
	pairix = [item for item in pair.columns if item.split(':')[0] == study]
	unpairix = [item for item in unpair.columns if item.split(':')[0] == study]
	
	if len(pairix) == 0:	
		print('No pairs for ' + study + '...skipping.')
		continue
	
	# Subset the data
	subpair = pair.ix[:,pairix].dropna()
	subunpair = unpair.ix[:,unpairix].dropna()
	
	# Get the tumor and normal indices for unpaired data
	tumidx = [item for item in subunpair.columns if item.split(':')[2] == 'Tumor']
	normidx = [item for item in subunpair.columns if item.split(':')[2] == 'Normal']
	
	if list(subpair.index) != list(subunpair.index):
		print('Indexes do not match for paired and unpaired data!')
		pdb.set_trace()
	
	# Store data
	res = pd.DataFrame(columns = ['PairFC', 'PairP', 'PairPAdj', 'UnpairFC', 'UnpairP', 'UnpairPAdj'] )
	
	# For each metabolite, do some testing
	for met in subpair.index:
	
		# Check to make sure there is some variation in the data
		pairstd = np.std(subpair.ix[met,:])
		unpairstd = np.std(subunpair.ix[met,:])
		
		if pairstd < 1e-10 or unpairstd < 1e-10:
			continue
		
		# Calculate FC
		res.ix[met,'PairFC'] = np.log2( np.mean(subpair.ix[met,:]) )
		res.ix[met,'UnpairFC'] = np.log2( np.mean(subunpair.ix[met,tumidx]) / np.mean(subunpair.ix[met,normidx]) )
		
		# Calculate p-values
		res.ix[met,'PairP'] = st.wilcoxon( np.log2( subpair.ix[met,:] ) )[1] # This tests symmetry about zero, so take a log since right now 1 means equal concentration
		res.ix[met,'UnpairP'] = 2*st.mannwhitneyu( subunpair.ix[met,tumidx],subunpair.ix[met,normidx] )[1] # Reports one-sided p, need to double
	
	# Adjust p-values
	res.ix[:,'PairPAdj'] = padjust.multipletests( res.ix[:,'PairP'], method = 'fdr_bh' )[1]
	res.ix[:,'UnpairPAdj'] = padjust.multipletests( res.ix[:,'UnpairP'], method = 'fdr_bh' )[1]
	
	# Make the plot data
	zeropair = np.where(res['PairPAdj'] > pthresh )[0]
	zerounpair = np.where(res['UnpairPAdj'] > pthresh )[0]
	
	res.ix[zeropair,'PairFC'] = 0
	res.ix[zerounpair,'UnpairFC'] = 0
	
	plt.figure(1)
	plt.plot( res['PairFC'],res['UnpairFC'],'o')
	plt.xlabel('Paired Data')
	plt.ylabel('Unpaired Data')
	plt.title(study)
	
	x1,x2,y1,y2 = plt.axis()
	minplot = np.min([x1,y1])
	maxplot = np.max([x2,y2])
	plt.plot([minplot, maxplot], [minplot,maxplot], 'k-', lw=2)
	
	plt.savefig('../../results/paired_vs_unpaired/' + study+ '.pdf',dpi = 300)
	plt.close()
	