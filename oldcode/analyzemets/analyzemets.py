# Script to test for differential abundance of metabolites

import os,sys, scipy as sp, numpy as np, pdb, csv, scipy.stats as st, statsmodels.sandbox.stats.multicomp as padjust

study = str(sys.argv[1])

#########################################################################################
def r2z(a,b,na,nb):
	na = np.float(na)
	nb = np.float(nb)
	
	# Find where a or b equal one
	noneidx = np.where((a!=1) & (b!=1))
	
	r2z = np.zeros(len(a))
	
	r2z[noneidx] = (0.5*np.log((1+a[noneidx])/(1-a[noneidx])) - 0.5*np.log((1+b[noneidx])/(1-b[noneidx])))/(np.sqrt(1/(na-3)+1/(nb-3)))
	
	return r2z
#########################################################################################
# Calculate differential abundance in metabolites

# Load data
d = np.load('../../data/metabolomics/' + study + '/' + study + '_importeddata.npz')
met = d['met']

# Automatically detect the type of sample
f = open('../../data/metabolomics/'+study+'/' + study + '_metdata.csv','rU')
r = csv.reader(f)
rowctr = 0
metnames = []
	
for row in r:
	
	if rowctr == 1:
		tumortype = ['Normal' if 'NORMAL' in item.upper() else 'Tumor' if 'TUMOR' in item.upper() else 'NA' for item in row[1:]]
		normidx = [item for item in range(len(tumortype)) if tumortype[item] == 'Normal'] 
		tumidx = [item for item in range(len(tumortype)) if tumortype[item] == 'Tumor']

	if rowctr > 1:
		metnames.append(row[0])
		
	rowctr += 1
f.close()

# Make sure that we have the same number of metabolites as rows in the data
if d['met'].shape[0] != len(metnames):
	print 'Error!'
	pdb.set_trace()

# For each metabolite, calculate differential abundance
# First column is average fold ratio
# Second column is p-value from Mann-Whitney U Test

diffab = np.ones((met.shape[0],3))

for i in range(met.shape[0]):
	
	# Calculate fold change
	diffab[i,0] = np.log2( np.average( met[i,tumidx] )/ np.average( met[i,normidx] ) )
	
	if np.std( met[i,tumidx] ) > 1e-10 or np.std( met[i,normidx] )> 1e-10: # Changed to OR, Sept 11,2015
		
		# Calculate Mann-Whitney p-value
		diffab[i,1] = 2*st.mannwhitneyu( met[i,tumidx], met[i,normidx] )[1]
	
		# Also do t-test
		diffab[i,2] = 2*st.ttest_ind( met[i,tumidx], met[i,normidx] )[1]

# Adjust p-values
diffab[:,1] = padjust.multipletests( diffab[:,1], method = 'fdr_bh' )[1]
diffab[:,2] = padjust.multipletests( diffab[:,2], method = 'fdr_bh' )[1]

# Save results
np.savez( '../../data/metabolomics/' + study + '/' + study + '_diffabund.npz', diffab = diffab )

# Also, write differential abundance results to file
f2 = open('../../data/diffabund/'+study+'_diffabund_all.csv','w')
w = csv.writer(f2)
w.writerow(['Metabolite Name','Log2 Fold Change','Mann-Whitney P Value FDR Adjusted','T-test P value FDR Adjusted'])
for metidx in range(len(metnames)):
	w.writerow([metnames[metidx]] + list(diffab[metidx,:]) )
f2.close()