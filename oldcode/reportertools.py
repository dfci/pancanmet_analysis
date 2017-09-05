# A library of useful functions for reporter metabolite analysis
import os, sys, scipy as sp, numpy as np, pickle, re, pdb, matplotlib.pyplot as plt, csv, itertools as it

#########################################################################################
# Simple function to cluster
#########################################################################################	
def cluster( rdsig, pltidx ):
	# A simple function that does hierarchical clustering using tools from scipy
	
	import scipy.cluster.hierarchy as sch
	import scipy.spatial.distance as spd
	
	d = spd.pdist(rdsig) 
	L = sch.linkage(d, method='complete')
	
	if pltidx == 1:
		Z = sch.dendrogram(L)
	else:
		Z = sch.dendrogram( L, no_plot = True )
	
	index = Z['leaves']
	
	return d,L,index

#########################################################################################
# Return percentile
#########################################################################################
def pctile( data, pct ):

	import scipy.stats as st
	val = st.scoreatpercentile(data,pct)
	
	return val
	
#########################################################################################
# Normalize colors
#########################################################################################
from matplotlib.colors import Normalize
class MidpointNormalize(Normalize):

	def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
		self.midpoint = midpoint
		Normalize.__init__(self, vmin, vmax, clip)

	def __call__(self, value, clip=None):
		# I'm ignoring masked values and all kinds of edge cases to make a
		# simple example...
		x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
		return np.ma.masked_array(np.interp(value, x, y))



#########################################################################################
# Arrange normal and tumor samples in pairs
#########################################################################################
def pairsamples( studylist, tissuelist, study, alldata, tndict, coldict ):
	import pandas as pd
	
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
	
	if len(pair_norm) == 0:
		normdata = 'Empty'
		tumdata = 'Empty'
		
		return tumdata,normdata
	
	# Construct the paired data
	normdata = alldata.ix[:,pair_norm]
	normdata = normdata.astype('float')
	tumdata = alldata.ix[:,pair_tum]
	tumdata = tumdata.astype('float')
	
	return tumdata, normdata