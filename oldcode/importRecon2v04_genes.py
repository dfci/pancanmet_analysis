# A script to import Recon2v.04 and find genes which use certain metabolites as substrates/products

import os, sys, numpy as np, scipy as sp, pandas as pd, pdb, pickle

from scipy.io import loadmat

# Import the Recon2 structure
recon = loadmat('../data/Recon2.v04.mat',struct_as_record=True)['modelR204'][0,0]

h2e = pickle.load( open( '/Users/ereznik/Documents/pancanmet/data/h2e.p', 'rb' ) )
e2h = {v:k for k, v in h2e.items()}

genenames = [str( recon['genes'][item][0][0].split('.')[0] ) for item in range(len(recon['genes'])) ]

uqgenenames = np.unique(genenames)

# Convert these names to hugo if possible
hugo = [e2h.setdefault(item,'NA') for item in uqgenenames]

g2write = pd.DataFrame( uqgenenames, columns = ['Entrez'] )
g2write['GeneName'] = hugo
g2write.to_csv('../data/Recon2v04_genenames.csv')