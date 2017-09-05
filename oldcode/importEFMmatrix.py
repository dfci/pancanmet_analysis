# Script to import Metallo S matrix and save as .mat

import os, sys, numpy as np, scipy as sp, csv, pdb, pandas as pd, scipy.io as io

d = pd.read_csv('/Users/ereznik/Documents/reporter/data/EFM/metallomodel/MFA_Metallo_Model_S.csv',index_col = 0, header = 0)
d = d.fillna(0)
d2 = d.as_matrix()

rev = pd.read_csv('/Users/ereznik/Documents/reporter/data/EFM/metallomodel/MFA_Metallo_Model_rev.csv',index_col = 0, header = 0)
rev = rev.as_matrix()

io.savemat( '/Users/ereznik/Documents/reporter/data/EFM/metallomodel/MFA_Metallo_Model_S.mat',{'S':d2, 'metnames' : d.index, 'rnames' : d.columns, 'rev' : rev} )