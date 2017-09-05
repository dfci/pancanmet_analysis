# A script to calculate differential abundance of all metabolites

import os, sys, numpy as np, pdb

studies = ['BLCA','BRCA','BRCATang','COAD','KICH','KIRC','LGG','OV','PAAD','PAADHussain1','PAADHussain2','PRAD','PRADLODA','STAD']

for s in studies:

	print(s)
	
	os.system('python analyzemets.py ' + s)