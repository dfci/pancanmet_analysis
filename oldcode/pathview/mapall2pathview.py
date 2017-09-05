# Script to map all data to pathview

import os, sys, pdb

studies = ['BLCA','BRCA','BRCATang','COAD','KICH','KIRC','LGG','OV','PAAD','PAADHussain1','PAADHussain2','PRAD', 'PRADLODA','STAD']

curdir = '/Users/ereznik/Documents/pancanmet/analysis/pathview/'

for study in studies:

	print(study)
	callstr = 'cd ' + curdir + '; python /Users/ereznik/Documents/pancanmet/analysis/pathview/map2pathview.py 5e-2 5e-2 ' + study
	print callstr
	
	os.system(callstr)
