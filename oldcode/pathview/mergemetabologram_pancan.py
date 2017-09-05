# A script to make metabolograms across studies 

import os, sys, numpy as np, pdb
from PIL import Image

usepdf = False

datadir = '/Users/ereznik/Documents/pancanmet/results/pathway_v2/'
studies = ['BLCA','BRCA','BRCATang','COAD','KIRC','OV','PAAD','PAADHussain1','PAADHussain2','PRAD','PRADLODA','STAD']

if usepdf:
	pways = ['_'.join( item.split('_')[1:3] ) for item in os.listdir(datadir + 'BRCA/metabolograms/') if item.endswith('.pdf')]
else:
	pways = ['_'.join( item.split('_')[1:3] ) for item in os.listdir(datadir + 'BRCA/metabolograms/') if item.endswith('.jpg')]
#studies = ['KIRC']

# For each item, find the corresponding items in the other directories and paste
pathwayctr = 0
pways = [pways[1]] # for testing

# Set a scaling factor
scaling = 3280,3280
for pway in pways:

	ctr = 0
	size = [np.nan]
		
	for study in studies:
	
		if study + '_' + pway in os.listdir( datadir + study + '/metabolograms/' ):
			
			# Convert the image to eps, then read it in
			if usepdf:
				os.system('pdftops -eps {0}'.format(datadir + study + '/metabolograms/' + study + '_' + pway))
			
				s = Image.open(datadir + study + '/metabolograms/' + study + '_' + pway[:-4] + '.eps')
			
				# Remove the .eps
				os.system('rm -f ' + datadir + study + '/metabolograms/' + study + '_' + pway[:-4] + '.eps')
			
			else:
				s = Image.open(datadir + study + '/metabolograms/' + study + '_' + pway[:-4] + '.jpg')
				
		else:
			ctr += 1
			'Skipping study ' + study
			continue

		# Detect the size
		if np.isnan(size[0]):
			size = s.size
			blank_image = Image.new("RGB", (size[0]*4, size[1]*3), 'white' )	

		blank_image.paste(s, ( size[0]*np.mod(ctr,4),size[1]*int(np.floor(ctr/4)) ) )
		
		ctr += 1
		print ctr,pathwayctr

	blank_image.thumbnail(scaling, Image.ANTIALIAS)
	blank_image.save( datadir +'merged/' + pway + '_metabolograms.pdf')


	
	
	
	