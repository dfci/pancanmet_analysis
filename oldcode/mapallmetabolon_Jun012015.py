# A script to map all metabolon data onto one consensus matrix

import os, sys, csv, numpy as np, scipy as sp, scipy.stats as st, matplotlib.pyplot as plt, pdb, pandas as pd


#####################################################################################
# Functions

def createnewentry( dictlist, current, nextentry, studymap, nextentryflag):
	# We have a totally new metabolite, update accordingly

	for item in range(len(dictlist)):
	
		for currentitem in current[item].split('|'):
			
			dictlist[item][ currentitem ] = nextentry
	
	nextentryflag = 1
	studymap.append(nextentry)
	
	return dictlist, current, nextentry, studymap, nextentryflag
#####################################################################################

#studies = ['BLCA','BRCA','BRCATang','COAD','KICH','KIRC','LGG','OV','PAAD','PAADHussein1','PAADHussein2','PRAD', 'PRADLODA','STAD']

#studies = ['BLCA','BRCA','BRCATang','COAD','KIRC','LGG','OV','PAAD','PAADHussain1','PAADHussain2','PRAD', 'PRADLODA','STAD'] # No KICH
studies = ['BLCA','BRCA'] # No KICH

# Initialize dictionaries
pubchemdict = dict()
keggdict = dict()
hmdbdict = dict()
chebidict = dict()
namedict = dict()
dictlist = [namedict, pubchemdict, keggdict, hmdbdict, chebidict]

# Initialize the counter for which row in the data we are adding next
nextentry = 0

# Initialize storage for the final indices to map to
finalmap = []

##########################################################################################
# Part 1: Go through each study and try to merge metabolites
##########################################################################################

# For each study, open the final metabolon id's file, assign a unique row to each metabolite
studyctr = 0
errorctr = 0

for study in studies:
	
	# Open the file
	f = open( '../data/studies/' + study + '/' + study + '_FinalMetIDs.csv', 'rU' )
	r = csv.reader( f, delimiter = ',' )
	rowctr = 0
	
	# The variable study map indicates to which row of the final dataset we will write the metabolites in the current dataset
	studymap = []
	
	# For each metabolite
	for row in r:
		
		# Set a flag indicating whether we have added a new metabolite
		nextentryflag = 0
		
		if rowctr == 0:
			rowctr += 1
			continue # Skip the header
		
		# Get the relevant data
		current = [row[0].title(),row[1],row[2],row[3],row[4]]
		key2add = [ [],[],[],[],[]]
		
		# Initialize an array storing the row each metabolite maps to
		map = []
		
		# For each dictionary
		for item in range(len(dictlist)):
			
			for currentitem in current[item].split('|'):
			
				if currentitem in dictlist[item].keys() and currentitem != '' and currentitem != 'NA' and currentitem != 'nan':
				
					# We have already mapped it, make sure it maps with all other keys
					map.append( dictlist[item][currentitem] )
		
				else:
					
					# This item is not in the dictionary, but we could potentially add it
					key2add[item].append( currentitem )
					
		# Check if we have consensus using the mappings
		if row[0] == 'allo-threonine' or row[0] == 'threonine':
			pdb.set_trace()
		if len(map) == 0:
			
			# We have a totally new metabolite, update accordingly
# 			for item in range(len(dictlist)):
# 			
# 				for currentitem in current[item].split('|'):
# 					
# 					dictlist[item][ currentitem ] = nextentry
# 			
# 			nextentryflag = 1
# 			studymap.append(nextentry)
		
			# Create a new entry
			dictlist, current, nextentry, studymap, nextentryflag = createnewentry( dictlist, current, nextentry, studymap, nextentryflag)
		
		else:
			# Check for consensus and that this metabolite hasn't already been mapped in this study
			if len( np.unique(map) ) == 1 and map[0] not in studymap:
				
				# We have perfect consensus, update the studymap
				studymap.append(map[0])
				
				# Add all keys in key2add
				for item in range(len(key2add)):
					for key in key2add[item]:
						dictlist[item][key] = map[0]
			
			else:
				invname = {v:k for k,v in namedict.items()}
				
				print 'Error, no consensus in mapping for metabolite ' + current[0]+ ' in ' + study + '!'
				
				nameoptions = [invname[item] for item in map]
				print nameoptions

				# Use name if possible
				if current[0] in namedict.keys():
					studymap.append( namedict[current[0]] )
					print 'Using name to match ' + current[0] + '\n'
				else:
					print 'No consensus and not in namedict, requesting user input for ' + current[0] + '!\n'
					
					usr_input = raw_input('Please enter index (first item is index 0) of metabolite name to use for ' + current[0] + ',given the following options, or type NEW to generate a new metabolite: \n' + '\n'.join( nameoptions ) + '\n' )
					
					if usr_input != 'NEW':
						studymap.append( map[int(usr_input)] )
					
					######################################################################
					# This is new as of May 29, 2015
					else:
						
						# We have a totally new metabolite, update accordingly
						for item in range(len(dictlist)):
							for currentitem in current[item].split('|'):
								dictlist[item][ currentitem ] = nextentry
						nextentryflag = 1
						studymap.append(nextentry)
					errorctr += 1
					######################################################################
		
		# Update nextentry
		nextentry = nextentry + nextentryflag
		
		# Update rowctr
		rowctr += 1

	# Update study counter
	studyctr += 1
		
	# Update finalmap
	finalmap.append( studymap )
	
	# Close the file
	f.close()

	#print(keggdict)
	#pdb.set_trace()
##########################################################################################
# Part 2: Return to each study and actually merge using finalmap
##########################################################################################
# First, check that all names match appropriately

# Initialize data array
alldata = np.zeros([ nextentry, len(studies) ], dtype = 'S100' )
studyctr = 0
mettypedict = dict()
submettypedict = dict()

for study in studies:
	
	# Open the file
	f = open( '../data/studies/' + study + '/' + study + '_FinalMetIDs.csv', 'rU' )
	r = csv.reader( f, delimiter = ',' )
	rowctr = 0
	
	# For each metabolite
	for row in r:
		if rowctr == 0:
			rowctr += 1
			continue
		alldata[ finalmap[studyctr][rowctr-1], studyctr ] = row[0]
		mettypedict[ row[0].title() ] = row[1]
		submettypedict[ row[0].title() ] = row[2]
	
		rowctr += 1
	
	studyctr += 1
	f.close()

# Print out
f = open('../data/merged_metabolomics/mergednames.csv', 'w')
w = csv.writer( f, delimiter = ',' )
w.writerow( studies )
for row in range(alldata.shape[0]):
	w.writerow( alldata[row,:] )
f.close()

# Merge into metabolite names
metnames = []
for row in range(alldata.shape[0]):
	notblank = [item for item in np.unique(alldata[row,:]) if item != '']
	if len(np.unique(notblank)) == 1:
		metnames.append( notblank[0].title() )
	else:
		print 'Not a unique metabolite name, using ' + notblank[0].title()
		metnames.append( notblank[0].title() )

# Assuming everything looks good, merge all the data
studyctr = 0
studyrow= []
tumorrow = []
patnames = []

for study in studies:
	
	# Open the file
	f = open( '../data/metabolomics/' + study + '/' + study + '_metdata.csv', 'rU' )
	r = csv.reader( f, delimiter = ',' )
	rowctr = 0
	
	# For each metabolite
	for row in r:
	
		if rowctr == 0:
			# Store patient names
			patnames_study = [':'.join( [study,row[item]] ) for item in np.arange(1,len(row)) ]

		if rowctr == 1:
			tumortype = ['Normal' if 'NORMAL' in item.upper() else item.title() for item in row[1:]]

			if studyctr == 0:
				
				# Initialize data array
				alldata = np.zeros( ( nextentry, len(tumortype) ) )
			
			else:
				tempdata = np.zeros( ( nextentry, len(tumortype) ) )
			
		if rowctr > 1:
			
			if studyctr == 0:
				alldata[ finalmap[studyctr][rowctr - 2], : ] = row[1:]
			
			else:
				tempdata[ finalmap[studyctr][rowctr - 2], : ] = row[1:]
				
		rowctr += 1
		
	# Stack data if necessary
	if studyctr != 0:
		alldata = np.hstack(( alldata, tempdata ))
	
	# Add data to the study row
	studyrow = studyrow + [study] * (len(row) - 1)
	
	# Add data to the tumor row
	tumorrow = tumorrow + tumortype
	
	# Merge patnames with type of tissue and 
	patnames = patnames + [':'.join([patnames_study[item],tumortype[item]]) for item in range(len(patnames_study))]
	
	studyctr += 1

# Convert into pandas dataframe
alldf = pd.DataFrame( alldata  )
studydf = pd.DataFrame( studyrow )
studydf = studydf.transpose()
tumordf = pd.DataFrame( tumorrow )
studydf = studydf.append( tumordf.transpose() )
studydf = studydf.append( alldf )
alldf = studydf
alldf.index = ['Study','TissueType'] + metnames

# Add a column corresponding to metabolite type and metabolite subtype
mettypename = ['',''] + [mettypedict[item] for item in metnames]
submettypename = ['',''] + [submettypedict[item] for item in metnames]

alldf[ alldf == 0 ] = np.nan
alldf[ alldf == '' ] = np.nan
alldf[ alldf == 'nan' ] = np.nan
print alldf.shape

# Finally, if there are any duplicate names, get merge the rows
dupidx = np.setdiff1d(np.arange(len(alldf.index)), np.unique(alldf.index, return_index=True)[1]) 
dup = [alldf.index[item] for item in dupidx ]
dropidx = []
for name in dup:
	idx = [item for item in range(len(alldf.index)) if alldf.index[item] == name]
	print alldf.index[idx[0]],alldf.index[idx[1]]
	
	# The last two columns are metabolite type names, so be careful
	newdata = [alldf.ix[idx[0],item] if ~np.isnan(alldf.ix[idx[0],item]) else alldf.ix[idx[1],item] for item in range( alldf.shape[1]-2 ) ]
	newdata = newdata + [alldf.ix[idx[0],alldf.shape[1]-2], alldf.ix[idx[0],alldf.shape[1]-1] ]
	alldf.ix[idx[0],:] = newdata
	dropidx.append( idx[1] )

# Drop duplicate rows
idx2keep = np.setdiff1d( range(alldf.shape[0]),dupidx )
alldf = alldf.ix[idx2keep,:]

# Add column names
alldf.columns = patnames

# Save to csv
alldf.to_csv( '../data/merged_metabolomics/alldata.csv', na_rep = 'NA')
		
##########################################################################################
# Part 3: Write out a list of IDs for each metabolite
##########################################################################################
invdictlist = []
for item in dictlist:
	invdictlist.append( {k:v for v,k in item.items()} )
	
# Write a file with the ids
f = open('../data/merged_metabolomics/merged_IDs.csv', 'w')
w = csv.writer( f, delimiter = ',' )
for i in range(alldf.shape[0]):
	w.writerow( [hotdict.get(i,'NA') for hotdict in invdictlist] )
f.close()