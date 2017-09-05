# File of conventions for naming of studies

# Use a common p-value threshold of 0.05
pthresh = 0.05

names2plot = c('BLCA' = 'Bladder','BRCA' = 'Breast\nTerunuma','BRCATang' = 'Breast\nTang',
               'KIRC' = 'Clear Cell\nKidney','LGG' = 'Glioma','OV' = 'Ovarian','PAAD' = 'Pancreatic\nKamphorst',
               'PAADHussain1' = 'Pancreatic\nZhang1','PAADHussain2' = 'Pancreatic\nZhang2',
               'PRAD' = 'Prostate\nSreekumar','PRADLODA' = 'Prostate\nPriolo')

sourcetissue = c('BLCA' = 'Bladder','BRCA' = 'Breast','BRCATang' = 'Breast',
               'KIRC' = 'Kidney','OV' = 'Ovary','PAAD' = 'Pancreas','LGG' = 'Brain',
               'PAADHussain1' = 'Pancreas','PAADHussain2' = 'Pancreas',
               'PRAD' = 'Prostate','PRADLODA' = 'Prostate')

sourcetissue_color = c('Bladder' = 'gold','Breast' = 'red','Kidney' = 'green','Ovary' = 'orange',
                      'Pancreas' = 'blue','Prostate' = 'purple','Brain' = 'black')
metclasses = c('glycerol 3-phosphate (g3p)' = 'Lipid',
               'kynurenine' = 'Amino Acid',
               'taurine' = 'Sulfonic Acid',
               'aspartate' = 'Amino Acid',
               'uracil' = 'Other',
               '10-nonadecenoate (19:1n9)' = 'Lipid',
               'acetylcarnitine' = 'Carnitine',
               'arginine' = 'Amino Acid',
               'asparagine' = 'Amino Acid',
               'butyrylcarnitine' = 'Carnitine',
               'carnitine' = 'Carnitine',
               'c-glycosyltryptophan*' = 'Other',
                 'deoxycarnitine' = 'Carnitine',
               'eicosenoate (20:1n9 or 11)' = 'Lipid',
               'hexanoylcarnitine' = 'Carnitine',
               'malate' =  'TCA Cycle',
               'myristate (14:0)' = 'Lipid',
               'n-acetylaspartate (naa)' = 'Other',
               'n-acetyl-aspartyl-glutamate (naag)'  = 'Other',
               'nicotinamide adenine dinucleotide (nad+)' = 'Other',
               'proline' = 'Amino Acid',
               'propionylcarnitine' = 'Carnitine',
               'urate' = 'Other',
               'uridine' = 'Nucleosides/Nucleotides',
               'alanine' = 'Amino Acid',
               "arachidonate (20:4n6)" = 'Other',
               'fumarate' ='TCA Cycle',
               'glycerol 2-phosphate' = 'Lipid',
               "glycerophosphorylcholine (gpc)" = 'Other',
               'inosine' = 'Nuclosides/Nucleotides',
               'lactate' = 'Glycolysis',
               'laurate (12:0)' = 'Lipid',
               "linoleate (18:2n6)" = 'Lipid',
               'lysine' = 'Amino Acid',
               'methionine' = 'Amino Acid',
               'myo-inositol' = 'Other')