# Visualize clinical analysis results

rm(list = ls())
setwd('/Users/ereznik/pancanmet_analysis/analysis/')
source('plottingconventions.R')
library(XLConnect)
library(ggplot2)
library(cowplot)
library(reshape2)
library(ggrepel)
pthresh = 5e-2
lowestp = 1e-5

# Read in clinical data
clin = read.csv('../data/merged_metabolomics/clinfeatures.csv',header = TRUE,row.names = 1,
                stringsAsFactors = FALSE)

# Read in clinical results
wb = loadWorkbook('../results/clinical/updated results stage and grade Jan 2017 max1.xlsx')
setMissingValue(wb, value = NA)
grade = readWorksheet(wb,header = TRUE,rownames = 1,sheet = 'Grade',startRow = 2)
gradenames = rownames(grade)
stage = readWorksheet(wb,header = TRUE,rownames = 1,sheet = 'Stage',startRow = 2)
stagenames = rownames(stage)

# Read in pathways, restricting to metabolic pathways with at least n (n=5, n=10?) members
pways = read.csv('../data/KEGG_pathwaymap.csv',header = TRUE,row.names = 1,check.names = FALSE)
pways = pways[1:84,]
pways = pways[which(rowSums(pways) >= 10),]

# Convert to numeric
grade = sapply(grade[,1:dim(grade)[2]], function(x){as.numeric(as.character(x))} )
rownames(grade) = gradenames

stage = sapply(stage[,1:dim(stage)[2]], function(x){as.numeric(as.character(x))} )
rownames(stage) = stagenames

for (dname in list('grade','stage')){
#for (dname in list('stage')){  
  d = get(dname)
  
  # Visualize this data
  s = d[,1:(dim(d)[2]-6)] # omit the last 6 columns, which report meta-p-values
  s[is.na(s)] = 1
  insig = which(s > pthresh, arr.ind = TRUE)
  s[insig] = 0
  s[s!=0] = 1
  
  # Separate tumor and normal columns
  tumcols = grep('tumor',colnames(s),ignore.case = TRUE)
  normcols = grep('normal',colnames(s),ignore.case = TRUE)
  
  # Separate out any meta-analysis columns
  metacols = grep('combined',colnames(s),ignore.case = TRUE)
  tumcols = setdiff(tumcols,metacols)
  normcols = setdiff(normcols,metacols)
  
  # Make a big dataframe to store results
  res = data.frame(d[,c('all.combined.adjusted','tumors.combined.adjusted','normals.combined.adjusted')])
  
  # Calculate metabolite sums
  msums = rowSums(s)
  msums_tumor = rowSums(s[,tumcols])
  msums_normal = rowSums(s[,normcols])
  signames = apply(s,1,function(x){y = which(x == 1);paste(colnames(s)[y],collapse = ';')})
  res[,'Significant'] = msums[rownames(res)]
  res[,'Significant_Tumor'] = msums_tumor[rownames(res)]
  res[,'Significant_normal'] = msums_normal[rownames(res)]
  res[,'SigNames'] = signames[rownames(res)]
  
  res = res[order(res$Significant,decreasing = TRUE),]
  
  # Assign a minimal non-zero p-value
  res[is.na(res$all.combined.adjusted),'all.combined.adjusted'] = 1
  res[is.na(res$tumors.combined.adjusted),'tumors.combined.adjusted'] = 1
  res[is.na(res$normals.combined.adjusted),'normals.combined.adjusted'] = 1
  res[res$all.combined.adjusted < lowestp,'all.combined.adjusted'] = lowestp
  res[res$tumors.combined.adjusted < lowestp,'tumors.combined.adjusted'] = lowestp
  res[res$normals.combined.adjusted < lowestp,'normals.combined.adjusted'] = lowestp
  
  # Calculate pathway scores, keeping only those metabolites which are associated in 2 or more tumors (not normals)
  sigmets = rownames(res)[which(res$tumors.combined.adjusted < pthresh)]
  ixmets = intersect(sigmets,colnames(pways)) # pathway file only has a subset of metabolites
  pwaysize = rowSums(pways)
  scores = rowSums( pways[,ixmets] ) / rowSums(pways)
  scores = data.frame(scores)
  scores$names = rownames(scores)
  scores$size = pwaysize[rownames(scores)]
  
  # Figure 1: P-values versus total number of sigmets
  fig1 = ggplot(res,aes(-log10(all.combined.adjusted),Significant)) + geom_point() + 
    xlab('-Log10(Meta-P)') + ylab('# of Significant Associations') +
    ggsave('../results/clinical/pvalues.pdf', height = 10, width = 10)
  
  # Figure 2: Pathway scores
  scores = scores[order(scores$scores),]
  scores$names = factor(scores$names,levels = scores$names)
  fig2 = ggplot(scores,aes(scores,names,size = size)) + geom_point() +
    ggsave('../results/clinical/pathways.pdf',height = 10,width = 10)
  
  # Figure 3: Associations across studies
  uqstudy = unique( sapply(colnames(s),function(x){strsplit(x,'\\.')[[1]][1]}) )
  res_study = data.frame('Tumor' = numeric(0),'Normal' = numeric(0))
  for (study in uqstudy){
    
    tumname = paste(study,'\\.Tumor',sep='')
    normname = paste(study,'\\.Normal',sep = '')
    res_study[study,'Tumor'] = length( grep(tumname,res[sigmets,'SigNames'],ignore.case = TRUE) )
    res_study[study,'Normal'] = length( grep(normname,res[sigmets,'SigNames'],ignore.case = TRUE) )
    #res_study[study,'All'] = length( grep(paste(study,'\\.',sep = ''),res$SigNames,ignore.case = TRUE) )
    res_study[study,'Study'] = study
  }
  res_study_melt = melt(res_study,id.vars = 'Study')
  
  # Add to res_study_melt the number of data points we have
  res_study_melt$SampleSize = 0
  for (rname in rownames(res_study_melt)){
    res_study_melt[rname,'SampleSize'] = length(which(clin$study == as.character(res_study_melt[rname,'Study']) & 
                                                        clin$type == as.character(res_study_melt[rname,'variable'])) )
  }
  res_study_melt$Name = names2plot[res_study_melt$Study]
  res_study_melt$Tissue = sourcetissue[res_study_melt$Study]
  res_study_melt = res_study_melt[-which(res_study_melt$SampleSize == 0),] # remove lines with 0 samples, i.e. LGG normal
  
  fig3 = ggplot(res_study_melt,aes(Study,value,fill = variable)) + 
    geom_bar(stat = 'identity',position = 'dodge') + theme_classic(base_size = 20) + 
    ylab('Number of Metabolites') + xlab(NULL) +
    theme(legend.title=element_blank()) +
    ggsave('../results/clinical/histogram.pdf',height = 6,width = 6)
  
  samplesize_cor = cor.test(res_study_melt$SampleSize,res_study_melt$value,method = 'spearman')
  fig4 = ggplot(res_study_melt,aes(SampleSize,value,color = Tissue,label = Name,shape = variable)) + 
    geom_point(size = 5,alpha = 0.7) + theme_classic(base_size = 10) + geom_text_repel(show.legend = FALSE) +
    xlab('Sample Size') + ylab('Number of Significantly Associated Metabolites') + 
    ggtitle(paste0(dname,' Spearman rho ',round(samplesize_cor$estimate,3),', p-value ',round(samplesize_cor$p.value,3))) + 
    scale_shape_discrete(name = 'Sample Type') + 
    scale_color_manual(values = sourcetissue_color)
    ggsave(paste0('../results/clinical/SampleSize_',dname,'.pdf'),height = 6,width= 6)
  if (dname == 'stage'){figstage = fig4}else{figgrade = fig4}
  
  # Draw full plot
#   rightcolumn = plot_grid(fig1,fig3,labels = c('B','C'),ncol = 1)
#   plot_grid(fig2,rightcolumn,ncol = 2,scale = c(0.5,1))
  print(fig1)
  print(fig2)
  print(fig3)
  print(fig4)

  
  write.csv(res,paste('../results/clinical/Visualized_',dname,'.csv',sep = ''))

  # Run through the results in res, and make a big dataframe of tumor/normal up down
  tumsign = data.frame( matrix(NA,dim(res)[1],length(names2plot)))
  colnames(tumsign) = names(names2plot)
  rownames(tumsign) = rownames(res)
  normsign = tumsign
  valmap = c('greater' = -1, 'less' = 1) # counterinuitive, this is how irina named things
  
  for (metabolite in rownames(d)){
  
    for (cname in colnames(d)){
      if (is.na(d[metabolite,cname])){next} 
      study = strsplit(cname,'\\.')[[1]][1]
      tissue = strsplit(cname,'\\.')[[1]][2]
      testtype = strsplit(cname,'\\.')[[1]][3]
      
      if (study %in% c('all','tumors','normals')){next} #skip these columns
      if (is.na(d[metabolite,"all.combined.adjusted"]) & is.na(d[metabolite,"tumors.combined.adjusted"])){next}
      if ( (d[metabolite,"all.combined.adjusted"] > pthresh | is.na(d[metabolite,"all.combined.adjusted"])) & 
           (d[metabolite,"tumors.combined.adjusted"] > pthresh |  is.na(d[metabolite,"tumors.combined.adjusted"])) ){next} #must be significant by metaanalysis
      
      if (d[metabolite,cname] > pthresh){
        value = 0
      }else{
        value = valmap[testtype]
      }
      
      # Add data to matrix. If value is not NA and not equal to zero, it was significant somewhere else, skip
      if (tissue == 'Normal' & (is.na(normsign[metabolite,study]) | normsign[metabolite,study]==0 )){
        normsign[metabolite,study] = value
      }
      if (tissue == 'Tumor' & ( is.na(tumsign[metabolite,study]) | tumsign[metabolite,study]==0 )){
        tumsign[metabolite,study] = value
      }
      
    }
    
  }
  
  # Identify studies and metabolites where we see consistent changes
  consistent_clin = which(tumsign == normsign & !is.na(tumsign) & tumsign != 0,arr.ind = TRUE)
  consistent_res = data.frame(consistent_clin,tumsign[consistent_clin],normsign[consistent_clin])
  colnames(consistent_res) = c('Row','Column','TumorAssociation','NormalAssociation')
  consistent_res$Metabolite = rownames(tumsign)[consistent_res$Row]
  consistent_res$Study = names2plot[ colnames(tumsign)[consistent_res$Column] ]
  
  write.csv(tumsign,paste0('../results/clinical/ClinicalTumor_SignedChanges_',dname,'.csv'))
  write.csv(normsign,paste0('../results/clinical/ClinicalNormal_SignedChanges_',dname,'.csv'))
  write.csv(consistent_res,paste0('../results/clinical/Tumor_Normal_Consistent_ClinicalChanges_',dname,'.csv'))
}

# Plot joint figure against sample size
pstagegrade = plot_grid(figstage+theme(legend.position = 'none') ,figgrade+theme(legend.position = 'none'),ncol=2,labels = 'AUTO')
save_plot('../results/clinical/SampleSize_JointStageGrade.pdf',pstagegrade,base_height = 4,base_width = 8)
