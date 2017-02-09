# Visualize clinical analysis results

rm(list = ls())
setwd('/Users/ereznik/pancanmet_analysis/analysis/')

library(XLConnect)
library(ggplot2)
library(cowplot)
library(reshape2)
pthresh = 5e-2
lowestp = 1e-5

# Read in clinical file
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

for (dname in list('grade')){
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
  sigmets = rownames(res)[which(res$Significant_Tumor >= 2)]
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
    res_study[study,'Tumor'] = length( grep(tumname,res$SigNames,ignore.case = TRUE) )
    res_study[study,'Normal'] = length( grep(normname,res$SigNames,ignore.case = TRUE) )
    #res_study[study,'All'] = length( grep(paste(study,'\\.',sep = ''),res$SigNames,ignore.case = TRUE) )
    res_study[study,'Study'] = study
  }
  res_study_melt = melt(res_study,id.vars = 'Study')
  fig3 = ggplot(res_study_melt,aes(Study,value,fill = variable)) + 
    geom_bar(stat = 'identity',position = 'dodge') +
    ggsave('../results/clinical/histogram.pdf',height = 10,width = 10)
  
  # Draw full plot
#   rightcolumn = plot_grid(fig1,fig3,labels = c('B','C'),ncol = 1)
#   plot_grid(fig2,rightcolumn,ncol = 2,scale = c(0.5,1))
  print(fig1)
  print(fig2)
  print(fig3)

  
  write.csv(res,paste('../results/clinical/Visualized_',dname,'.csv',sep = ''))
}
