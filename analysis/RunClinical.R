# Clinical analysis v2 (edited by Ed August 2, 2016)
setwd('/Users/ereznik/pancanmet_analysis/analysis/')
d<-read.csv("../data/merged_metabolomics/merged_metabolomics.csv")
cl<-read.csv("../data/merged_metabolomics/clinfeatures.csv")


#start matching IDs between clinical and metabolomics info
cl$type<-as.character(cl$type)
cl$type[cl$type=="Met" ]<-"Metastasis"
cl$sample<-gsub("-",".",cl$sample)
cl$sample<-gsub(" ",".",cl$sample)
cl$sample[substr(cl$sample,nchar(cl$sample),nchar(cl$sample))=="."] <-substr(cl$sample,1,nchar(cl$sample)-1)[substr(cl$sample,nchar(cl$sample),nchar(cl$sample))=="."]
cl$id<-paste(cl$study,cl$sample,cl$type,sep=".")


#There is no clinical information at all on BLCA, PAAD, PAADHussain1 PAADHussain2
d<-d[,!grepl("BLCA",names(d)) & !grepl("PAADHussain1",names(d)) &!grepl("PAADHussain2",names(d)) &!grepl("PAAD",names(d)) ]


ov<-intersect(names(d),cl$id)
d<-d[,c(1,match( ov,names(d)))]
rownames(d)<-d[,1]
d<-d[,-1]
cl<-cl[match( ov,cl$id),]


d<-d[,!is.na(cl$stage) | !is.na(cl$grade)]
cl<-cl[!is.na(cl$stage) | !is.na(cl$grade),]

d<-d[!apply(is.na(d),1,all),] #remove metabolites with no data
d<-d[apply(d,1,function(x){length(table(x))})>1,]#remove metabolites with single value


#end matching IDs between clinical and metabolomics info


library(clinfun)
library(interval)

#We compute a one-sided composite p-values applying Fisher's method to one-sided p-values from individual cohorts. 
#Then, we combine these two into a two-sided p-value using formula 2*min(p1,p2)

fishersMethod = function(x) {if(all(is.na(x)) | sum(!is.na(x))==1) return (NA) 
                             else pchisq(-2 * sum(log(x[!is.na(x)])),df=2*length(x[!is.na(x)]),lower=FALSE)}

combine.pvalues<-
  function(x) {
    x1<-x[1:(length(x)/2)] #left sided p-values
    x2<-x[-(1:(length(x)/2))]#right sided p-values
    if(sum(!is.na(x1))<2) return (NA) else 
    {
      p1<-fishersMethod(x1)
      p2<-fishersMethod(x2)
    return(2*min(p1,p2))
  }}


cl$cohort<-paste(cl$study,cl$type)
study.list.for.stage<-unique(cl$cohort[!is.na(cl$stage)])
study.list.for.grade<-unique(cl$cohort[!is.na(cl$grade)])


#create  matrices of one-sided p-values for stage and grade
ns<-length(study.list.for.stage)
ng<-length(study.list.for.grade)
pmat.stage<-pmat.grade<-matrix(NA,nrow=nrow(d),ncol=20)
set.seed(100)
for (i in 1:nrow(d))
{
  if (any(!is.na(cl$stage) & !is.na(d[i,]))) 
  {
    for (j in 1:length(study.list.for.stage))
    { print(c(i,j));
      w<-!is.na(cl$stage) & !is.na(d[i,]) & cl$cohort==study.list.for.stage[j]
      if (any(w))
        if (mean(d[i,w]==min(d[i,w],na.rm=T))<=0.8) 
        {
       
          if (mean(d[i,w]==min(d[i,w],na.rm=T))>0.2 ) 
            {pmat.stage[i,j]<-ictest(Surv(as.numeric(d[i,w]),as.numeric(d[i,w]), ifelse(as.numeric(d[i,w])!=min(d[i,w]),1,2),type="interval")~cl$stage[w],exact = TRUE,alternative="less")$p.value
             pmat.stage[i,j+ns]<-ictest(Surv(as.numeric(d[i,w]),as.numeric(d[i,w]), ifelse(as.numeric(d[i,w])!=min(d[i,w]),1,2),type="interval")~cl$stage[w],exact = TRUE,alternative="greater")$p.value
          } else
          {
          pmat.stage[i,j]<- jonckheere.test(as.numeric(d[i,w]), cl$stage[w],nperm=5000, alternative="increasing")$p.value
        pmat.stage[i,j+ns]<- jonckheere.test(as.numeric(d[i,w]), cl$stage[w],nperm=5000, alternative="decreasing")$p.value
        }
        }
    }
  }
  
 if (any(!is.na(cl$grade) & !is.na(d[i,]))) 
  {
    for (j in 1:length(study.list.for.grade))
    {
      w<-!is.na(cl$grade) & !is.na(d[i,]) & cl$cohort==study.list.for.grade[j]
     if (any(w))
       if (mean(d[i,w]==min(d[i,w],na.rm=T))<=0.8) 
       {      if (mean(d[i,w]==min(d[i,w],na.rm=T))>0.2 ) 
         {pmat.grade[i,j]<-ictest(Surv(as.numeric(d[i,w]),as.numeric(d[i,w]), ifelse(as.numeric(d[i,w])!=min(d[i,w]),1,2),type="interval")~cl$grade[w],exact = TRUE,alternative="less")$p.value
          pmat.grade[i,j+ng]<-ictest(Surv(as.numeric(d[i,w]),as.numeric(d[i,w]), ifelse(as.numeric(d[i,w])!=min(d[i,w]),1,2),type="interval")~cl$grade[w],exact = TRUE,alternative="greater")$p.value
       } else
       {
              pmat.grade[i,j]<- jonckheere.test(as.numeric(d[i,w]), cl$grade[w],nperm=5000, alternative="increasing")$p.value
         pmat.grade[i,j+ng]<- jonckheere.test(as.numeric(d[i,w]), cl$grade[w],nperm=5000, alternative="decreasing")$p.value
       }
    }
  }
}
}
pmat.stage<-as.data.frame(pmat.stage)
pmat.grade<-as.data.frame(pmat.grade)


names(pmat.stage)<-c(paste(study.list.for.stage,"less",sep="-"),paste(study.list.for.stage,"greater",sep="-") )
names(pmat.grade)<-c(paste(study.list.for.grade,"less",sep="-"),paste(study.list.for.grade,"greater",sep="-") )
rownames(pmat.stage)<-rownames(pmat.grade)<-rownames(d)

pmat.stage<-pmat.stage[,!apply(is.na(pmat.stage),2,all)]
pmat.grade<-pmat.grade[,!apply(is.na(pmat.grade),2,all)]

save(pmat.stage,pmat.grade,file="one.sided.RData")


pmat.stage$all.combined<-apply(pmat.stage,1,combine.pvalues)
pmat.stage$all.combined.adjusted<-p.adjust(pmat.stage$all.combined,method="fdr")

pmat.stage$tumors.combined<-apply(pmat.stage[,grepl("Tumor",names(pmat.stage))],1,combine.pvalues)
pmat.stage$tumors.combined.adjusted<-p.adjust(pmat.stage$tumors.combined,method="fdr")

pmat.stage$normals.combined<-apply(pmat.stage[,grepl("Normal",names(pmat.stage))],1,combine.pvalues)
pmat.stage$normals.combined.adjusted<-p.adjust(pmat.stage$normals.combined,method="fdr")


pmat.grade$all.combined<-apply(pmat.grade,1,combine.pvalues)
pmat.grade$all.combined.adjusted<-p.adjust(pmat.grade$all.combined,method="fdr")

pmat.grade$tumors.combined<-apply(pmat.grade[,grepl("Tumor",names(pmat.grade))],1,combine.pvalues)
pmat.grade$tumors.combined.adjusted<-p.adjust(pmat.grade$tumors.combined,method="fdr")

pmat.grade$normals.combined<-apply(pmat.grade[,grepl("Normal",names(pmat.grade))],1,combine.pvalues)
pmat.grade$normals.combined.adjusted<-p.adjust(pmat.grade$normals.combined,method="fdr")

write.csv(cbind(pmat.stage,pmat.grade),file="updated results stage and grade.csv")


#plots in order of singificance
pdf("updated boxplots in order of singificance with stage.pdf",width=11,height=8.5)

for (i in sort.list(pmat.stage$tumors.combined.adjusted))
  if (!is.na(pmat.stage$tumors.combined.adjusted[i]))
  {
    par(mfrow=c(3,3))
    for (j in study.list.for.stage) 
      try(boxplot(as.numeric(d[i,cl$cohort==j])~as.character(cl$stage)[cl$cohort==j],main=paste(rownames(d)[i],j)))
    title(paste("Tumors only adj.p=",pmat.stage$tumors.combined.adjusted[i]),line=3)
  }
dev.off()


pdf("updated boxplots in order of singificance with grade.pdf",width=11,height=8.5)

for (i in sort.list(pmat.grade$tumors.combined.adjusted))
  if (!is.na(pmat.grade$tumors.combined.adjusted[i]))
  {
    par(mfrow=c(3,3))
    for (j in study.list.for.grade) 
      try(boxplot(as.numeric(d[i,cl$cohort==j])~as.character(cl$grade)[cl$cohort==j],main=paste(rownames(d)[i],j)))
    title(paste("Tumors only adj.p=",pmat.grade$tumors.combined.adjusted[i]),line=3)
  }
dev.off()