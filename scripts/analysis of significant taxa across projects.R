rm(list=ls())
path='/path/'
setwd(path)
library(reshape2)
library(dplyr)
library(coin)
library(nlme)
library(dplyr)
library(lme4)



feature.raw=read.table('Union_raw.prop.proj7.txt',sep = '\t',header = TRUE)
feature.cutoff.edit=feature.raw %>% select_if(colnames(feature.raw) %in% colnames(feature.cutoff))

meta=read.table('MS_union.meta.txt',sep = '\t',header = TRUE)

work.mat=feature.cutoff.edit
work.mat$Group=meta$Group[match(rownames(work.mat),meta$sampleID)]
work.mat$Project=meta$Project[match(rownames(work.mat),meta$sampleID)]
work.mat$Region=meta$Region[match(rownames(work.mat),meta$sampleID)]
work.mat$KIT=meta$KIT[match(rownames(work.mat),meta$sampleID)]
work.mat$Platform=meta$Platform[match(rownames(work.mat),meta$sampleID)]
work.mat$Country=meta$Country[match(rownames(work.mat),meta$sampleID)]
work.mat$Gender=meta$Gender[match(rownames(work.mat),meta$sampleID)]
work.mat$Age_group=meta$Age_group[match(rownames(work.mat),meta$sampleID)]
work.mat$BMI_group=meta$BMI_group[match(rownames(work.mat),meta$sampleID)]

# MS & Ctrl together
# Group test
# 267 257
w.group.pvalue={}
for(i in 1:187){
  temp=wilcox_test(work.mat[,i]~as.factor(work.mat$Group)|as.factor(work.mat$Project))
  w.group.pvalue[i]=pvalue(temp)
}

w.group.pvalue.adj=p.adjust(w.group.pvalue,method = 'fdr')
wilcox.p.value=data.frame(cbind(w.group.pvalue,w.group.pvalue.adj))
rownames(wilcox.p.value)=colnames(feature.cutoff.edit)

w.gender.pvalue={}
for(i in 1:187){
  temp=kruskal_test(work.mat[,i]~as.factor(work.mat$Gender)|as.factor(work.mat$Project))
  w.gender.pvalue[i]=pvalue(temp)
}

w.gender.pvalue.adj=p.adjust(w.gender.pvalue,method = 'fdr')
wilcox.p.value=cbind(wilcox.p.value,w.gender.pvalue,w.gender.pvalue.adj)

#BMI, for colmeans(taxa)==0, NAN
BMI.mat=subset(work.mat,!is.na(work.mat$BMI_group))
w.BMI.pvalue={}
for(i in 1:187){
  temp=kruskal_test(BMI.mat[,i]~as.factor(BMI.mat$BMI_group)|as.factor(BMI.mat$Project))
  w.BMI.pvalue[i]=pvalue(temp)
}

w.BMI.pvalue.adj=p.adjust(w.BMI.pvalue,method = 'fdr')
wilcox.p.value=cbind(wilcox.p.value,w.BMI.pvalue,w.BMI.pvalue.adj)

w.region.pvalue={}
for(i in 1:187){
  temp=kruskal_test(work.mat[,i]~as.factor(work.mat$Region)|as.factor(work.mat$Project))
  w.region.pvalue[i]=pvalue(temp)
}

w.region.pvalue.adj=p.adjust(w.region.pvalue,method = 'fdr')
wilcox.p.value=cbind(wilcox.p.value,w.region.pvalue,w.region.pvalue.adj)


w.KIT.pvalue={}
for(i in 1:187){
  temp=kruskal_test(work.mat[,i]~as.factor(work.mat$KIT)|as.factor(work.mat$Project))
  w.KIT.pvalue[i]=pvalue(temp)
}

w.KIT.pvalue.adj=p.adjust(w.KIT.pvalue,method = 'fdr')
wilcox.p.value=cbind(wilcox.p.value,w.KIT.pvalue,w.KIT.pvalue.adj)


w.country.pvalue={}
for(i in 1:187){
  temp=kruskal_test(work.mat[,i]~as.factor(work.mat$Country)|as.factor(work.mat$Project))
  w.country.pvalue[i]=pvalue(temp)
}

w.country.pvalue.adj=p.adjust(w.country.pvalue,method = 'fdr')
wilcox.p.value=cbind(wilcox.p.value,w.country.pvalue,w.country.pvalue.adj)


w.platform.pvalue={}
for(i in 1:187){
  temp=wilcox_test(work.mat[,i]~as.factor(work.mat$Platform)|as.factor(work.mat$Project))
  w.platform.pvalue[i]=pvalue(temp)
}

w.platform.pvalue.adj=p.adjust(w.platform.pvalue,method = 'fdr')
wilcox.p.value=cbind(wilcox.p.value,w.country.pvalue,w.country.pvalue.adj)

write.table(wilcox.p.value,'./meta_analysis/pvalue.wilcox.ctrl&ms.fdr.txt',sep = '\t',quote = FALSE)

# lme method
lme.mat=log(feature.cutoff.edit+1)
lme.mat$Group=meta$Group[match(rownames(lme.mat),meta$sampleID)]
lme.mat$Project=meta$Project[match(rownames(lme.mat),meta$sampleID)]
lme.mat$Region=meta$Region[match(rownames(lme.mat),meta$sampleID)]
lme.mat$KIT=meta$KIT[match(rownames(lme.mat),meta$sampleID)]
lme.mat$Platform=meta$Platform[match(rownames(lme.mat),meta$sampleID)]
lme.mat$Country=meta$Country[match(rownames(lme.mat),meta$sampleID)]

lme.mat$Gender=meta$Gender[match(rownames(lme.mat),meta$sampleID)]
lme.mat$Age_group=meta$Age_group[match(rownames(lme.mat),meta$sampleID)]
lme.mat$BMI_group=meta$BMI_group[match(rownames(lme.mat),meta$sampleID)]

fit_lme=list()
l.age.pvalue={}
l.gender.pvalue={}
l.group.pvalue={}

# group test
# Control      MS 
# 267     257

for(i in 1:187){
  taxa=lme.mat[,i]
  fit_lme[[i]]=lme(taxa~Group,random = ~1|Project, data = lme.mat)
  temp=data.frame(anova(fit_lme[[i]]))
  l.group.pvalue[i]=temp[2,4]
  
}
l.group.pvalue.adj=p.adjust(l.group.pvalue,method = 'fdr')
lme.p.value.group=data.frame(cbind(l.group.pvalue,l.group.pvalue.adj))
rownames(lme.p.value.group)=colnames(feature.cutoff.edit)
write.table(lme.p.value.group,'./meta_analysis/pvalue.lme.group.ctrl&ms.txt',sep = '\t',quote = FALSE)

#BMI test
#      lean      obese overweight
#       146         71         68

lme.BMI.mat=subset(lme.mat,!is.na(lme.mat$BMI_group))
lme.BMI.mat.edit=lme.BMI.mat[,-which(colMeans(lme.BMI.mat[,c(1:187)])==0)]

fit_lme=list()
l.bmi.pvalue={}

for(i in 1:186){
  taxa=lme.BMI.mat.edit[,i]
  fit_lme[[i]]=lme(taxa~BMI_group,random = ~1|Project, data = lme.BMI.mat.edit)
  temp=data.frame(anova(fit_lme[[i]]))
  l.bmi.pvalue[i]=temp[2,4]
}


l.bmi.pvalue.adj=p.adjust(l.bmi.pvalue,method = 'fdr')

lme.p.value.bmi=data.frame(cbind(l.bmi.pvalue,l.bmi.pvalue.adj))
rownames(lme.p.value.bmi)=colnames(lme.BMI.mat.edit[,c(1:186)])
write.table(lme.p.value.bmi,'pvalue.lme.bmi.ctrl&ms.txt',sep = '\t',quote = FALSE)

#Gender test
#female    male unknown
#320     147      57
fit_lme=list()
l.gender.pvalue={}
for(i in 1:187){
  taxa=lme.mat[,i]
  fit_lme[[i]]=lme(taxa~Gender,random = ~1|Project, data = lme.mat)
  temp=data.frame(anova(fit_lme[[i]]))
  l.gender.pvalue[i]=temp[2,4]
  
}
l.gender.pvalue.adj=p.adjust(l.gender.pvalue,method = 'fdr')
lme.p.value.gender=data.frame(cbind(l.gender.pvalue,l.gender.pvalue.adj))
rownames(lme.p.value.gender)=colnames(feature.cutoff.edit)
write.table(lme.p.value.gender,'pvalue.lme.gender.ctrl&ms.txt',sep = '\t',quote = FALSE)

#platform test
#     454 Illumina 
#     119      405

fit_lme=list()
l.platform.pvalue={}
for(i in 1:187){
  taxa=lme.mat[,i]
  fit_lme[[i]]=lme(taxa~Platform,random = ~1|Project, data = lme.mat)
  temp=data.frame(anova(fit_lme[[i]]))
  l.platform.pvalue[i]=temp[2,4]
  
}
l.platform.pvalue.adj=p.adjust(l.platform.pvalue,method = 'fdr')
lme.p.value.platform=data.frame(cbind(l.platform.pvalue,l.platform.pvalue.adj))
rownames(lme.p.value.platform)=colnames(feature.cutoff.edit)
write.table(lme.p.value.platform,'pvalue.lme.platform.ctrl&ms.txt',sep = '\t',quote = FALSE)

# country test
#   China Japan   USA 
#      68    52   404
fit_lme=list()
l.country.pvalue={}
for(i in 1:187){
  taxa=lme.mat[,i]
  fit_lme[[i]]=lme(taxa~Country,random = ~1|Project, data = lme.mat)
  temp=data.frame(anova(fit_lme[[i]]))
  l.country.pvalue[i]=temp[2,4]
  
}
l.country.pvalue.adj=p.adjust(l.country.pvalue,method = 'fdr')
lme.p.value.country=data.frame(cbind(l.country.pvalue,l.country.pvalue.adj))
rownames(lme.p.value.country)=colnames(feature.cutoff.edit)
write.table(lme.p.value.country,'pvalue.lme.country.ctrl&ms.txt',sep = '\t',quote = FALSE)

# region test
# V1-V2 V1-V3 V3-V4 V3-V5    V4 
# 52    50    68    67       287 
fit_lme=list()
l.region.pvalue={}
for(i in 1:187){
  taxa=lme.mat[,i]
  fit_lme[[i]]=lme(taxa~Region,random = ~1|Project, data = lme.mat)
  temp=data.frame(anova(fit_lme[[i]]))
  l.region.pvalue[i]=temp[2,4]
  
}
l.region.pvalue.adj=p.adjust(l.region.pvalue,method = 'fdr')
lme.p.value.region=data.frame(cbind(l.region.pvalue,l.region.pvalue.adj))
rownames(lme.p.value.region)=colnames(feature.cutoff.edit)
write.table(lme.p.value.region,'pvalue.lme.region.ctrl&ms.txt',sep = '\t',quote = FALSE)

# KIT test
# In_house    MoBio   QIAamp 
# 52          404       68

fit_lme=list()
l.kit.pvalue={}
for(i in 1:187){
  taxa=lme.mat[,i]
  fit_lme[[i]]=lme(taxa~KIT,random = ~1|Project, data = lme.mat)
  temp=data.frame(anova(fit_lme[[i]]))
  l.kit.pvalue[i]=temp[2,4]
  
}
l.kit.pvalue.adj=p.adjust(l.kit.pvalue,method = 'fdr')
lme.p.value.kit=data.frame(cbind(l.kit.pvalue,l.kit.pvalue.adj))
rownames(lme.p.value.kit)=colnames(feature.cutoff.edit)
write.table(lme.p.value.kit,'pvalue.lme.kit.ctrl&ms.txt',sep = '\t',quote = FALSE)

# combine lme results
lme.result=ls(pattern = 'lme.p.value.*')
for(df in lme.result){
  temp=get(df)
  temp$Genus=rownames(temp)
  assign(df,temp)
  rm(temp)
}

temp.lme=get(lme.result[1])
for(df in lme.result[-1]){
  temp.lme=full_join(temp.lme,get(df),by='Genus')
}

rownames(temp.lme)=temp.lme$Genus
temp.lme=temp.lme[,-3]
write.table(temp.lme,'pvalue.lme.ctrl&ms.fdr.txt',sep = '\t',quote = FALSE)
lme.p.value=temp.lme
