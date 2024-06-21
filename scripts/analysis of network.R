rm(list=ls())
path='/path'
setwd(path)

library(devtools)
install_github("zdk123/SpiecEasi")
library(SpiecEasi)
library(dplyr)
library(reshape2)
library(qgraph)


feature=read.table('Union_cutoff.0.1.proj7.prop.txt',sep = '\t')
meta=read.table('MS_union.meta.txt',sep = '\t',header = TRUE)


feature$Group=meta$Group[match(rownames(feature),meta$sampleID)]
feature.ms=subset(feature,feature$Group=='MS')
zero.ms=which(colSums(feature.ms[,-188])==0)
feature.ms=feature.ms[,-zero.ms]
feature.ms=feature.ms[,-171]
write.table(feature.ms,'sparcc.ms.prop.txt',sep = '\t',quote = FALSE)

feature.ctrl=subset(feature,feature$Group=='Control')
zero.ctrl=which(colSums(feature.ctrl[,-188])==0)
feature.ctrl=feature.ctrl[,-188]
write.table(feature.ctrl,'sparcc.ctrl.prop.txt',sep = '\t',quote = FALSE)


ms.pvalue=read.table('sparcc_ms.pvalue.mat.txt',sep = '\t')
ms.cors=read.table('sparcc_ms.correlation.mat.txt',sep = '\t')

ctrl.pvalue=read.table('sparcc_ctrl.pvalue.mat.txt',sep = '\t')
ctrl.cors=read.table('sparcc_ctrl.correlation.mat.txt',sep = '\t')

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

reorder_cor_and_p <- function(cormat, pmat){
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
  pmat <- pmat[hc$order, hc$order]
  list(r = cormat, p = pmat)
}


mat.all=reorder_cor_and_p(ctrl.cors,ctrl.pvalue)
mat.p=mat.all$p
mat.r=mat.all$r

mat.pvalue=get_upper_tri(mat.p)
mat.cor=get_upper_tri(mat.r)

mat.pvalue.pair=melt(as.matrix(mat.pvalue))
mat.cor.pair=melt(as.matrix(mat.cor))

mat.pvalue.pair.edit=na.omit(mat.pvalue.pair)
colnames(mat.pvalue.pair.edit)[3]='pvalue'
mat.cor.pair.edit=na.omit(mat.cor.pair)
colnames(mat.cor.pair.edit)[3]='correlation'

joined.mat=left_join(mat.pvalue.pair.edit,mat.cor.pair.edit,by=c('Var1','Var2'))
joined.mat$fdr=p.adjust(joined.mat$pvalue,method = 'fdr')
joined.mat.edit=subset(joined.mat,joined.mat$Var1!=joined.mat$Var2)
write.table(joined.mat.edit,'sparcc.ctrl.mat.no.cutoff.txt',sep = '\t',quote = FALSE)

# cutoff---pvalue<0.05, abs(correlation)>0.4
#sparcc.plot=joined.mat %>% filter(fdr<0.05)%>%filter(abs(correlation)>0.4)
sparcc.plot=joined.mat %>% filter(fdr<0.05)
sparcc.plot.edit=subset(sparcc.plot,abs(sparcc.plot$correlation)>0.15)
write.table(sparcc.plot.edit,'sparcc.ctrl.mat.p_0.05_cor_0.15_mat.txt',sep = '\t',quote = FALSE)

qgraph(sparcc.plot.edit[,c(1,2,4)], 
       title='MS.network_FULL',
       layout='spring',
       theme='Reddit',
       posCol='light green', 
       negCol ='red',
       directed=FALSE,
       curveALL=TRUE,
       label.cex=0.6,
       label.scale=FALSE)


#--------------------------------------------------------------------------#
full_ms=read.table('sparcc.ms.mat.p_0.05_cor_0.15_mat.txt',sep = '\t')
full_ctrl=read.table('sparcc.ctrl.mat.p_0.05_cor_0.15_mat.txt',sep = '\t')

ctrl.share.ms=semi_join(full_ctrl,full_ms,by=c('Var1','Var2'))
ctrl.uniq.ms=anti_join(full_ctrl,full_ms,by=c('Var1','Var2'))
colnames(ctrl.uniq.ms)[c(1:2)]=c('Var2','Var1')
ctrl.share.ms.p2=semi_join(ctrl.uniq.ms,full_ms,by=c('Var1','Var2'))
ctrl.uniq.ms.p2=anti_join(ctrl.uniq.ms,full_ms,by=c('Var1','Var2'))
ctrl.share.ms.plot=rbind(ctrl.share.ms,ctrl.share.ms.p2)

write.table(ctrl.share.ms.plot,'network_ctrl_share_ms.table.txt',sep = '\t',quote = FALSE)
write.table(ctrl.uniq.ms.p2,'network_ctrl_uniq.table.txt',sep = '\t',quote = FALSE)

qgraph(ctrl.uniq.ms.p2[,c(1,2,4)], 
       title='ctrl.network_uniq',
       layout='spring',
       theme='Reddit',
       posCol='light green', 
       negCol ='red',
       directed=FALSE,
       curveALL=TRUE,
       label.cex=0.6,
       label.scale=FALSE)

ms.share.ctrl=semi_join(full_ms,full_ctrl,by=c('Var1','Var2'))
ms.uniq.ctrl=anti_join(full_ms,full_ctrl,by=c('Var1','Var2'))
colnames(ms.uniq.ctrl)[c(1:2)]=c('Var2','Var1')
ms.share.ctrl.p2=semi_join(ms.uniq.ctrl,full_ctrl,by=c('Var1','Var2'))
ms.uniq.ctrl.p2=anti_join(ms.uniq.ctrl,full_ctrl,by=c('Var1','Var2'))
ms.share.ctrl.plot=rbind(ms.share.ctrl,ms.share.ctrl.p2)

write.table(ms.share.ctrl.plot,'network_ms_share_ctrl.table.txt',sep = '\t',quote = FALSE)
write.table(ms.uniq.ctrl.p2,'network_ms_uniq.table.txt',sep = '\t',quote = FALSE)

qgraph(ms.share.ctrl.plot[,c(1,2,4)], 
       title='ms.network_share',
       layout='spring',
       theme='Reddit',
       posCol='light green', 
       negCol ='red',
       directed=FALSE,
       curveALL=TRUE,
       label.cex=0.6,
       label.scale=FALSE)
