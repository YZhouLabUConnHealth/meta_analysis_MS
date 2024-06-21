library(log4r)
library(vegan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(labdsv)
library(spaa)
library(MASS)


get_log=function(name,message){
  log_file=paste0(name,".log")
  file_logger=logger("WARN",appenders = file_appender(log_file))
  warn(file_logger,paste(message))
  
}


get_feature=function(name){
  print(list.files())
  feature.name=readline("Please input the feature table:\n")
  temp.obj=read.table(feature.name,header = TRUE,sep = '\t',check.names = FALSE,row.names = 1)
  col.count=ncol(temp.obj)
  row.count=nrow(temp.obj)
  message_1=paste(feature.name,"has imported successfully")
  message_col=paste("It has",col.count,"colunms")
  message_row=paste("It has",row.count,"rows")
  get_log(name,message_1)
  get_log(name,message_col)
  get_log(name,message_row)
  return (temp.obj)
}

get_rarefy_all=function(name,object){
  threshold=min(colSums(object))
  get_log(name,paste("Do rarefy by",threshold))
  temp.rarefy=rrarefy(object,threshold)
  proportion.matrix=prop.table(as.matrix(temp.rarefy),2)*100
  
  return(proportion.matrix)
}

input_convert=function(input){
  temp1=strsplit(input," ")
  temp2=as.character(temp1[[1]])
  return(temp2)
  
}
#this function is for extracting matirx
extract_matrix=function(column,object){
  if(column=="0"){
    return(object)
  }
  temp.object=object %>% dplyr::select(-one_of(column))
  return(temp.object)
}

get_rarefy_split=function(name,object){
  temp_del=readline("Please input the column number you want to remove,0 means keep all")
  col_del=input_convert(temp_del)
  message_1.1="Please input the column number you want to remove,0 means keep all"
  message_1.2=col_del
  get_log(name,message_1.1)
  get_log(name,message_1.2)
  #print(col_del)
  #temp.object=object %>% select(-one_of(col_del))
  temp.object=extract_matrix(col_del,object)
  #get_log(name,temp.object)
  print(colSums(temp.object)[order(colSums(temp.object),decreasing = TRUE)])
  
  
  
  
  temp_filter=readline("Please input the column number you want to tranfer proportion directly, use space to input multiples and zero to keep all")
  message_2.1="Please input the column number you want to tranfer proportion directly, use space to input multiples and enter to stop"
  col_filter=input_convert(temp_filter)
  message_2.2=col_filter
  get_log(name,message_2.1)
  get_log(name,message_2.2)
  
   
  #temp.rarefy=temp.object %>% select(-one_of(col_filter))
  if(col_filter=="0"){
    #temp.rarefy=extract_matrix(col_filter,temp.object)
    temp.rarefy=temp.object
    temp.rarefied=get_rarefy_all(name,temp.rarefy)
    #prop.matrix=temp.rarefied
    #proportion.matrix=prop.table(as.matrix(prop.matrix),2)*100
    return(proportion.matrix)
    
  }else{
    temp.rarefy=extract_matrix(col_filter,temp.object)
    temp.rarefied=get_rarefy_all(name,temp.rarefy)
    temp.prop=temp.object %>% dplyr::select(one_of(col_filter))
    prop.matrix=prop.table(as.matrix(temp.prop),2)*100
    proportion.matrix=cbind(prop.matrix,temp.rarefied)
    
    #prop.matrix=cbind(temp.prop,temp.rarefied)
    #proportion.matrix=prop.table(as.matrix(prop.matrix),2)*100
    return(proportion.matrix)
    
  }
    
}



get_rarefy=function(name,object){
  sum_col=colSums(object)[order(colSums(object),decreasing = TRUE)]
  message_1="To rarefy, here is the colsums"
  print(message_1)
  get_log(name,message_1)
  print(sum_col)
  choice=readline("Do you use all samples to do the rarefy directly ?")
  message_2="Do you use all samples to do the rarefy directly?"
  get_log(name,message_2)
  message_3=choice
  get_log(name,message_3)
  temp.object=switch(choice,"yes"=get_rarefy_all(name,object),"no"=get_rarefy_split(name,object))
  #print(temp.object)
  return(temp.object)
  
}


fill_NA=function(object){
  for(i in 1:nrow(object)){
    for(j in 1:ncol(object)){
      if(object[i,j]=='NA'){
        label=switch(j,'k_','p_','c_','o_','f_','g_')
        new.label=paste0(label,"unclassified_",object[i,j-1])
        object[i,j]=new.label
      }
    }
  }
  return(object)
  }


tidy_taxa=function(object){
  for(i in 1:ncol(object)){
    temp=strsplit(as.character(object[,i]),'_')
    label=switch(i,'k_','p_','c_','o_','f_','g_')
    p1=as.character(lapply(temp,'[',3))
    p2=strsplit(p1,";")
    p3=as.character(lapply(p2,'[',1))
    p3[which(p3=='')]=NA
    p3=ifelse(is.na(p3),'NA',paste0(label,p3))
    object[,i]=p3
  }
  
  return(object)
}



#which(as.numeric(as.character(temp1$V7))!='NA')
#temp1$V7[which(as.numeric(as.character(temp.col))!='NA')]='NA'

get_taxa=function(name){
  print(list.files())
  taxa.name=readline("Please input the taxa file name:")
  message_1="Please input the taxa file name:"
  get_log(name,message_1)
  get_log(name,taxa.name)
  taxa.table=read.table(taxa.name,fill = TRUE )
  temp.taxa=taxa.table[-1,]
  rownames(temp.taxa)=temp.taxa[,1]
  temp.taxa=temp.taxa[,-c(1,8,9)]
  #revome numeric values
  temp.col=temp.taxa$V7
  temp.taxa$V7[which(as.numeric(as.character(temp.col))!='NA')]='NA'
  
  taxa.matrix=tidy_taxa(temp.taxa)
  taxa.matrix=fill_NA(taxa.matrix)
  colnames(taxa.matrix)=c('kingdom','phylum','class','order','family','genus')
  write.table(taxa.matrix,paste0(name,'_taxa_new.txt'),sep = '\t',row.names = TRUE,quote = FALSE)
  
  return(taxa.matrix)
  
}

#object1 should be taxa table
merge_table=function(name,object1,object2){
  temp.matrix=merge(object1,object2,by=0)
  rownames(temp.matrix)=temp.matrix[,1]
  #print(temp.matrix)
  genus.matrix=temp.matrix[,-c(1:6)]
  genus.target.matrix=aggregate(.~genus,genus.matrix,sum)
  write.table(genus.target.matrix,paste0(name,'_genus_matrix.txt'),sep = '\t',row.names = FALSE,quote = FALSE)
  
  family.matrix=temp.matrix[,-c(1:5,7)]
  family.target.matrix=aggregate(.~family,family.matrix,sum)
  write.table(family.target.matrix,paste0(name,"_family_matrix.txt"),sep = '\t',row.names = FALSE,quote = FALSE)
  
  order.matrix=temp.matrix[,-c(1:4,6,7)]
  order.target.matrix=aggregate(.~order,order.matrix,sum)
  write.table(order.target.matrix,paste0(name,'_order_matrix.txt'),sep = '\t',row.names = FALSE,quote = FALSE)
  
  class.matrix=temp.matrix[,-c(1:3,5:7)]
  class.target.matrix=aggregate(.~class,class.matrix,sum)
  write.table(class.target.matrix,paste0(name,'_class_matrix.txt'),sep = '\t',row.names = FALSE,quote = FALSE)
  
  phylum.matrix=temp.matrix[,-c(1:2,4:7)]
  phylum.target.matrix=aggregate(.~phylum,phylum.matrix,sum)
  write.table(phylum.target.matrix,paste0(name,'_phylum_matrix.txt'),sep = '\t',row.names = FALSE,quote = FALSE)
  
#  return(target.matrix)
  
}

remove_na_meta=function(meta){

temp.meta=meta %>%mutate_all(na_if,"")


temp.meta.row=subset(temp.meta,rowSums(is.na(temp.meta))!=ncol(temp.meta))

temp.meta.col=subset(temp.meta.row,colSums(is.na(temp.meta.row))!=nrow(temp.meta.row))

return(temp.meta.col)
}




get_meta=function(name){
  meta.name=readline("Please input the meta file:")
  message_1="Please input the meta file:"
  get_log(name,message_1)
  get_log(name,meta.name)
  #should provide multiple separate choice
  meta.matrix=read.table(meta.name,header=TRUE,sep = '\t')
  #meta.matrix=read.table(meta.name,header=TRUE,sep = ',')
  #meta.matrix=remove_na_meta(meta.matrix)
  print(meta.matrix)
  print(colnames(meta.matrix))
  col_group=readline("Please select the target columns:")
  col_temp=input_convert(col_group)
  message_2.1='Please select the target columns:'
  get_log(name,message_2.1)
  get_log(name,col_temp)
  
  
  #the first colunm should be sampleID, the second is the group info.
  group.matrix=meta.matrix %>%dplyr::select(one_of(col_temp))
  #add extract samples based on otu feature table from meta table.  
  
  
  
  

  group.matrix$group=ifelse(group.matrix[,2] %in% c("Remission","Active","MS",'ms','RRMS','CIS','SPMS','PPMS','RRMS','Multiple sclerosis in remission'),'MS','Control')
  colnames(group.matrix)=c('sampleID','disease','group')
  write.table(group.matrix,paste0(name,'_meta_1.txt'),sep = '\t',quote = FALSE)
  #print(group.matrix)
  return(group.matrix)
  
  
}


get_group=function(object,meta,label){

  group.list=intersect(colnames(object),meta$sampleID[meta$group==label])
  return(group.list)
}





get_top25=function(object,meta){
  #rownames(object)=object[,1]
  #object=object[,-1]
  order.object=object[order(rowMeans(object),decreasing = TRUE),]
  top25=order.object[c(1:25),]
  top.list=rownames(top25)
  
  ctrl.list=get_group(top25,meta,"Control")
  print(ctrl.list)
  #write.table(ctrl.list,'ctrl.ist.txt')
  #ctrl.group=top25[,ctrl.list]
  #ctrl.group=top25 %>% dplyr::select(one_of(ctrl.list))
  print(match(ctrl.list,colnames(top25)))
  ctrl.group=top25[,match(ctrl.list,colnames(top25))]
  print(dim(ctrl.group))
  ctrl.plot=data.frame(rowMeans(ctrl.group))
  ctrl.plot$group='Control'
  ctrl.plot$taxa=rownames(ctrl.plot)
  colnames(ctrl.plot)[1]='value'
  
  ms.list=get_group(top25,meta,"MS")
  print(ms.list)
  #ms.group=top25[,ms.list]
  #ms.group=top25 %>% dplyr::select(one_of(ms.list))
  ms.group=top25[,match(ms.list,colnames(top25))]
  print(dim(ms.group))
  ms.plot=data.frame(rowMeans(ms.group))
  ms.plot$group='MS'
  ms.plot$taxa=rownames(ms.plot)
  colnames(ms.plot)[1]='value'
  
  plot25=rbind(ms.plot,ctrl.plot)
  
  COLORS=c("#00FFFF","#000000","#0000FF","#FF00FF","#808080","#008000","#00FF00","#800000","#000080","#808000","#800080","#ff0000","#C0C0C0","#008080","#FFB6C1","#DAA520","#98FB98","#5900FF","#FFDD00","#CAE1FF","#FFE4C4","#33A1C9","#B200FF","#0019FF","#FF9400")
  Genus=factor(plot25$taxa,levels = top.list)
  
  fig.25=ggplot(plot25,aes(x=plot25$group,y=plot25$value,fill=Genus))+geom_bar(stat='identity')+ylab('Relative abundance %')+scale_fill_manual(values=COLORS)+xlab("Group")
  ggsave('top25.png',width = 1920/72,height = 1080/72,dpi = 72)
  #return(fig.25)
  #fig.25
}

get_alpha=function(object,meta){
  #rownames(object)=object[,1]
  #object=object[,-1]
  div.matrix=t(object)
  alpha.div=diversity(div.matrix,'shannon')
  group=meta$group[match(rownames(div.matrix),meta$sampleID)]
  alpha.matrix=cbind(alpha.div,group)
  colnames(alpha.matrix)=c('value','group')
  alpha.matrix=data.frame(alpha.matrix)
  fig.alpha=ggplot(alpha.matrix,aes(x=group,y=as.numeric(value)))+geom_boxplot()+geom_jitter(shape=16,position = position_jitter(w=0.2,h=0.2))+ggtitle('Alpha Diversity at baseline' )+ylab("value")
  #ggsave('alpha_div.pdf')
  #return(fig.alpha)
  #fig.alpha
  ggsave('alpha.png',width = 1920/72,height = 1080/72,dpi = 72)
}

get_richness=function(object,meta){
  #rownames(object)=object[,1]
  #object=object[,-1]
  div.matrix=t(object)
  group=meta$group[match(rownames(div.matrix),meta$sampleID)]
  rich.matrix=specnumber(div.matrix)
  rich.matrix=cbind(rich.matrix,group)
  rich.matrix=data.frame(rich.matrix)
  colnames(rich.matrix)=c('value','group')
  fig.rich=ggplot(rich.matrix,aes(x=group,y=as.numeric(value)))+geom_boxplot()+geom_jitter(shape=16,position = position_jitter(w=0.2,h=0.2))+ggtitle('Richness at baseline' )+ylab("value")
  #ggsave(fig.rich,'rich.pdf')
  #return(fig.rich)
  #fig.rich
  ggsave('rich.png',width = 1920/72,height = 1080/72,dpi = 72)
  
}

get_beta=function(object,meta){
  #rownames(object)=object[,1]
  #object=object[,-1]
  div.matrix=t(object)
  beta.div=dsvdis(div.matrix,index = "bray/curtis",diag = TRUE)
  beta.div.paired=dist2list(beta.div)
  beta.div.paired$col.var=meta$group[match(beta.div.paired$col,meta$sampleID)]
  beta.div.paired$row.var=meta$group[match(beta.div.paired$row,meta$sampleID)]
  beta.div.paired$dist=ifelse(beta.div.paired$col.var==beta.div.paired$row.var,'within','between')
  
  ctrl.beta=subset(beta.div.paired,beta.div.paired$col.var=='Control')
  ctrl.test.res=t.test(ctrl.beta$value~ctrl.beta$dist)
  ctrl.pvalue=ctrl.test.res$p.value
  ctrl.beta.fig=ggplot(ctrl.beta,aes(x=dist,y=value))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16,position = position_jitter(w=0.2,h=0.2))+ggtitle(paste('Control group Beta Diversity with pvalue',ctrl.pvalue ))
  #ggsave('control_beta.pdf')
  ggsave('control_beta.png',width = 1920/72,height = 1080/72,dpi = 72)
  ms.beta=subset(beta.div.paired,beta.div.paired$col.var=='MS')
  ms.test.res=t.test(ms.beta$value~ms.beta$dist)
  ms.pvalue=ms.test.res$p.value
  ms.beta.fig=ggplot(ms.beta,aes(x=dist,y=value))+geom_boxplot(outlier.shape = NA)+geom_jitter(shape=16,position = position_jitter(w=0.2,h=0.2))+ggtitle(paste('MS group Beta Diversity at baseline with pvalue',ms.pvalue ))
  #ggsave('ms_beta.pdf')
  #ms.beta.fig
  ggsave('ms_beta.png',width = 1920/72,height = 1080/72,dpi = 72)
}



get_pca=function(object,meta){
  #rownames(object)=object[,1]
  #object=object[,-1]
  div.matrix=t(object)
  group=meta$group[match(rownames(div.matrix),meta$sampleID)]
  
  pca.dis=vegdist(div.matrix)
  pca.mdso=isoMDS(pca.dis,k=2)
  pca.col=ifelse(group=='MS',"blue","red")
  
  pca.fig=plot(pca.mdso$points,col=pca.col,main = "NMDS",xlab ="NMDS1",ylab = "NMDS2",pch=19 )
  lg=c("MS","Ctrl")
  legend("topleft",col = c("blue","red"),lg,pch = 19)
  
  pca.fig
  
  #ggsave('pca.png',width = 1920/72,height = 1080/72,dpi = 72)
  #ggsave('pca.png')
  
}

get_ggpca=function(object,meta){
  div.matrix=t(object)
  group=meta$group[match(rownames(div.matrix),meta$sampleID)]
  
  pca.dis=vegdist(div.matrix)
  pca.mdso=isoMDS(pca.dis,k=2)
  pca.col=ifelse(group=='MS',"blue","red")

  pca.matrix=data.frame(pca.mdso$data)
  pca.matrix=cbind(pca.matrix,pca.col)
  #pca.matrix=data.frame(pca.matrix)
  print(pca.matrix)
  write.table(pca.matrix,'pca.matrix.txt')
  colnames(pca.matrix)=c('NMDS1','NMDS2','Group')
  ggpca.fig=ggplot(pca.matrix,aes(x='NMDS1',y='NMDS2',color='Group'))+geom_points()
  
  ggsave('ggpca.png',width = 1920/72,height = 1080/72,dpi = 72)
}


get_sharelist=function(object1,object2){
  feature.col=colnames(object1)
  meta.row=as.character(object2$sampleID)
  share.list=intersect(feature.col,meta.row)
  
  return(share.list)
}

get_start=function(){
  name=readline("Please input the project name")
  path=readline("Please input the path:")
  setwd(path)
  
  message_1="Please input the project name:"
  message_2="Please input the path:"

  log_file=paste0(name,".log")
  if(file.exists(log_file)){
    file.remove(log_file)
  }
  get_log(name,message_1)
  get_log(name,name)
  get_log(name,message_2)
  get_log(name,path)
  
  
  taxa.matrix=get_taxa(name)
  feature.matrix=get_feature(name)
  meta.matrix=get_meta(name)
  meta.matrix$sampleID=as.character(meta.matrix$sampleID)

  
  share.list=get_sharelist(feature.matrix,meta.matrix)
  share.count=length(share.list)
  get_log(name,paste(share.count,"shared samples below used for further analysis"))
  get_log(name,share.list)
  
  feature.matrix2=feature.matrix %>% dplyr::select(one_of(share.list))
  write.table(feature.matrix2,paste0(name,"-feature_1.txt"),sep = '\t',quote = FALSE)

  meta.matrix2=meta.matrix[match(meta.matrix$sampleID,share.list),]
  meta.matrix2=meta.matrix[match(meta.matrix[,1],share.list),]
  meta.matrix2=na.omit(meta.matrix2)

  print("This is the meta.matrix2 table")
  print(meta.matrix2)
  write.table(meta.matrix2,"meta_2.txt",sep = '\t',quote = FALSE)
  
  feature.matrix3=get_rarefy(name,feature.matrix2)
  write.table(feature.matrix3,paste0(name,"_prop_matrix.txt"),sep = '\t',row.names = TRUE,quote = FALSE)
  
  
  share.list.new=get_sharelist(feature.matrix3,meta.matrix)
  meta.matrix3=meta.matrix[match(meta.matrix$sampleID,share.list.new),]
  meta.matrix3=na.omit(meta.matrix3)
  write.table(meta.matrix3,paste0(name,'_meta_3.txt'),sep = '\t',quote = FALSE,row.names = FALSE)
  
  #taxa.matrix=get_taxa(name)
  #genus.matrix=merge_table(name,taxa.matrix,feature.matrix3)
  merge_table(name,taxa.matrix,feature.matrix3)
   #write.table(genus.matrix,'genus_matrix.txt',sep = '\t',quote = FALSE,row.names = FALSE)
  #print(temp4.matrix)
  #print("GO on processing")
#  rownames(genus.matrix)=genus.matrix[,1]
#  genus.matrix=genus.matrix[,-1]
  
  
 # get_top25(genus.matrix,meta.matrix)
 # get_alpha(genus.matrix,meta.matrix)
 # get_richness(genus.matrix,meta.matrix)
 # get_beta(genus.matrix,meta.matrix)
  #get_pca(genus.matrix,meta.matrix)
  #get_ggpca(genus.matrix,meta.matrix)
  
  
}

