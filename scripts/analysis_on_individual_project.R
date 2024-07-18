library(log4r)
library(vegan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(labdsv)
library(spaa)
library(MASS)
library(ggpubr)

get_log=function(name,message){
  log_file=paste0("analysis_result/",name,".log")
  file_logger=logger("WARN",appenders = file_appender(log_file))
  warn(file_logger,paste(message))
  
}



get_rarefy_all=function(name,object){
  threshold=min(rowSums(object))
  #get_log(name,paste("Do rarefy by",threshold))
  print(paste("Do rarefy by",threshold))
  temp.rarefy=rrarefy(object,threshold)
  #write.table(temp.rarefy,'analysis_result/temp.rarefy.txt',sep = '\t',quote = FALSE)
  proportion.matrix=prop.table(as.matrix(temp.rarefy),1)*100
  
  return(proportion.matrix)
}

input_convert=function(input){
  temp1=strsplit(input," ")
  temp2=as.character(temp1[[1]])
  return(temp2)
  
}

#this function is for extracting matirx
extract_matrix=function(row,object){
  if(row=="0"){
    return(object)
  }
  #temp.object=object %>% dplyr::select(-one_of(column))
  temp.object=object[-row,]
  return(temp.object)
}

get_rarefy_split=function(name,object){
  temp_del=readline("Please input the row number you want to remove,0 means keep all")
  col_del=input_convert(temp_del)
  message_1.1="Please input the rows number you want to remove,0 means keep all"
  message_1.2=col_del
  #get_log(name,message_1.1)
  #get_log(name,message_1.2)
  #print(col_del)
  #temp.object=object %>% select(-one_of(col_del))
  temp.object=extract_matrix(col_del,object)
  write.table(temp.object,paste0('analysis_result/',name,'_otu_counts_matrix.txt'),sep = '\t',quote = FALSE)
  #get_log(name,temp.object)
  print(colSums(temp.object)[order(colSums(temp.object),decreasing = TRUE)])
  
  
  
  
  temp_filter=readline("Please input the row number you want to tranfer proportion directly, use space to input multiples and zero to keep all")
  message_2.1="Please input the column row you want to tranfer proportion directly, use space to input multiples and enter to stop"
  col_filter=input_convert(temp_filter)
  message_2.2=col_filter
  #get_log(name,message_2.1)
  #get_log(name,message_2.2)
  
  
  #temp.rarefy=temp.object %>% select(-one_of(col_filter))
  if(col_filter=="0"){
    #temp.rarefy=extract_matrix(col_filter,temp.object)
    temp.rarefy=temp.object
    temp.rarefied=get_rarefy_all(name,temp.rarefy)
    #prop.matrix=temp.rarefied
    proportion.matrix=prop.table(as.matrix(prop.matrix),1)*100
    write.table(proportion.matrix,paste0(name,'_proportion.table.txt'),sep = '\t',quote = FALSE)
    return(proportion.matrix)
    
  }else{
    temp.rarefy=extract_matrix(col_filter,temp.object)
    temp.rarefied=get_rarefy_all(name,temp.rarefy)
    temp.prop=temp.object %>% dplyr::select(one_of(col_filter))
    prop.matrix=prop.table(as.matrix(temp.prop),1)*100
    proportion.matrix=rbind(prop.matrix,temp.rarefied)
    write.table(proportion.matrix,paste0(name,'_proportion.table.txt'),sep = '\t',quote = FALSE)
    #prop.matrix=cbind(temp.prop,temp.rarefied)
    #proportion.matrix=prop.table(as.matrix(prop.matrix),2)*100
    return(proportion.matrix)
    
  }
  
}



get_rarefy=function(name,object){
  sum_col=colSums(object)[order(colSums(object),decreasing = TRUE)]
  message_1="To rarefy, here is the colsums"
  print(message_1)
  #get_log(name,message_1)
  print(sum_col)
  choice=readline("Do you use all samples to do the rarefy directly ?")
  message_2="Do you use all samples to do the rarefy directly?"
  #get_log(name,message_2)
  message_3=choice
  #get_log(name,message_3)
  temp.object=switch(choice,"yes"=get_rarefy_all(name,object),"no"=get_rarefy_split(name,object))
  #print(temp.object)
  return(temp.object)
  
}




get_meta=function(name){
  meta.name=readline("Please input the meta file:")
  message_1="Please input the meta file:"
  #get_log(name,message_1)
  #get_log(name,meta.name)
  #should provide multiple separate choice
  meta.matrix=read.table(meta.name,header=TRUE,sep = '\t')
  #meta.matrix=read.table(meta.name,header=TRUE,sep = ',')
  #meta.matrix=remove_na_meta(meta.matrix)
 # print(meta.matrix)
  #print(colnames(meta.matrix))
  #col_group=readline("Please select the target columns:")
  #col_temp=input_convert(col_group)
  #message_2.1='Please select the target columns:'
  #print(message_2.1)
  #get_log(name,message_2.1)
  #get_log(name,col_temp)
  
  
  #the first colunm should be sampleID, the second is the group info.
  #group.matrix=meta.matrix %>%dplyr::select(one_of(col_temp))
  #add extract samples based on otu feature table from meta table.  
  
  
  
  
  
  #group.matrix$group=ifelse(group.matrix[,2] %in% c("Remission","Active","MS",'ms','RRMS','CIS','SPMS','PPMS','RRMS','Multiple sclerosis in remission'),'MS','Control')
  #colnames(group.matrix)=c('sampleID','disease','group')
  #write.table(group.matrix,paste0('analysis_result/',name,'_meta_1.txt'),sep = '\t',quote = FALSE)
  #print(group.matrix)
  return(meta.matrix)
  
  
}



get_group=function(object,meta,label){
  
  group.list=intersect(colnames(object),meta$sampleID[meta$group==label])
  return(group.list)
}

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
  write.table(taxa.matrix,paste0('analysis_result/',name,'_taxa_new.txt'),sep = '\t',row.names = TRUE,quote = FALSE)
  
  return(taxa.matrix)
  
}

get_feature=function(name){
  print(list.files())
  feature.name=readline("Please input the feature table:\n")
  temp.obj=readRDS(feature.name)
  col.count=ncol(temp.obj)
  row.count=nrow(temp.obj)
  message_1=paste(feature.name,"has imported successfully")
  message_col=paste("It has",col.count,"colunms")
  message_row=paste("It has",row.count,"rows")
  print(message_1)
  print(message_col)
  print(message_row)
  #get_log(name,message_1)
  #get_log(name,message_col)
  #get_log(name,message_row)
  return (temp.obj)
}



get_sharelist=function(object1,object2){
  feature.row=rownames(object1)
  meta.row=as.character(object2$sampleID)
  share.list=intersect(feature.row,meta.row)
  
  return(share.list)
}

get_start=function(){
  name=readline("Please input the project name")
  path=readline("Please input the path:")
  setwd(path)
  
  message_1="Please input the project name:"
  message_2="Please input the path:"
  
  #save_Path='analysis_result/'
  
  
  log_file=paste0(name,".log")
  if(file.exists(log_file)){
    file.remove(log_file)
  }
  get_log(name,message_1)
  get_log(name,name)
  get_log(name,message_2)
  get_log(name,path)
  
  
  #taxa.matrix=get_taxa(name)
  feature.matrix=get_feature(name)
  meta.matrix=get_meta(name)
  #meta.matrix$sampleID=as.character(meta.matrix$sampleID)
  
  
  share.list=get_sharelist(feature.matrix,meta.matrix)
  #share.count=length(share.list)
  #get_log(name,paste(share.count,"shared samples below used for further analysis"))
  #get_log(name,share.list)
  
  #feature.matrix2=feature.matrix %>% dplyr::select(one_of(share.list))
  feature.matrix2=data.frame(feature.matrix[share.list,])
  #write.table(feature.matrix2,paste0('analysis_result/',name,"-feature_1.txt"),sep = '\t',quote = FALSE)
  
  
  
  feature.matrix3=get_rarefy(name,feature.matrix2)
  write.table(feature.matrix3,paste0('analysis_result/',name,"_genus_prop_matrix.txt"),sep = '\t',row.names = TRUE,quote = FALSE)
  
  
  share.list.new=get_sharelist(feature.matrix3,meta.matrix)
  meta.matrix2=subset(meta.matrix,meta.matrix$sampleID %in% share.list.new)
  
  write.table(meta.matrix2,paste0('analysis_result/',name,'_meta_2.txt'),sep = '\t',quote = FALSE,row.names = TRUE)
  
  
  genus.matrix=feature.matrix3
  
  
}

