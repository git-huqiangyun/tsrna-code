
##TCGA肿瘤样本去重
del_dup_sample <- function(data,col_rename =T){
  
  data <- data[,sort(colnames(data))]
  pid = gsub("(.*?)-(.*?)-(.*?)-.*","\\1-\\2-\\3",colnames(data))
  library(dplyr)
  
  reps = pid[duplicated(pid)]##重复样本
  if(length(reps)== 0 ){
    message("There were no duplicate patient samples for this data")
  }else{
    data <- data[,!duplicated(pid)]
  }
  if(col_rename == T){
    colnames(data) <- gsub("(.*?)-(.*?)-(.*?)-.*","\\1-\\2-\\3",colnames(data))
  }
  return(data)
}
