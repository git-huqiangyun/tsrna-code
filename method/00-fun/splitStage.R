


splitStage <- function(stage){###用来切割"stage_event"列的数据
  clin_se <- data.frame()
  for(i in stage){
    stageT <- ifelse(grepl("T",i),gsub("(T.*?)N.*","\\1",i),"")
    if(grepl("M",i)){
      stageN <- ifelse(grepl("N",i),gsub(paste0(stageT,"(.*?)M.*"),"\\1",i),"")
      stageM <- gsub(paste0(stageT,stageN,"(.*?)"),"\\1",i)
      
    }
    else{
      stageN <- ifelse(grepl("N",i),gsub(paste0(stageT,"(.*?)"),"\\1",i),"")
      stageM <- ""
    }
    
    #stageM <- ifelse(grepl("M",i),gsub(paste0(stageT,stageN,"(.*?)"),"\\1",i),"")
    df <- data.frame(stageT = stageT,
                     stageN = stageN,
                     stageM = stageM)
    clin_se <- rbind(clin_se,df)
  }
  clin_se$stageT <- gsub("[a-d]","",clin_se$stageT)
  clin_se$stageN <- gsub("[a-d]","",clin_se$stageN)
  clin_se$stageM <- gsub("[a-d]","",clin_se$stageM)
  clin_se$stageM <- substr(clin_se$stageM,1,2)
  return(clin_se)
}
get_stage_substr <- function(str,char){#get_stage_substr("M0T2bN0","T")="T2"
  if(grepl(char,str)){
    satr <- str_locate(str,char)
    satr <- as.numeric(satr[1,1])
    stage <- substr(str,satr,satr+1)
  }else{stage <- NA}
  return(stage)
}

splitStage2 <- function(stage){###用来切割"stage_event"列的数据
  clin_se <- data.frame()
  for(i in stage){
    stageT <- get_stage_substr(str=i,char="T")
    stageN <- get_stage_substr(str=i,char="N")
    stageM <- get_stage_substr(str=i,char="M")
    df <- data.frame(stageT = stageT,
                     stageN = stageN,
                     stageM = stageM)
    clin_se <- rbind(clin_se,df)
  }
  return(clin_se)
}