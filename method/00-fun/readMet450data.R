


readMet450data = function(dataPath = "GDCdata",project = NULL,json = NULL){

  message(paste0("===========================MedBioInfoCloud:Starting==========================="))

  ###从json文件获取信息
  if(is.null(project)& is.null(json))stop("Please set  parameter:project/json")
  if(!is.null(json)){
    metadata_json <- rjson::fromJSON(file=json)
    json_info <- do.call(rbind, lapply(1:length(metadata_json), function(i){
      TCGA_Barcode <- metadata_json[[i]][["associated_entities"]][[1]][["entity_submitter_id"]]
      json_File_Info <- data.frame(cases = TCGA_Barcode)
      rownames(json_File_Info) <- metadata_json[[i]][["file_name"]]
      return(json_File_Info)
    }))
  }else if(!is.null(project)){
    query <- GDCquery(
      project = project,
      data.category = "DNA Methylation",
      data.type = "Methylation Beta Value",
      workflow.type = "SeSAMe Methylation Beta Estimation"
    )
    isServeOK()
    if(missing(query)) stop("Please set query parameter")
    json_info <- query[,1][[1]]
    json_info <- json_info[,c("file_name","cases")]
    rownames(json_info) <- json_info$file_name
    source <- ifelse(query$legacy,"legacy","harmonized")
    # files <- file.path(
    #   query$results[[1]]$project, source,
    #   gsub(" ","_",query$results[[1]]$data_category),
    #   gsub(" ","_",query$results[[1]]$data_type),
    #   gsub(" ","_",query$results[[1]]$file_id),
    #   gsub(" ","_",query$results[[1]]$file_name)
    # )
    # files <- file.path("GDCdata", files)
  }

  filepath = dir(path = dataPath,
                 pattern = "methylation_array.sesame.level3betas.txt$",
                 full.names = T,
                 recursive = T)
  # wd <- filepath[1]
  if(length(filepath)!=0){
    beta <- lapply(filepath,function(wd){
      tempPath <- unlist(strsplit(wd,"/"))
      filename <- tempPath[length(unlist(strsplit(wd,"/")))]
      message(paste0("WeChat:MedBioInfoCloud:Reading:\n",filename))
      oneSampBeta <- read.table(wd,header = F,sep = "\t")
      SampBeta <- data.frame(value = oneSampBeta[,2])
      rownames(SampBeta) <- oneSampBeta[,1]
      colnames(SampBeta) <- json_info[filename,"cases"]
      return(SampBeta)
    })
    si <- unique(unlist(lapply(beta, function(x){nrow(x)})))
    if(length(si) == 1){
      data <- do.call(cbind,beta)
    }else{

      message(paste0("部分文件函数不一样，处理需要更长的时间"))
      mergdat <- NULL
      for(f in 1:length(beta)){
        oneSampBeta <- dplyr::mutate(beta[[f]],cg = rownames(beta[[f]]),.before = 1)
        if(is.null(mergdat)){
          mergdat <- oneSampBeta
        }else{
          mergdat <- merge(mergdat,oneSampBeta,by = "cg",all.x = T)
        }
      }
      rownames(mergdat) <- mergdat$cg
      data <- mergdat[,-1]
    }
    message(paste0("=====================MedBioInfoCloud:Finish======================="))
    return(data)
  }else{message("No related files were maped")}
}

Met450prepare <- function(query){
  isServeOK()
  if(missing(query)) stop("Please set query parameter")
  json_info <- query[,1][[1]]
  json_info <- json_info[,c("file_name","cases")]
  rownames(json_info) <- json_info$file_name
  # source <- ifelse(query$legacy,"legacy","harmonized")
  files <- file.path(
    query$results[[1]]$project,
    gsub(" ","_",query$results[[1]]$data_category),
    gsub(" ","_",query$results[[1]]$data_type),
    gsub(" ","_",query$results[[1]]$file_id),
    gsub(" ","_",query$results[[1]]$file_name)
  )
  filepath <- file.path("GDCdata", files)
  if(length(filepath)!=0){
    beta <- lapply(filepath,function(wd){
      tempPath <- unlist(strsplit(wd,"/"))
      filename <- tempPath[length(unlist(strsplit(wd,"/")))]
      message(paste0("WeChat:MedBioInfoCloud:Reading:\n",filename))
      oneSampBeta <- read.table(wd,header = F,sep = "\t")
      SampBeta <- data.frame(value = oneSampBeta[,2])
      rownames(SampBeta) <- oneSampBeta[,1]
      colnames(SampBeta) <- json_info[filename,"cases"]
      return(SampBeta)
    })
    si <- unique(unlist(lapply(beta, function(x){nrow(x)})))
    if(length(si) == 1){
      data <- do.call(cbind,beta)
    }else{

      message(paste0("Parts of the files have variable numbers of rows and take longer to process"))
      mergdat <- NULL
      for(f in 1:length(beta)){
        oneSampBeta <- dplyr::mutate(beta[[f]],cg = rownames(beta[[f]]),.before = 1)
        if(is.null(mergdat)){
          mergdat <- oneSampBeta
        }else{
          mergdat <- merge(mergdat,oneSampBeta,by = "cg",all.x = T)
        }
      }
      rownames(mergdat) <- mergdat$cg
      data <- mergdat[,-1]
      }

    message(paste0("=====================MedBioInfoCloud:Finish======================="))
    return(data)
  }else{message("No related files were maped")}

}


