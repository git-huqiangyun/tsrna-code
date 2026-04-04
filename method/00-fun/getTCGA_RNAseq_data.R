
getTCGA_RNAseq_data = function(dataPath,project = NULL,json = NULL,data_type){
  
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
    query <- TCGAbiolinks::GDCquery(
      project = project,
      data.category = "Transcriptome Profiling",
      data.type = "Gene Expression Quantification",
      workflow.type = "STAR - Counts"
    )
    json_info <- query[,1][[1]]
    
    json_info <- json_info[,c("file_name","cases")]
    rownames(json_info) <- json_info$file_name
    
  }
  
  ##读取转录组数据
  filepath = dir(path = dataPath,
                 pattern = "counts.tsv$",
                 full.names = T,
                 recursive = T)
  if(length(filepath)!=0){
    exp <- lapply(filepath,function(wd){
      tempPath <- unlist(strsplit(wd,"/"))
      filename <- tempPath[length(unlist(strsplit(wd,"/")))]
      message(paste0("WeChat:MedBioInfoCloud:Reading:\n",filename))
      oneSampExp <- read.table(wd,comment.char = "#",header = T,sep = "\t")
      oneSampExp = oneSampExp[-c(1:4),]
      if(wd == filepath[1]){
        oneSampExp = oneSampExp[,c("gene_id","gene_name","gene_type",data_type)]
        colnames(oneSampExp) <- c("gene_id","gene_name","gene_type",
                                  json_info[filename,"cases"])
        rownames(oneSampExp) <- oneSampExp[,"gene_id"]
      }else{
        rnm = oneSampExp[,"gene_id"]
        oneSampExp = data.frame(value = oneSampExp[,data_type])
        colnames(oneSampExp) <- json_info[filename,"cases"]
        rownames(oneSampExp) <-  rnm
      }
      return(oneSampExp)
    })
    data <- do.call(cbind,exp)
    return(data)
  }else{
    message("Please check the datapath")
  }
}
