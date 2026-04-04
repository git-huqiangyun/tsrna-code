readRNAseq_STAR_Counts2 <- function(query,data_type = c("unstranded","tpm_unstrand")){
  isServeOK()
  if(missing(query)) stop("Please set query parameter")
  
  message(paste0("===========================MedBioInfoCloud:Starting==========================="))
  json_info <- query[,1][[1]]
  json_info <- json_info[,c("file_name","cases")]
  rownames(json_info) <- json_info$file_name
  
  
  source <- ifelse(query$legacy,"legacy","harmonized")
  files <- file.path(
    query$results[[1]]$project, source,
    gsub(" ","_",query$results[[1]]$data_category),
    gsub(" ","_",query$results[[1]]$data_type),
    gsub(" ","_",query$results[[1]]$file_id),
    gsub(" ","_",query$results[[1]]$file_name)
  )
  filepath <- file.path("GDCdata", files)
  
  if(length(filepath)!=0){

    exp <- lapply(filepath,function(wd){
      tempPath <- unlist(strsplit(wd,"/"))
      filename <- tempPath[length(unlist(strsplit(wd,"/")))]
      message(paste0("WeChat:MedBioInfoCloud===========Reading:\n",filename))
      oneSampExp <- read.table(wd,comment.char = "#",header = T,sep = "\t")
      oneSampExp = oneSampExp[-c(1:4),]

      oneSampExp <- dplyr::mutate(oneSampExp,cases = json_info[filename,"cases"],.before =1)
      return(oneSampExp)
      })
    
    df <- do.call(rbind,exp)
    STAR_types <- c("unstranded",
                   "stranded_first",
                   "stranded_second",
                   "tpm_unstranded",
                   "fpkm_unstranded",
                   "fpkm_uq_unstranded")
    data <- list()
   
    # dt = "tpm_unstrand"
    for(dt in data_type){
      if(dt %in% STAR_types){
        dat <- df[,c(colnames(df)[c(1:4)],dt)]
        # dat <- reshape2::dcast(dat,gene_id + gene_name + gene_type ~ cases)
        dat <- tidyr::spread(dat,cases,dt,fill = NA)
        data[[dt]] <- dat
      }
    }
    message(paste0("=====================MedBioInfoCloud:Finish======================="))
    return(data)
  }
}


