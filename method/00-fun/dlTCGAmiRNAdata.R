###用于读入TCGA数据库miRNA数据(Isoform Expression Quantification)
dlTCGAmiRNAdata1 <- function(data){#data是一个数据框，第一列是miRNA的名称，第二列是值，第三列是barcode
  colnames(data) <- c("mirName","value","barcode")
  temp <- lapply(unique(data$barcode), FUN = function(r){
    assign(r,data[data$barcode == r,] %>% dplyr::group_by(mirName) %>%
             dplyr::summarize(value= sum(value)))
    
  })
  names(temp) <- unique(data$barcode)
  mirdata <- NULL
  for(id in names(temp)){
    data <- temp[[id]]
    colnames(data) <- c("mirName",id)
    if(is.null(mirdata)){
      mirdata <- data
    }else(mirdata <- merge(mirdata,data,by="mirName",all = TRUE))
    
  }
  rownames(mirdata) <- mirdata$mirName
  mirdata <- mirdata[,-1]
  mirdata[is.na(mirdata)] <- 0
  return(mirdata)
}
dlTCGAmiRNAdata2 <- function(data){
  colnames(data) <- c("mirName","value","barcode")
  mirdata <- NULL
  barcode <- unique(data$barcode)
  for(id in barcode){
    oneSampExp <- data[data[,"barcode"]== id,] %>%
      group_by(mirName) %>%
      summarize(value= sum(value))
    oneSampExp$barcode <- id
    if(is.null(mirdata)){mirdata <- oneSampExp
    }else{mirdata <- rbind(mirdata,oneSampExp)}
  }
  mirdata <- spread(mirdata,barcode,value,fill =0)
  mirdata <- as.data.frame(mirdata)
  rownames(mirdata) <- mirdata$mirName
  mirdata <- mirdata[,-1]
  return(mirdata)
}