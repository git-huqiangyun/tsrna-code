readRNAseq_STAR_Counts <- function(query,data_type = c("unstranded","tpm_unstranded")){
  isServeOK()
  if(missing(query)) stop("Please set query parameter")
  
  message(paste0("===========================MedBioInfoCloud:Starting==========================="))
  
  source <- ifelse(query$legacy,"legacy","harmonized")
  files <- file.path(
    query$results[[1]]$project, source,
    gsub(" ","_",query$results[[1]]$data_category),
    gsub(" ","_",query$results[[1]]$data_type),
    gsub(" ","_",query$results[[1]]$file_id),
    gsub(" ","_",query$results[[1]]$file_name)
  )
  files <- file.path("GDCdata", files)
  cases <- ifelse(
    grepl("TCGA|TARGET|CGCI-HTMCP-CC",query$results[[1]]$project %>% unlist()),
    query$results[[1]]$cases,
    query$results[[1]]$sample.submitter_id
  )
  # read files that has 4 not necessary rows, and has several columns
  # gene_id gene_name gene_type
  # unstranded stranded_first stranded_second tpm_unstranded fpkm_unstranded
  x <- plyr::alply(files,1, function(f) {
    data.table::fread(f)
  }, .progress = "time")
  
  df <- data.table::rbindlist(
    x, use.names = TRUE, idcol = "case_barcode"
  )
  if(!missing(cases))  {
    df$case_barcode <- factor(
      cases[df$case_barcode %>% as.numeric()],
      levels = cases
    )
  }
  # this part changes the cases order if not a factor
  df <- data.table::dcast(
    data = df,
    formula = gene_id + gene_name + gene_type ~ case_barcode,
    value.var = colnames(df)[-c(1:4)]
  )
  
  df <- as.data.frame(df)
  df <- df[grep("^ENSG",df$gene_id),]
  
  
  STAR_types <- c("unstranded",
                  "stranded_first",
                  "stranded_second",
                  "tpm_unstranded",
                  "fpkm_unstranded",
                  "fpkm_uq_unstranded")
  
  data = list()
  for(dt in data_type){
    if(dt %in% STAR_types){
      dat <- df[,(colnames(df)[c(1:3,grep(paste0("^",dt),colnames(df)))])]
      colnames(dat) <- gsub(paste0(dt,"_"),"",colnames(dat))
      data[[dt]] <- dat
    }
  }
  message(paste0("=====================MedBioInfoCloud:Finish======================="))
  return(data)
 
  
}


