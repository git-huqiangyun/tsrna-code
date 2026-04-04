library(TCGAbiolinks)
library(stringr)
# library(dplyr)   # 如果使用 dplyr::arrange 建议加载

# 设置工作目录（根据实际路径调整）
getwd()
setwd("D:/2026毕业论文/")

# 定义输出目录
output <- "./01_data/"
if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

# 定义所有要下载的癌症类型
cancer_types <- c("TCGA-ACC", "TCGA-BLCA", "TCGA-CESC", "TCGA-CHOL", "TCGA-COAD",
                  "TCGA-DLBC", "TCGA-ESCA", "TCGA-HNSC", "TCGA-KICH", "TCGA-KIRC", "TCGA-KIRP", 
                  "TCGA-LGG", "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-MESO", "TCGA-OV", "TCGA-PAAD",
                  "TCGA-PCPG", "TCGA-PRAD", "TCGA-READ", "TCGA-SARC", "TCGA-SKCM", "TCGA-STAD",
                  "TCGA-TGCT", "TCGA-THCA", "TCGA-THYM", "TCGA-UCEC", "TCGA-UCS", "TCGA-UVM") # 根据您的研究添加



# 全局设置
options(stringsAsFactors = FALSE)

# 加载自定义函数（只需要执行一次）
source("D:/2026毕业论文/00-fun/splitStage.R")

# 循环处理每一种癌症
for (proj in cancer_types) {
  cat("\n========== 正在处理:", proj, "==========\n")
  
  # 尝试下载并处理，遇到错误跳过当前癌种
  tryCatch({
    
    # 1. 查询临床数据
    query <- GDCquery(project = proj,
                      data.category = "Clinical",
                      data.type = "Clinical Supplement",
                      data.format = "bcr xml")
    
    # 2. 下载数据
    GDCdownload(query)
    # GDCdownload(query, method = "api", files.per.chunk = 20)# 下载速度过慢或出错使用
    
    # 3. 初始化存储列表
    clinicData <- list()
    
    ### ---------- follow_up ----------
    clini_follow_up <- GDCprepare_clinic(query, clinical.info = "follow_up")
    if (!is.null(clini_follow_up) && nrow(clini_follow_up) > 0) {
      clini_follow_up <- dplyr::arrange(clini_follow_up, bcr_patient_barcode)
      clini_follow_up <- clini_follow_up[, c("bcr_patient_barcode",
                                             "vital_status",
                                             "days_to_last_followup",
                                             "days_to_death")]
      clini_follow_up$vital_status <- as.vector(clini_follow_up$vital_status)
      clini_follow_up <- clini_follow_up[clini_follow_up$vital_status != "", ]
      clini_follow_up$vitalStat <- ifelse(clini_follow_up$vital_status == "Alive", 0, 1)
      clini_follow_up$surTime <- ifelse(clini_follow_up$vital_status == "Alive",
                                        clini_follow_up$days_to_last_followup,
                                        clini_follow_up$days_to_death)
      clini_follow_up <- clini_follow_up[!is.na(clini_follow_up$surTime), ]
      clini_follow_up <- dplyr::arrange(clini_follow_up, bcr_patient_barcode, desc(surTime))
      clini_follow_up <- clini_follow_up[!duplicated(clini_follow_up$bcr_patient_barcode), ]
      clinicData[["follow_up"]] <- clini_follow_up
    }
    
    ### ---------- stage_event ----------
    clini_stage_event <- GDCprepare_clinic(query, clinical.info = "stage_event")
    if (!is.null(clini_stage_event)) {
      clini_stage_event <- clini_stage_event[, c("bcr_patient_barcode",
                                                 "pathologic_stage",
                                                 "tnm_categories")]
      clini_stage_event$pathologic_stage <- gsub("[A-D]", "", clini_stage_event$pathologic_stage)
      clini_stage_event <- cbind(clini_stage_event,
                                 splitStage2(as.vector(clini_stage_event$tnm_categories)))
      clini_stage_event <- dplyr::arrange(clini_stage_event, bcr_patient_barcode)
      clini_stage_event <- clini_stage_event[!duplicated(clini_stage_event$bcr_patient_barcode), ]
      clinicData[["stage_event"]] <- clini_stage_event
    }
    
    ### ---------- radiation ----------
    clini_radiation <- GDCprepare_clinic(query, clinical.info = "radiation")
    if (!is.null(clini_radiation) && nrow(clini_radiation) > 0) {
      clinicData[["radiation"]] <- clini_radiation
    }
    
    ### ---------- drug ----------
    clini_drug <- GDCprepare_clinic(query, clinical.info = "drug")
    if (!is.null(clini_drug)) {
      clinicData[["drug"]] <- clini_drug
    }
    
    ### ---------- patient ----------
    clini_patient <- GDCprepare_clinic(query, clinical.info = "patient")
    if (!is.null(clini_patient) && nrow(clini_patient) > 0) {
      # 提取生存信息
      sur <- clini_patient[, c("bcr_patient_barcode",
                               "vital_status",
                               "days_to_last_followup",
                               "days_to_death")]
      sur$vital_status <- as.vector(sur$vital_status)
      sur <- sur[sur$vital_status != "", ]
      sur$vitalStat <- ifelse(sur$vital_status == "Alive", 0, 1)
      sur$surTime <- ifelse(sur$vital_status == "Alive",
                            sur$days_to_last_followup,
                            sur$days_to_death)
      sur <- dplyr::arrange(sur, bcr_patient_barcode, desc(surTime))
      sur <- sur[!duplicated(sur$bcr_patient_barcode), ]
      clinicData[["SurvivalData"]] <- sur
      
      clinicData[["patient"]] <- clini_patient
    }
    
    ### ---------- new_tumor_event (可选) ----------
    # clini_new_tumor_event <- GDCprepare_clinic(query, clinical.info = "new_tumor_event")
    # if (!is.null(clini_new_tumor_event)) clinicData[["new_tumor_event"]] <- clini_new_tumor_event
    
    # 4. 保存数据，文件名包含癌症类型
    save_file <- paste0(output, proj, "-clindata.Rdata")
    save(clinicData, file = save_file)
    cat("已保存:", save_file, "\n")
    
  }, error = function(e) {
    cat("处理", proj, "时发生错误，已跳过。错误信息：", e$message, "\n")
  })
}

 cat("\n全部处理完成！\n")