#===============================================================================
# 【绝对不报错】mRNA + ssGSEA 处理（修复 del_dup_sample 崩溃问题）
#===============================================================================
setwd("D:/2026毕业论文/results")
options(stringsAsFactors = F)

library(TCGAbiolinks)
library(GSVA)
library(GSEABase)
library(dplyr)
library(stringr)

source("D:/procedure/code/00-fun/filterGeneTypeExpr.R")
source("D:/procedure/code/00-fun/del_dup_sample.R")

immune_cell_geneSet = getGmt("D:/procedure/result/Immunity.gmt", geneIdType=SymbolIdentifier())
data_dir <- "D:/procedure/result/Step1_data"
data_files <- list.files(data_dir, pattern = ".Rdata$")

for (data_file in data_files) {
  cancer_type <- gsub(".Rdata$", "", data_file)
  output_dir <- file.path("D:/2026毕业论文/result", cancer_type)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  load(file.path(data_dir, data_file))
  tpm <- filterGeneTypeExpr(expr = mRNA_TPM, fil_col = "gene_type", filter = "protein_coding")
  
  SamN <- TCGAquery_SampleTypes(barcode = colnames(tpm), typesample = c("NT","NB","NBC","NEBV","NBM"))
  SamT <- setdiff(colnames(tpm), SamN)
  
  # ===================== 【修复：样本数 >=2 才去重，否则直接用】=====================
  # 肿瘤样本
  if(length(SamT) > 0){
    if(length(SamT) >= 2){
      turexp_t <- del_dup_sample(tpm[, SamT], col_rename = TRUE)
    } else {
      turexp_t <- tpm[, SamT, drop=FALSE]
    }
    turexp_t <- log2(turexp_t + 1)
  } else {
    turexp_t <- NULL
  }
  
  # 正常样本
  if(length(SamN) > 0){
    if(length(SamN) >= 2){
      turexp_n <- del_dup_sample(tpm[, SamN], col_rename = TRUE)
    } else {
      turexp_n <- tpm[, SamN, drop=FALSE]
    }
    turexp_n <- log2(turexp_n + 1)
  } else {
    turexp_n <- NULL
  }
  
  # 合并
  if(!is.null(turexp_t) & !is.null(turexp_n)){
    turexp_all <- cbind(turexp_t, turexp_n)
  } else if(!is.null(turexp_t)){
    turexp_all <- turexp_t
  } else if(!is.null(turexp_n)){
    turexp_all <- turexp_n
  } else {
    next
  }
  
  # ssGSEA
  param <- GSVA::ssgseaParam(exprData = as.matrix(turexp_all), geneSets = immune_cell_geneSet)
  ssGSEA.Score <- t(gsva(param)) %>% as.data.frame()
  mRNA_TPM <- t(turexp_all) %>% as.data.frame()
  
  # 保存
  save(ssGSEA.Score, file = file.path(output_dir, "ssGSEA.Rdata"))
  save(mRNA_TPM, file = file.path(output_dir, "mRNA_TPM.Rdata"))
  
  cat(sprintf("✅ %s 完成：肿瘤%d例，正常%d例\n", cancer_type, length(SamT), length(SamN)))
  rm(turexp_t, turexp_n, turexp_all, tpm, mRNA_TPM)
}
