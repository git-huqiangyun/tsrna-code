##===========

setwd("D:/procedure/result")
options(stringsAsFactors = F)
library(TCGAbiolinks)
library(GSVA)
library(GSEABase)
library(dplyr)

# 加载必要的函数
source("D:/procedure/code/00-fun/filterGeneTypeExpr.R")
source("D:/procedure/code/00-fun/del_dup_sample.R")

# 加载免疫细胞基因集
immune_cell_geneSet = getGmt("D:/procedure/result/Immunity.gmt",
                             geneIdType=SymbolIdentifier())

# 获取所有数据文件
data_dir <- "D:/procedure/result/Step1_data"
data_files <- list.files(data_dir, pattern = ".Rdata$")

# 循环处理每个数据文件
for (data_file in data_files) {
  # 提取癌症类型（例如：TCGA-LUAD）
  cancer_type <- gsub(".Rdata$", "", data_file)
  
  # 创建输出目录
  output_dir <- file.path("D:/2026毕业论文/result", cancer_type)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 加载数据
  load(file.path(data_dir, data_file))
  
  # 处理表达数据
  tpm <- filterGeneTypeExpr(expr = mRNA_TPM,
                            fil_col = "gene_type",
                            filter = "protein_coding")
  
  # 筛选样本
  SamN <- TCGAquery_SampleTypes(barcode = colnames(tpm),
                                typesample = c("NT","NB","NBC","NEBV","NBM"))
  SamT <- setdiff(colnames(tpm), SamN)
  
  # 处理肿瘤样本表达数据
  turexp <- del_dup_sample(tpm[, SamT], col_rename = TRUE)
  turexp <- log2(turexp + 1)
  
  # 计算ssGSEA分数
  param <- GSVA::ssgseaParam(exprData = as.matrix(turexp), 
                             geneSets = immune_cell_geneSet)
  ssGSEA_Score <- gsva(param)
  ssGSEA.Score <- t(ssGSEA_Score) %>% as.data.frame()
  
  # 计算GSVA分数
  gsva.Score <- gsva(GSVA::ssgseaParam(exprData = as.matrix(turexp), 
                                       geneSets = immune_cell_geneSet)) %>% 
    t() %>% 
    as.data.frame()
  
  # 保存结果
  mRNA_TPM <- as.data.frame(turexp)
  save(ssGSEA.Score, file = file.path(output_dir, "ssGSEA.Rdata"))
  save(mRNA_TPM, file = file.path(output_dir, "mRNA_TPM.Rdata"))
  
  # 打印进度
  cat(sprintf("Processed %s: %d tumor samples, %d normal samples\n",
              cancer_type, length(SamT), length(SamN)))
}
