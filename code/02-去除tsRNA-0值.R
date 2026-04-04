#=========
library(TCGAbiolinks)
library(dplyr)
source("D:/procedure/code/00-fun/del_dup_sample.R")

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
  tpm <- tsRNA_RPM
  
  # 筛选样本
  SamN <- TCGAquery_SampleTypes(barcode = colnames(tpm),
                                typesample = c("NT","NB","NBC","NEBV","NBM"))
  SamT <- setdiff(colnames(tpm), SamN)
  
  # 处理肿瘤样本表达数据
  turexp <- del_dup_sample(tpm[, SamT], col_rename = TRUE)
  turexp <- log2(turexp + 1)
  
  # 转置并筛选方差大于0的基因
  tsRNA_RPM <- as.data.frame(t(turexp))
  tsRNA_RPM <- tsRNA_RPM[, apply(tsRNA_RPM, 2, var) > 0]
  
  # 打印处理信息
  cat(sprintf("\n处理 %s:\n", cancer_type))
  cat("总样本数:", nrow(tpm), "\n")
  cat("处理后基因数:", ncol(tsRNA_RPM), "\n")
  
  # 保存结果
  output_file <- file.path(output_dir, "tsRNA_RPM.Rdata")
  save(tsRNA_RPM, file = output_file)
  cat("结果已保存至:", output_file, "\n")
}
