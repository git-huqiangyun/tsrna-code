#===============================================================================
# 【修复版】泛癌合并（自动对齐所有tsRNA，永不报错）
#===============================================================================
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)

root <- "D:/2026毕业论文/result"
cancer_dirs <- list.dirs(root, recursive = FALSE)

# 用来存所有数据
df_list <- list()

for (cdir in cancer_dirs) {
  cancer <- basename(cdir)
  f <- file.path(cdir, "tsRNA_RPM.Rdata")
  if (!file.exists(f)) next
  
  load(f)
  df <- tsRNA_RPM
  
  # 加入样本ID、癌种、样本类型
  df$SampleID <- rownames(df)
  df$CancerType <- cancer
  df$SampleType <- ifelse(str_sub(df$SampleID,14,15) %in% paste0("0",1:9), "Tumor","Normal")
  
  df_list[[cancer]] <- df
}

# ===================== 【核心修复】自动合并所有列，补齐缺失值 =====================
pan_cancer_data <- bind_rows(df_list)

# 去掉完全空的列
pan_cancer_data <- pan_cancer_data[, !apply(is.na(pan_cancer_data), 2, all)]

# 保存最终矩阵
save(pan_cancer_data, file = file.path(root, "泛癌_tsRNA_肿瘤+正常_最终矩阵.Rdata"))
cat("✅ 泛癌合并成功！\n")
