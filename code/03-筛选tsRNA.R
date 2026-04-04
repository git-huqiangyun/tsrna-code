#============
#基于相关性系数abs>0.8, p<0.05
setwd("D:/procedure/result")
library(WGCNA)
library(dplyr)

# 获取所有癌症类型目录
cancer_dirs <- list.dirs("D:/2026毕业论文/result", recursive = FALSE)

# 循环处理每个癌症类型
for (cancer_dir in cancer_dirs) {
  # 提取癌症类型名称
  cancer_type <- basename(cancer_dir)
  
  # 构建文件路径
  tsRNA_file <- file.path(cancer_dir, "tsRNA_RPM.Rdata")
  ssGSEA_file <- file.path(cancer_dir, "ssGSEA.Rdata")
  
  # 检查文件是否存在
  if (file.exists(tsRNA_file) && file.exists(ssGSEA_file)) {
    # 加载数据
    load(tsRNA_file)
    load(ssGSEA_file)
    
    # 计算相关性矩阵
    cor_matrix <- cor(tsRNA_RPM, ssGSEA.Score, method = "pearson", use = "pairwise.complete.obs")
    
    # 计算样本数
    n_samples <- nrow(tsRNA_RPM)
    
    # 计算P值
    p_matrix <- corPvalueStudent(cor_matrix, n_samples)
    
    # 筛选显著相关的tsRNA
    significant <- (abs(cor_matrix) > 0.2) & (p_matrix < 0.05)
    selected_tsRNA <- rownames(cor_matrix)[rowSums(significant) > 0]
    
    # 处理筛选后的tsRNA数据
    tsRNA_RPM_filtered <- tsRNA_RPM[, selected_tsRNA] %>% 
      t() %>% 
      as.data.frame()
    
    # 创建输出目录
    output_dir <- file.path("D:/2026毕业论文/result", cancer_type)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # 保存结果
    output_file <- file.path(output_dir, "tsRNA_RPM_filtered.Rdata")
    save(tsRNA_RPM_filtered, file = output_file)
    
    # 打印处理信息
    cat(sprintf("\n处理 %s:\n", cancer_type))
    cat("总tsRNA数量:", ncol(tsRNA_RPM), "\n")
    cat("筛选后的tsRNA数量:", length(selected_tsRNA), "\n")
    cat("结果已保存至:", output_file, "\n")
  } 
}


