#============
setwd("D:/procedure/result")
library(GSEABase)
library(dplyr)

# 加载免疫细胞基因集
immune_cell_geneSet = getGmt("D:/procedure/result/Immunity.gmt",
                             geneIdType=SymbolIdentifier())

# 获取所有基因集信息
cat("GeneSetCollection结构检查:\n")
cat("类别:", class(immune_cell_geneSet), "\n")
cat("基因集数量:", length(immune_cell_geneSet), "\n")
cat("基因集名称:", paste(names(immune_cell_geneSet), collapse = ", "), "\n")

# 合并所有免疫相关基因
all_immune_genes <- unique(unlist(lapply(immune_cell_geneSet, geneIds)))

# 获取所有癌症类型目录
cancer_dirs <- list.dirs("D:/2026毕业论文/result", recursive = FALSE)
cancer_dirs
# 循环处理每个癌症类型
for (cancer_dir in cancer_dirs) {
  
  cancer_type <- basename(cancer_dir)
  
  # 构建mRNA_TPM文件路径
  tpm_file <- file.path(cancer_dir, "mRNA_TPM.Rdata")
  
  # 检查文件是否存在
  if (file.exists(tpm_file)) {
    # 加载数据
    load(tpm_file)
    
    # 筛选免疫相关基因
    immune_genes_in_matrix <- intersect(all_immune_genes, rownames(mRNA_TPM))
    
    # 打印匹配信息
    cat(sprintf("\n处理 %s:\n", cancer_type))
    cat("总基因数:", nrow(mRNA_TPM), "\n")
    cat("匹配的免疫相关基因数:", length(immune_genes_in_matrix), "\n")
    
    immune_mRNA_TPM <- mRNA_TPM[rownames(mRNA_TPM) %in% immune_genes_in_matrix, ] %>% 
      t() %>% 
      as.data.frame()
    
    # 保存结果
    output_file <- file.path(cancer_dir, "immune_mRNA_TPM.Rdata")
    save(immune_mRNA_TPM, file = output_file)
  }
}
