#===========
# 加载必要的包
library(dplyr)

setwd("D:/procedure/result/")
load("D:/2026毕业论文/result3/merged_tsRNA_RPM.Rdata")
output_dir <- "D:/2026毕业论文/results/result-数据集划分"
clinical_dir <- "D:/2026毕业论文/处理的临床数据" 

# 确保输出目录存在
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# 获取临床数据
clindata_files <- list.files(clinical_dir, pattern = "-survival\\.Rdata$", full.names = TRUE)


# 提取 vitalStat 和 surTime，以 submitter_id 作为行名
clindata_list <- lapply(clindata_files, function(file) {
  load(file)
  clindata <- survival_data
  required_cols <- c("submitter_id", "vitalStat", "surTime")
  if (!all(required_cols %in% colnames(clindata))) {
    stop("文件 ", file, " 缺少 ", paste(required_cols[!required_cols %in% colnames(clindata)], collapse = ", "), " 列！")
  }
  
  # 提取 vitalStat 和 surTime 列，以 submitter_id 作为行名
  result_df <- clindata[, c("vitalStat", "surTime")]
  colnames(result_df) <- c("OS", "OS.Time")
  rownames(result_df) <- clindata$submitter_id  # 设置 submitter_id 为行名
  
  return(result_df)
})

names(clindata_list) <- gsub("-survival\\.Rdata$", "", basename(clindata_files))

# 合并 
merged_clindata <- do.call(rbind, clindata_list)

# 步骤 5: 保存结果
save(merged_clindata, file = file.path(output_dir, "merged_clindata_vitalStat_surTime.Rdata"))

common_samples <- intersect(rownames(merged_df), rownames(merged_clindata))
#cat("Number of common samples:", length(common_samples), "\n")
merged_df_all <- cbind( merged_clindata[common_samples, ], merged_df[common_samples, ])

na_summary <- colSums(is.na(merged_df_all))
print("Number of NA values per column:")
print(na_summary[na_summary > 0])  # 只显示有缺失值的列
merged_df_all[is.na(merged_df_all)] <- 0
# 过滤 OS.Time 为 0 的行
merged_df_all <- merged_df_all[merged_df_all$OS.Time != 0, ]
colnames(merged_df_all) <- gsub("-", "_", colnames(merged_df_all))
#划分测试和训练集
dataSet <- merged_df_all
set.seed(123)  # 设置随机种子以确保可重复性
trainSet <- sample_frac(dataSet,size= 0.7,replace = F) ## 通过sample_frac将数据集按照7:3分开。
testSet <- dataSet[!rownames(dataSet) %in% rownames(trainSet), ] 

save(dataSet, trainSet,testSet, file = file.path(output_dir, "dataSet.Rdata"))
