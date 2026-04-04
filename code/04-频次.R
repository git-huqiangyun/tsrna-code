# 加载必要的R包
# 设置工作目录到输入文件夹
setwd("D:/procedure/result/")

# 创建输出目录
output_dir <- "D:/2026毕业论文/result3"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# 获取所有文件夹
folders <- list.dirs("D:/2026毕业论文/result", recursive = FALSE)

# 创建一个列表来存储所有行名
all_row_names <- list()

# 遍历每个文件夹
for (folder in folders) {
  # 构建Rdata文件的完整路径
  rdata_file <- file.path(folder, "tsRNA_RPM_filtered.Rdata")
  
  # 检查文件是否存在
  if (file.exists(rdata_file)) {
    # 加载Rdata文件
    load(rdata_file)
    
    # 假设加载的数据框名为tsRNA_RPM_filtered
    # 获取行名并存储
    all_row_names[[basename(folder)]] <- rownames(tsRNA_RPM_filtered)
  }
}

# 保存所有行名到Rdata文件
save(all_row_names, file = file.path(output_dir, "all_tsRNA_row_names.Rdata"))

# 计算每个行名出现的频次
# 将所有行名合并成一个向量
all_unique_row_names <- unique(unlist(all_row_names))

# 创建一个数据框来存储频次统计
frequency_df <- data.frame(
  row_name = all_unique_row_names,
  frequency = sapply(all_unique_row_names, function(x) {
    sum(sapply(all_row_names, function(y) x %in% y))
  })
)
high_frequency_df <- frequency_df[frequency_df$frequency >= 4, ]
save(high_frequency_df, file = file.path(output_dir, "high_frequency_df.Rdata"))
# 按频次降序排序
frequency_df <- frequency_df[order(frequency_df$frequency, decreasing = TRUE), ]

