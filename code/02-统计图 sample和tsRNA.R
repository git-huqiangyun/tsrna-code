# ===================
library(ggplot2)  
library(dplyr)    
library(tidyr)   

# 设置工作目录到输入文件夹
setwd("D:/procedure/result/")
folders <- list.dirs("D:/2026毕业论文/result", full.names = TRUE, recursive = FALSE)

# 创建空的数据框来存储所有文件夹的信息
all_data <- data.frame(
  TCGA_name = character(),
  sample = numeric(),
  tsRNA = numeric()
)

# 处理每个文件夹
for (folder in folders) {
  # 获取当前文件夹名称
  folder_name <- basename(folder)
  load(file.path(folder, "tsRNA_RPM.Rdata"))
  row_count <- nrow(tsRNA_RPM)  
  col_count <- ncol(tsRNA_RPM) 
  
  # 将数据添加到汇总数据框
  all_data <- rbind(all_data, data.frame(
    TCGA_name = folder_name,
    sample = row_count,
    tsRNA = col_count
  ))
  
  # 打印处理进度
  #cat("已处理文件夹:", folder_name, "\n")
  #cat("行数:", row_count, "\n")
  #cat("列数:", col_count, "\n\n")
}

# 创建第一个图：样本数量统计
p1 <- ggplot(all_data, aes(x = TCGA_name, y = sample)) +
  geom_bar(stat = "identity", fill = "#66C2A5", width = 0.7) +
  geom_text(aes(label = sample), vjust = -0.5, size = 3) +  
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 10),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),  
    plot.title = element_text(hjust = 0.5, size = 12)
  ) +
  labs(
    x = "TCGA ",
    y = "Number of Samples",
    title = "Sample Size TCGA Projects"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# 创建第二个图：tsRNA数量统计
p2 <- ggplot(all_data, aes(x = TCGA_name, y = tsRNA)) +
  geom_bar(stat = "identity", fill = "#FC8D62", width = 0.7) +
  geom_text(aes(label = tsRNA), vjust = -0.5, size = 3) + 
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 10),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    plot.title = element_text(hjust = 0.5, size = 12)
  ) +
  labs(
    x = "TCGA Project",
    y = "Number of tsRNAs",
    title = "tsRNA Count TCGA Projects"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# 保存图形
ggsave(file.path("D:/2026毕业论文/picture-result", "sample_distribution.pdf"), p1, width = 12, height = 6)
ggsave(file.path("D:/2026毕业论文/picture-result", "tsRNA_distribution.pdf"), p2, width = 12, height = 6)



