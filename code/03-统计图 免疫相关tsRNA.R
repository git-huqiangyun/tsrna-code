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
  tsRNA = numeric()
)

# 处理每个文件夹
for (folder in folders) {
  # 获取当前文件夹名称
  folder_name <- basename(folder)
  load(file.path(folder, "tsRNA_RPM_filtered.Rdata"))
  row_count <- nrow(tsRNA_RPM_filtered)  
  
  # 将数据添加到汇总数据框
  all_data <- rbind(all_data, data.frame(
    TCGA_name = folder_name,
    tsRNA = row_count
  ))
  
  # 打印处理进度
  #cat("已处理文件夹:", folder_name, "\n")
  #cat("行数:", row_count, "\n")
  #cat("列数:", col_count, "\n\n")
}

my_colors = c("#D92B03","#F1997B","#F38B2F","#A5405E","#F4C288","#088C00","#025939","#214EA7"
              ,"#8B511F","#DB6C76","#03C088","#7552A7","#DCA0DD","#F2E851","#F1B543","#BF7533"
              ,"#A38277","#592E13","#F2CDCF","#6E86A5","#B7DAFE","#7BB1E3","#00C892","#038766"
              ,"#396251","#E8DC6C",    "#F4A63A","#CFBCD4","#A77D9A","#DF562C","#8A4B43","#1963B3",
              "#B12222","#FF4400","#FEA600","#FFE100","#DEB887","#F5E7BC",  "#2E8A57","#9BCD31"
              ,"#20B3AA","#B0E1E7","#789CF2","#3B63E6")
# 创建第二个图：tsRNA数量统计
p2 <- ggplot(all_data, aes(x = TCGA_name, y = tsRNA, fill = TCGA_name)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = tsRNA), vjust = -0.5, size = 3) + 
  theme_bw() +
  scale_fill_manual(values = my_colors, name = "TCGA Project") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 12),
    legend.position = "bottom",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6),    # 缩小图例文字
    legend.spacing.x = unit(0.2, 'cm'),      # 调整图例项水平间距
    legend.key.size = unit(0.4, 'cm')        # 调整图例色块大小
  ) +
  guides(fill = guide_legend(
    nrow = 2,                # 分为2行
    title.position = "top",  # 标题位置
    label.position = "bottom"# 标签位置
  )) +
  labs(
    y = "tsRNA_RPM_filtered",
    title = "tsRNA Count TCGA Projects"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# 保存图形
ggsave(file.path("D:/2026毕业论文/picture-result", "tsRNA_RPM_filtered.pdf"), p2, width = 12, height = 6)



