# 加载必要的包
library(ggplot2)
library(dplyr)

# 假设 UniVar 已存在，筛选P<0.05的变量（可选）
sig_genes <- UniVar[UniVar$P.Value < 0.05, ]

# 若变量太多，可选取前20个（例如按P值排序）
sig_genes <- sig_genes[order(sig_genes$P.Value), ]  # 按P值升序
sig_genes <- head(sig_genes, 20)                    # 取前20

# 按HR值排序（使图形更易读）
sig_genes <- sig_genes[order(sig_genes$Hazard.Ratio), ]
sig_genes$characteristics <- factor(sig_genes$characteristics, 
                                    levels = sig_genes$characteristics)


output_dir <- "D:/2026毕业论文/results/result-单因素Cox回归/"
output_file <- file.path(output_dir, "coxResult.txt")
# 创建森林图
p_forest <- ggplot(sig_genes, aes(x = Hazard.Ratio, y = characteristics)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_errorbar(aes(xmin = HR.95L, xmax = HR.95H), 
                width = 0.2, linewidth = 0.8, orientation = "y") +  # height 改为 width
  geom_point(size = 3, color = "steelblue") +
  scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8),
                labels = c("0.25", "0.5", "1", "2", "4", "8")) +
  labs(x = "Hazard Ratio (95% CI)", y = NULL, 
       title = "单因素Cox回归森林图") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 10),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank())
# 查看图形
print(p_forest)

# 保存
ggsave(filename = file.path(output_dir, "Univariate_Cox_forestplot.pdf"), 
       plot = p_forest, width = 8, height = 6, device = cairo_pdf)
