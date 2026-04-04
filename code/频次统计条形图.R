# 加载高频tsRNA数据
load(file.path(output_dir, "high_frequency_df.Rdata"))

# 筛选频次≥5的tsRNA
high_freq_ge5 <- high_frequency_df[high_frequency_df$frequency >= 5, ]

# 按频次升序排列（绘图时y轴从上到下为从高到低）
high_freq_ge5 <- high_freq_ge5[order(high_freq_ge5$frequency), ]

# 创建水平条形图
p <- ggplot(high_freq_ge5, aes(x = frequency, y = reorder(row_name, frequency))) +
  geom_col(aes(fill = frequency)) +
  scale_fill_gradient(low = "lightblue", high = "darkblue")+
  geom_text(aes(label = frequency), hjust = -0.2, size = 3, color = "black") +
  labs(x = "出现频次", y = "tsRNA", 
       title = "跨癌种高频tsRNA频次分布 (频次≥5)") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 9, family = "sans"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1)))

# 保存图形（使用cairo_pdf支持中文）
pdf(file.path(output_dir, "tsRNA_frequency_barplot_ge5.pdf"), 
          width = 8, height = 6,family = "GB1")  # 高度可根据实际条目数调整（约15条时6英寸合适）
print(p)
dev.off()

# 同时保存为PNG
ggsave(file.path(output_dir, "tsRNA_frequency_barplot_ge5.png"), 
       plot = p, width = 8, height = 6, dpi = 300, bg = "white")
