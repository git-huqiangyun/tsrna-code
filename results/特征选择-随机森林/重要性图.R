
imp <- rsf$importance   # 命名向量

# 定义你想展示的基因
# my_genes <- c("tsRNA_Asn_i_0144", "tsRNA_Glu_5_0049",
#              "tsRNA_Gly_5_0056","tsRNA_iMet_i_0021", "tsRNA_Glu_3_0002")

my_genes <- final_genes #随机森林筛选的基因

# 提取这些基因的重要性（确保基因名存在于 imp 中）
imp_selected <- imp[my_genes]

# 转换为数据框
imp_df <- data.frame(
  gene = names(imp_selected),
  importance = as.numeric(imp_selected)
)

# 按重要性排序（可选，若你想按基因列表顺序排列，可设置因子水平）
imp_df$gene <- factor(imp_df$gene, levels = my_genes)

# 绘制条形图
p_imp <- ggplot(imp_df, aes(x = importance, y = gene)) +
  geom_col(fill = "steelblue") +
  labs(x = "Variable Importance (VIMP)", y = NULL,
       title = "核心tsRNA的随机森林重要性") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
print(p_imp)
png("importance_plot.png", width = 2400, height = 1800, res = 300)
print(p_imp)
dev.off()







