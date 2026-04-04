# 安装并加载所需包
if (!require("ggVennDiagram")) install.packages("ggVennDiagram")
if (!require("ggplot2")) install.packages("ggplot2")
library(ggVennDiagram)
library(ggplot2)

# 请确保这些变量已正确加载
load("D:/2026毕业论文/results/特征选择-随机森林/rfsrc_selected_genes.Rdata")
load("D:/2026毕业论文/results/result-lasso回归/lasso_geneids.Rdata")
load("D:/2026毕业论文/results/result-单因素Cox回归/coxResult.Rdata")

# 将三个基因集放入列表
gene_lists <- list(
  RandomForest = final_genes,
  LASSO = lasso_geneids,
  Cox = sigFactors
)

# 绘制韦恩图
venn_plot <- ggVennDiagram(gene_lists, 
                           label_alpha = 0,          # 标签背景透明
                           set_size = 4,              # 集合名称大小
                           label_size = 4,            # 数字标签大小
                           edge_size = 1) +           # 边框粗细
  scale_fill_gradient(low = "white", high = "skyblue") +  # 填充颜色
  theme(legend.position = "none")                     # 隐藏图例

# 显示图形
print(venn_plot)

# 保存图形（可选）
ggsave("Venn_three_methods.pdf", plot = venn_plot, width = 6, height = 5)
ggsave("Venn_three_methods.png", plot = venn_plot, width = 6, height = 5, dpi = 300)

# 三者共有
common_three <- Reduce(intersect, gene_lists)

# 随机森林和LASSO共有（不含Cox）
common_rf_lasso <- intersect(final_genes, lasso_geneids)
# 再排除Cox中的基因
common_rf_lasso_only <- setdiff(common_rf_lasso, sig_cox)
