# 清空环境
rm(list = ls())

# 加载必要的包
library(survival)
library(timeROC)
library(ggplot2)

# -------------------- 数据准备 --------------------
fp <- paste0("D:/2026毕业论文/results/result-数据集划分/", "dataSet.Rdata")
load(fp)  # 载入 dataSet

# 定义六个目标基因
genes <- c("tsRNA_Val_5_0028", "tsRNA_Asn_i_0144", "tsRNA_Gly_5_0056",
           "tsRNA_iMet_i_0021", "tsRNA_Tyr_5_0013", "tsRNA_Cys_5_0065")

# 提取所需列，并将OS.Time转换为月（如果原数据是天，则除以30.44取整）
selected_data <- dataSet[, c("OS", "OS.Time", genes)]
selected_data$OS.Time <- round(selected_data$OS.Time / 30.44, 0)  # 天转月，四舍五入取整

# 检查数据是否有缺失，如有可考虑删除或填充（此处简单删除含NA的行）
selected_data <- na.omit(selected_data)

# -------------------- 定义输出目录 --------------------
opt <- "D:/2026毕业论文/results/ROC_curves/"
if (!dir.exists(opt)) dir.create(opt, recursive = TRUE)

# 定义感兴趣的时间点（月）：1年=12月，3年=36月，5年=60月
time_points <- c(12, 36, 60)
time_labels <- c("1年", "3年", "5年")

# -------------------- 循环处理每个基因 --------------------
for (gene in genes) {
  
  # 检查基因列是否存在
  if (!gene %in% colnames(selected_data)) {
    warning(paste("基因", gene, "不存在于数据集中，跳过"))
    next
  }
  
  cat("正在处理基因:", gene, "\n")
  
  # 使用 tryCatch 防止某个基因计算失败导致整个循环中断
  tryCatch({
    
    # ---------- 计算时间依赖ROC ----------
    roc_obj <- timeROC(T = selected_data$OS.Time,
                       delta = selected_data$OS,
                       marker = selected_data[[gene]],
                       cause = 1,                 # 感兴趣的事件（通常1表示死亡）
                       times = time_points,
                       iid = TRUE)                 # 计算置信区间（如果需要）
    
    # 提取AUC值并保留两位小数
    auc_vals <- round(roc_obj$AUC, 2)
    
    # 创建绘图数据框
    roc_data <- data.frame(
      FPR = c(roc_obj$FP[, 1], roc_obj$FP[, 2], roc_obj$FP[, 3]),
      TPR = c(roc_obj$TP[, 1], roc_obj$TP[, 2], roc_obj$TP[, 3]),
      Time = factor(rep(time_labels, each = nrow(roc_obj$FP)),
                    levels = time_labels)
    )
    
    # 构造图例标签（包含AUC值）
    legend_labels <- paste0(time_labels, " (AUC = ", auc_vals, ")")
    
    # ---------- 绘制ROC曲线 ----------
    p <- ggplot(roc_data, aes(x = FPR, y = TPR, color = Time)) +
      geom_line(size = 1) +                        # 曲线
      geom_abline(linetype = "dashed", color = "gray50", size = 0.5) + # 对角线
      scale_color_manual(values = c("1年" = "#E69F00",   # 橙色
                                    "3年" = "#56B4E9",   # 蓝色
                                    "5年" = "#009E73"),  # 绿色
                         labels = legend_labels) +
      labs(title = paste0(gene, " 的时间依赖ROC曲线"),
           x = "False Positive Rate",
           y = "True Positive Rate",
           color = "时间点") +
      theme_classic(base_size = 12) +
      theme(legend.position = c(0.7, 0.3),           # 图例位置
            legend.background = element_rect(fill = "white", color = NA),
            plot.title = element_text(hjust = 0.5))
    
    # ---------- 保存为PDF ----------
    pdf(file.path(opt, paste0(gene, "_ROC.pdf")), width = 6, height = 5, family = "GB1")
    print(p)
    dev.off()
    
    cat("  完成，AUC =", paste(auc_vals, collapse = ", "), "\n")
    
  }, error = function(e) {
    warning(paste("基因", gene, "ROC计算失败:", e$message))
  })
}

cat("所有基因ROC分析完成！结果保存在:", opt, "\n")