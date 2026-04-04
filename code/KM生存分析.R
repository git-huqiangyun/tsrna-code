rm(list = ls())
# 加载必要的包
library(survival)
library(survminer)

fp <- paste0("D:/2026毕业论文/results/result-数据集划分/","dataSet.Rdata")
load(fp)# dataSet,testSet,trainSet

# 定义六个基因
genes <- c("tsRNA_Val_5_0028", "tsRNA_Asn_i_0144", "tsRNA_Gly_5_0056", 
           "tsRNA_iMet_i_0021", "tsRNA_Tyr_5_0013", "tsRNA_Cys_5_0065")
# 选择所需列：OS, OS.Time 和六个基因
selected_data <- dataSet[, c("OS", "OS.Time", genes)]
# 转换天数到月（四舍五入到整数）
selected_data$OS.Time <- round(dataSet$OS.Time / 30.44, 0)

# 查看前几行（类似截图）
head(selected_data)

# 创建输出文件夹
# ==================== 定义输出目录 ====================
opt <- "D:/2026毕业论文/results/KM_curves/"
ifelse(dir.exists(opt),"The folder already exists.",dir.create(opt,recursive = T))

# 循环处理每个基因
for (gene in genes) {
  
  # 确保基因列存在且为数值
  if (!gene %in% colnames(selected_data)) {
    warning(paste("基因", gene, "不存在于数据集中，跳过"))
    next
  }
  
  # ---------- 使用最佳截断值分组----------
  # 利用 surv_cutpoint 寻找最佳分割点（基于最大化 log-rank 统计量）
  res.cut <- surv_cutpoint(selected_data,
                           time = "OS.Time",
                           event = "OS",
                           variables = gene)
  
  # 如果最佳分割点有效，则分类；否则回退到中位数
  if (!is.null(res.cut$cutpoint)) {
    res.cat <- surv_categorize(res.cut)
    selected_data$group <- res.cat[[gene]]  # 自动分为 "high" 和 "low"
  } else {
    selected_data$group <- ifelse(selected_data[[gene]] > median(selected_data[[gene]], na.rm = TRUE),
                             "High", "Low")
  }
  
  # 检查分组是否包含两个水平
  if (length(unique(selected_data$group)) < 2) {
    warning(paste("基因", gene, "分组后只有一个水平，无法绘制KM曲线"))
    next
  }
  
  # 拟合生存曲线
  fit <- survfit(Surv(OS.Time, OS) ~ group, data = selected_data)
  
  # 计算 log-rank 检验的 p 值（用于图例）
  surv_diff <- survdiff(Surv(OS.Time, OS) ~ group, data = selected_data)
  pval <- 1 - pchisq(surv_diff$chisq, df = 1)
  
  # 绘制KM曲线（样式与图片一致）
  plot_title <- paste0(gene, " 的 Kaplan-Meier 曲线")
  p <- ggsurvplot(fit,
                  data = selected_data,
                  pval = TRUE,                     # 显示 p 值
                  pval.method = FALSE,              
                  conf.int = FALSE,                 # 不显示置信区间
                  risk.table = TRUE,                # 显示风险表
                  legend.labs = c("高表达", "低表达"),
                  legend.title = "表达水平",
                  palette = c("orange", "green"),    # 蓝色高表达，黄色低表达
                  xlab = "时间（月）",
                  ylab = "生存概率",
                  title = plot_title,
                  ggtheme = theme_classic(),        # 简洁主题
                  tables.theme = theme_classic(), 
                  risk.table.height = 0.25,)
  
  # 保存为 PDF 文件
  pdf(file.path(opt, paste0(gene, "_KM.pdf")), width = 6, height = 5,family = "GB1")
  print(p, newpage = FALSE)
  dev.off()
  cat("已完成基因", gene, "\n")
}
