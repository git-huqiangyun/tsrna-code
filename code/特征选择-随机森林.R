rm(list = ls())
setwd("D:/procedure/result/")
library(TCGAbiolinks)
library(tidyverse) ## 加载R包
library(survivalsvm) ## 加载R包
library(randomForestSRC)
fp <- paste0("D:/2026毕业论文/results/result-数据集划分/","dataSet.Rdata")
load(fp)# dataSet
colnames(trainSet) <- gsub("-", "_", colnames(trainSet))
head(trainSet)[,1:5]

# ==================== 定义输出目录 ====================
opt <- "D:/2026毕业论文/results/特征选择-随机森林/"
ifelse(dir.exists(opt),"The folder already exists.",dir.create(opt,recursive = T))

#### Step II 迭代特征选取 ####
## Chr I 随机森林（生存分析) ##
# trainSet <- trainSet[,-c(1:c(grep("OS",colnames(trainSet))-1))] ## 只留下生存信息和基因表达量信息，进行特征选取
# testSet <- testSet[,-c(1:c(grep("vitalStat",colnames(testSet))-1))]

rsf <- rfsrc(Surv(OS.Time, OS) ~ . ,
             data = trainSet,
             ntree = 100, ## tree的个数
             nsplit = 1,
             importance = TRUE) ## 随机拆分数（非负整数），可提高运算速度，默认值为0
# 设置树的个数
plot(get.tree(rsf, 5))
# plot(rsf)

print.rfsrc(rsf) ## 输出结果信息

# ---------- 提取 VIMP ----------
vimp <- rsf$importance
# 按降序排序
vimp_sorted <- sort(vimp, decreasing = TRUE)
# 取前 N 个（例如前 10）
topN <- 15
vimp_top <- names(vimp_sorted)[1:min(topN, length(vimp_sorted))]

# ---------- 提取最小深度 ----------
md <- max.subtree(rsf)$order      
# 第一列是 depth（最小深度）
colnames(md) <- c("depth", "") 
head(md)
depth <- md[, "depth"]            # 如果列名不是 "depth"，可改为 md[, 1]

# 最小深度越小，变量越重要。选择深度小于某个阈值的变量
# 取深度最小的前 N 个
# md_top <- names(sort(depth))[1:min(topN, length(depth))]
depth_threshold <- median(depth)   # 或 mean(depth)
md_selected <- names(depth)[depth < depth_threshold]

# ---------- 结合两种方法 ----------
# 取交集（两种方法都认为重要的变量）
final_genes <- intersect(vimp_top, md_selected)

# ---------- 查看结果 ----------
print(final_genes)

# ==================== 保存结果 ====================
save(final_genes, file = paste0(opt, 'rfsrc_selected_genes.Rdata'))
writeLines(final_genes, con = paste0(opt, 'rfsrc_selected_genes.txt'))
