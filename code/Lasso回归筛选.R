rm(list = ls())

library(survival)
library(glmnet) ## 加载R包

fp <- paste0("D:/2026毕业论文/results/result-数据集划分/","dataSet.Rdata")

load(fp)# dataSet,testSet,dataSet

head(dataSet)[,1:5]

#### Step II 迭代特征选取 ####
dataSet <- dataSet[dataSet$OS.Time>0,]
colnames(dataSet)[1:8]
dataSet$OS.Time <- dataSet$OS.Time/30
gs_expr <- as.matrix(dataSet[,-c(1:c(grep("^OS.Time$",colnames(dataSet))))]) ## 根据R包的要求，将数据需要筛选的部分提取转换为矩阵,lasso不需要生存时间和生存状态，所以多删除两列
response <- data.matrix(Surv(dataSet$OS.Time,dataSet$OS))
fit <- glmnet(gs_expr, response,family = "cox" ,maxit = 5000)# 

dat = gs_expr
set.seed(123) 
cvfit = cv.glmnet(x = dat,
                  Surv(dataSet$OS.Time,dataSet$OS), 
                  nfold= 5,#10倍交叉验证，非必须限定条件.
                  family = "cox") 
coef.min = coef(cvfit, s = "lambda.min")  ## lambda.min & lambda.1se 取一个
active.min = which(coef.min != 0 ) ## 找出那些回归系数没有被惩罚为0的
(lasso_geneids <- colnames(dat)[active.min]) ## 提取基因名称

# ==================== 定义输出目录 ====================
opt <- "D:/2026毕业论文/results/result-lasso回归"

# 检查目录是否存在，若不存在则递归创建
if (!dir.exists(opt)) {
  dir.create(opt, recursive = TRUE)
}

# ==================== 保存图片 ====================
plot(cvfit)      # 绘制交叉验证图
plot(fit, label = TRUE)   # 绘制系数路径图

# 保存为 PNG
png("coefficient_path.png", width = 800, height = 600, res = 120)
plot(fit, label = TRUE)   # 绘制系数路径图
dev.off()                  # 关闭图形设备




# ==================== 保存基因列表 ====================
save(lasso_geneids, file = file.path(opt, "lasso_geneids.Rdata"))
writeLines(lasso_geneids, con = file.path(opt, "lasso_geneids_genes.txt"))
