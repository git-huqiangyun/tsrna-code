rm(list = ls())
library(plyr)
fp <- paste0("D:/2026毕业论文/results/result-数据集划分/","dataSet.Rdata")

load(fp)# dataSet,testSet,trainSet

head(trainSet)[,1:5]

cong <- colnames(trainSet)[-c(1:grep("OS.Time",colnames(trainSet)))]

#创建输出地址
output_dir <- "D:/2026毕业论文/results/result-单因素Cox回归/"
output_file <- file.path(output_dir, "coxResult.txt")

BaSurv <- Surv(time = trainSet$OS.Time, ## 生存时间
               event = trainSet$OS) ## 生存状态

UniCox <- function(x){ ## 构建一个R function 便于后期调用
  FML <- as.formula(paste0('BaSurv~',x)) ## 构建生存分析公式
  GCox <- coxph(FML, data = trainSet) ## Cox分析
  GSum <- summary(GCox) ## 输出结果
  HR <- round(GSum$coefficients[,2],2) ## 输出HR值
  PValue <- round(GSum$coefficients[,5],3) ## 输出P值
  CI <- paste0(round(GSum$conf.int[,3:4],2),collapse = "-") ## 输出HR的执行区间
  Unicox <- data.frame("characteristics" = x, ## 返回结果，并构建数据框
                       "Hazard Ratio" = HR,
                       "CI95" = CI,
                       HR.95L = GSum$conf.int[,"lower .95"],
                       HR.95H = GSum$conf.int[,"upper .95"],
                       "P Value" = PValue)
  return(Unicox)
}
UniCox(colnames(trainSet)[15]) ## 试验一下 function是否存在错误？

VarNames <- cong ## 输出需要分析的变量名字
UniVar <- lapply(VarNames,UniCox) ## 批量做Cox分析
UniVar <- ldply(UniVar,data.frame) ## 将结果整理为一个数据框
(sigFactors <- UniVar$characteristics[which(UniVar$P.Value < 0.05)] %>% as.character()) ## 筛选其中P值<0.2的变量纳入多因素cox分析。

write.table(UniVar,file=output_file,sep="\t",row.names=F,quote=F)
save(UniVar,sigFactors,file = paste0(output_dir,"coxResult.Rdata"))
writeLines(sigFactors,con = paste0(output_dir,"sigFactors.txt"))
