processed_mRNA_TPM <- function(mRNA_TPM){
  source("D:/代码/00-fun/filterGeneTypeExpr.R")
  source("D:/代码/00-fun/del_dup_sample.R")
  tpm <- mRNA_TPM
  tpm <- filterGeneTypeExpr(expr = tpm,
                          fil_col = "gene_type",
                          filter = 	"protein_coding")#	protein_coding
  ##过滤不表达的基因
  tpm <- tpm[apply(tpm,1,var) !=0,]
  ##正常组织样本ID
  SamN <- TCGAquery_SampleTypes(barcode = colnames(tpm),
                              typesample = c("NT","NB","NBC","NEBV","NBM"))

  ##肿瘤组织样本ID
  SamT <- setdiff(colnames(tpm),SamN)

  ###去除重复样本
  tur_exp <- del_dup_sample(tpm[,SamT],col_rename = T)

  return(as.data.frame(tur_exp))
}


