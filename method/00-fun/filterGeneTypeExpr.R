

#filter :"protein_coding","lncRNA","miRNA","misc_RNA","snRNA","miRNA","scRNA",....
filterGeneTypeExpr <- function(expr,fil_col = "gene_type",filter = FALSE){
  
  ##Delete unexpressed genes(rows)in all samples from expr.
  dat <- expr[apply(expr[,-c(1:3)], 1, var)!=0,]
  ##rowSums
  dat <- dplyr::mutate(dat,Sums = rowSums(dat[,-c(1:3)]),.before = 4)
  dat <- dplyr::arrange(dat,gene_name,desc(Sums))
  dat <- dat[!duplicated(dat$gene_name),]
  rownames(dat) <- dat$gene_name
  
  if (filter == FALSE) {
    return(dat[,-c(1:4)])
    
  }else{
    dat <- dat[dat[,fil_col] == filter,-c(1:4)]
    return(dat)
  }
}



