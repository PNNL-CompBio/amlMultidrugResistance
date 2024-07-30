##RUn differential expression analysis based on harmonized proteomics from step 0

library(limma)


proteomics <- read.csv('global_proteomics_matrix.csv')
rownames(proteomics) <- proteomics$X
colnames(proteomics) <- sub("^X", "", colnames(proteomics))
columns <- colnames(proteomics)
parental_cols<- columns[grep("Parental", colnames(proteomics))]


meta <- read.csv('meta.csv')
rownames(meta) <- proteomics$name


design <- model.matrix(~0 + factor(meta$name))
colnames(design)<-lapply(colnames(design), function(x) strsplit(x, ")")[[1]][2])
contrast_matrix <- limma::makeContrasts(
  D_only-none_Parental, 
  G_Early-none_Parental, 
  G_Late-none_Parental, 
  
  GD_Early-none_Parental, 
  GD_Late-none_Parental, 
  
  #GV_Early-none_Parental, 
  #GV_Late-none_Parental, 
  
  #GVD_Early-none_Parental, 
  #GVD_Late-none_Parental, 
  levels = design
) 
dim(contrast_matrix)
fit <- limma::lmFit(proteomics, design)
fit <- limma::contrasts.fit(fit, contrast_matrix)
fit <- limma::eBayes(fit)

limma::write.fit(fit, file = "diff_exp_d_only.csv",  adjust = "BH", sep = ",")


for (i in c('D_only', 'G_Early', 'G_Late','GD_Early', 'GD_Late')){#'GV_Early', 'GV_Late', , 'GVD_Early', 'GVD_Late')){
  constrast <- paste(i, ' - none_Parental', sep='')
  res <- limma::topTable(fit, coef=constrast, number=Inf, sort.by="none")
  write.table(res, gsub(' ', '', paste('diff_ex', constrast, '.csv')), sep=',')
  print(constrast)
}

###Phospho results
phospho <- read.csv('phospho_proteomics_matrix.csv')
rownames(phospho) <- phospho$X
colnames(phospho) <- sub("^X", "", colnames(phospho))
columns <- colnames(phospho)
parental_cols<- columns[grep("Parental", colnames(phospho))]


meta <- read.csv('meta.csv')
rownames(meta) <- phospho$name

design <- model.matrix(~0 + factor(meta$name))
colnames(design)<-lapply(colnames(design), function(x) strsplit(x, ")")[[1]][2])
contrast_matrix <- limma::makeContrasts(
  D_only-none_Parental, 
  G_Early-none_Parental, 
  G_Late-none_Parental, 
  
  GD_Early-none_Parental, 
  GD_Late-none_Parental, 
  
  #GV_Early-none_Parental, 
  #GV_Late-none_Parental, 
  
  #GVD_Early-none_Parental, 
  #GVD_Late-none_Parental, 
  levels = design
) 
dim(contrast_matrix)
phospho_fit <- limma::lmFit(phospho, design)
phospho_fit <- limma::contrasts.fit(phospho_fit, contrast_matrix)
phospho_fit <- limma::eBayes(phospho_fit)

limma::write.fit(phospho_fit, file = "diff_exp_phospho.csv",  adjust = "BH", sep = ",")


for (i in c('D_only', 'G_Early', 'G_Late','GD_Early', 'GD_Late')){ #'GV_Late','GV_Early' , 'GVD_Early', 'GVD_Late')){
  #  print(i)
  constrast <- paste(i, ' - none_Parental', sep='')
  res <- limma::topTable(phospho_fit, coef=constrast, number=Inf, sort.by="none")
  write.table(res, gsub(' ', '', paste('diff_ex_phospho_', constrast, '.csv')), , sep=',')
  print(constrast)
}