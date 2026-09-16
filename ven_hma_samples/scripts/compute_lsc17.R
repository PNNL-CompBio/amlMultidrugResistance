##### NOTE ##########

## Obtained from https://github.com/biodev/beataml2_manuscript/blob/main/R/preprocess.R
## This is the function to compute lsc17 scores used in the Dan + Jeff publication https://pmc.ncbi.nlm.nih.gov/articles/PMC9378589/pdf/nihms-1822675.pdf
## paper name is: Integrative Analysis of Drug Response and Clinical Outcome in Acute Myeloid Leukemia
library(data.table)

compute.lsc17 <- function(exprs){
   
   #from https://www.nature.com/articles/nature20598#Tab8
   #(DNMT3B×0.0874) + (ZBTB46×−0.0347) + (NYNRIN×0.00865) + (ARHGAP22×−0.0138) + (LAPTM4B×0.00582) + (MMRN1×0.0258) + (DPYSL3×0.0284) + 
   #(KIAA0125×0.0196) + (CDK6×−0.0704) + (CPXM1×−0.0258) + (SOCS2×0.0271) + (SMIM24×−0.0226) + (EMP1×0.0146) + (NGFRAP1×0.0465) + 
   #(CD34×0.0338) + (AKR1C3×−0.0402) + (GPR56×0.0501). 
   #As above- and below-median scores in the training cohort were associated with adverse and favourable cytogenetic risk, respectively, a median threshold was used to discretize scores into high and low groups.
   
   lsc17.coefs <- c(DNMT3B=0.0874, ZBTB46=-0.0347, NYNRIN=0.00865, ARHGAP22=-0.0138, LAPTM4B=0.00582, MMRN1=0.0258, DPYSL3=0.0284, 
                    KIAA0125=0.0196, CDK6=-0.0704, CPXM1=-0.0258, SOCS2=0.0271, SMIM24=-0.0226, EMP1=0.0146, NGFRAP1=0.0465, CD34=0.0338, 
                    AKR1C3=-0.0402,GPR56=0.0501)
   
   common.genes <- intersect(names(lsc17.coefs), colnames(exprs))#16
   
   if (length(common.genes) == 16){
      stopifnot(all(setdiff(names(lsc17.coefs), colnames(exprs)) == "SMIM24"))
      
      #smim24 has an ensembl id of ENSG00000095932 so maps to C19orf77
      
      names(lsc17.coefs)[names(lsc17.coefs) == "SMIM24"] <- "C19orf77"
      
      common.genes <- intersect(names(lsc17.coefs), colnames(exprs))#17
      
      stopifnot(length(common.genes) == 17)
   }
   
   # lsc17.score <- exprs[,names(lsc17.coefs)] %*% lsc17.coefs
   lsc17.score <- exprs[, common.genes] %*% lsc17.coefs[common.genes]
   
   #should compute median per cohort
   lsc.dt <- data.table(ptid=rownames(lsc17.score), LSC17=lsc17.score[,1])
   
   lsc.dt
   
}
