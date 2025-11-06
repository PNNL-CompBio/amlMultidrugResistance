# figure s1: DIA vs. TMT and bead vs. flow: CD14+ vs. CD34+
library(ggplot2)
# load sorted proteomics signature
# sig.paths <- list("DIA: bead" = "analysis/combined24-27/DIA_2batches_noOutliers_noMSC/Sort Type_Bead/CD14_Pos_vs_Neg/global/Differential_expression/Differential_expression_results.csv",
#                   "DIA: flow" = "analysis/combined24-27/DIA_2batches_noOutliers_noMSC/Sort Type_Flow/CD14_Pos_vs_Neg/global/Differential_expression/Differential_expression_results.csv",
#                   "TMT: bead" = "analysis/combined24-27/TMT_noOutliers_noMSC/Sort Type_Bead/CD14_Pos_vs_Neg/global/Differential_expression/Differential_expression_results.csv",
#                   "TMT: flow" = "analysis/combined24-27/TMT_noOutliers_noMSC/Sort Type_Flow/CD14_Pos_vs_Neg/global/Differential_expression/Differential_expression_results.csv")
sig.paths <- list("DIA: bead" = "analysis/combined24-27/using_cellType-sortType-patient_factors/DIA_2batches_noOutliers_noMSC/Sort Type_Bead/CD14_Pos_vs_Neg/global/Differential_expression/Differential_expression_results.csv",
                  "DIA: flow" = "analysis/combined24-27/using_cellType-sortType-patient_factors/DIA_2batches_noOutliers_noMSC/Sort Type_Flow/CD14_Pos_vs_Neg/global/Differential_expression/Differential_expression_results.csv",
                  "TMT: bead" = "analysis/combined24-27/using_cellType-sortType-patient_factors/TMT_noOutliers_noMSC/Sort Type_Bead/CD14_Pos_vs_Neg/global/Differential_expression/Differential_expression_results.csv",
                  "TMT: flow" = "analysis/combined24-27/using_cellType-sortType-patient_factors/TMT_noOutliers_noMSC/Sort Type_Flow/CD14_Pos_vs_Neg/global/Differential_expression/Differential_expression_results.csv")

# import signatures and filter
sigs <- list()
venn.list <- list()
setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/")
for (i in names(sig.paths)) {
  if (file.exists(sig.paths[[i]])) { # no TMT flow for some reason
    sigs[[i]] <- read.csv(sig.paths[[i]]) 
    sigs[[i]] <- na.omit(sigs[[i]][sigs[[i]]$adj.P.Val <= 0.05,c("Gene","Log2FC")])
    venn.list[[i]] <- sigs[[i]]$Gene
  }
}
ggvenn::ggvenn(venn.list, show_percentage=FALSE, set_name_size=5, text_size=5)
ggplot2::ggsave("mono_vs_prog_signature_vennDiagram_diffMethods.pdf",width=7, height=7)

# run correlations between signatures
sig.df <- data.table::rbindlist(sigs, use.names=TRUE, idcol="Signature")
sig.df <- reshape2::dcast(sig.df, Gene ~ Signature, value.var="Log2FC")
corr <- data.frame()
for (i in names(sigs)) {
  otherSigs <- names(sigs)[names(sigs) != i]
  temp.input <- sig.df[,c("Gene",i,otherSigs)]
  temp.corr <- DMEA::rank_corr(temp.input,plots=FALSE)$result
  temp.corr[nrow(temp.corr)+1,] <- c(i,1,rep(NA, ncol(temp.corr)-2))
  temp.corr$Signature <- i
  corr <- rbind(corr, temp.corr)
}
corr$Pearson.est <- as.numeric(corr$Pearson.est)
corr$Pearson.p <- as.numeric(corr$Pearson.p)
corr$Pearson.q <- as.numeric(corr$Pearson.q)
write.csv(corr, "mono_vs_prog_signature_correlations_withSelfCorr_diffMethods.csv", row.names=FALSE)
ggplot(corr, aes(x=Drug, y=Signature, fill=Pearson.est)) + geom_tile() + 
  scale_fill_gradient2(limits=c(-1,1), low="blue", mid="grey", high="red")+labs(fill="Pearson r")+
  theme_minimal(base_size=16) + theme(axis.text=element_text(vjust=1, hjust=1, angle=45), axis.title=element_blank())
ggsave("mono_vs_prog_signature_correlations_withSelfCorr_diffMethods.pdf", width=4, height=2.5)

# import signatures and filter
sigs <- list()
venn.list <- list()
setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/")
for (i in names(sig.paths)) {
  if (file.exists(sig.paths[[i]])) {
    sigs[[i]] <- read.csv(sig.paths[[i]])
    sigs[[i]] <- na.omit(sigs[[i]])
    venn.list[[i]] <- sigs[[i]]$Gene 
  }
}
ggvenn::ggvenn(venn.list, show_percentage=FALSE, set_name_size=5, text_size=5)
ggsave("mono_vs_prog_signature_vennDiagram_diffMethods_unfiltered.pdf",width=7, height=7)

# run correlations between signatures
sig.df <- data.table::rbindlist(sigs, use.names=TRUE, idcol="Signature")
sig.df <- reshape2::dcast(sig.df, Gene ~ Signature, value.var="Log2FC")
corr <- data.frame()
for (i in names(sigs)) {
  otherSigs <- names(sigs)[names(sigs) != i]
  temp.input <- sig.df[,c("Gene",i,otherSigs)]
  temp.corr <- DMEA::rank_corr(temp.input,plots=FALSE)$result
  temp.corr[nrow(temp.corr)+1,] <- c(i,1,rep(NA, ncol(temp.corr)-2))
  temp.corr$Signature <- i
  corr <- rbind(corr, temp.corr)
}
corr$Pearson.est <- as.numeric(corr$Pearson.est)
corr$Pearson.p <- as.numeric(corr$Pearson.p)
corr$Pearson.q <- as.numeric(corr$Pearson.q)
write.csv(corr, "mono_vs_prog_signature_correlations_withSelfCorr_diffMethods_unfiltered.csv", row.names=FALSE)
ggplot(corr, aes(x=Drug, y=Signature, fill=Pearson.est)) + geom_tile() + 
  scale_fill_gradient2(limits=c(-1,1), low="blue", mid="grey", high="red")+labs(fill="Pearson r")+
  theme_minimal(base_size=16) + theme(axis.text=element_text(vjust=1, hjust=1, angle=45), axis.title=element_blank())
ggsave("mono_vs_prog_signature_correlations_withSelfCorr_diffMethods_unfiltered.pdf", width=4, height=2.5)

