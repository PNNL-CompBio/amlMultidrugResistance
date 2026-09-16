# running mono vs. prog sorted sig on exp 28

# load crosstab
exp28 <- read.table("~/Downloads/exp28_log2_sample_centered_crosstab_complete_features 1.txt", sep="\t")
exp28 <- as.data.frame(t(exp28))
exp28$sampleID <- rownames(exp28)
exp28 <- exp28[,c("sampleID",colnames(exp28)[colnames(exp28)!="sampleID"])]

# just check how many missing
exp28notfull <- read.table("~/Downloads/exp28_log2_sample_centered_crosstab 1.txt", sep="\t")
exp28notfull <- as.data.frame(t(exp28notfull))
exp28notfull$sampleID <- rownames(exp28notfull)
exp28notfull <- exp28notfull[,c("sampleID",colnames(exp28notfull)[colnames(exp28notfull)!="sampleID"])]

# load sigs
sigs2 <- readRDS(file.path("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/",
                           "Monocyte_vs_progenitor_signatures_beadOnly_2025-07-14/topMixVenSensSigs_2025-07-15.rds"))
sigs <- list("Sorted" = read.csv(file.path("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/",
                                           "analysis/combined24-27/using_cellType-sortType-patient_factors/DIA_2batches_noOutliers_noMSC_Bead/no_filter/CD14_Pos_vs_Neg/global/Differential_expression/Differential_expression_results_max_5_percent_FDR.csv")),
                  "Sorted_25" = sigs2$`Sorted: 25 proteins`)

# perform WV
wv.results <- list()
wv.dfs <- list()
unused.dfs <- list()
for (i in names(sigs)) {
  wv.results[[i]] <- DMEA::WV(exp28,sigs[[i]])
  wv.dfs[[i]] <- wv.results[[i]]$scores
  unused.dfs[[i]] <- wv.results[[i]]$unused.weights
}
wv.df <- data.table::rbindlist(wv.dfs, use.names=TRUE, idcol="Signature")
write.csv(wv.df, "Exp28_scoredUsingSortedSigs.csv", row.names=FALSE)

unused <- data.table::rbindlist(unused.dfs, use.names=TRUE, idcol="Signature", fill=TRUE)
write.csv(unused, "Exp28_unusedSortedProteins.csv", row.names=FALSE)

sig.df <- data.table::rbindlist(sigs, use.names=TRUE, idcol="Signature", fill=TRUE)
write.csv(sig.df, "sortedProteinSignatures.csv", row.names=TRUE)

sorted.full <- sigs$Sorted
sorted.unused <- unused[unused$Signature=="Sorted",]
sorted.used <- sorted.full[!(sorted.full$Gene %in% sorted.unused$Gene),]

sorted.25 <- sigs$Sorted_25
sorted.25.unused <- unused[unused$Signature=="Sorted_25",]
sorted.25.used <- sorted.25[!(sorted.25$Gene %in% sorted.25.unused$Gene),]

sorted.used.all <- list("Sorted" = sorted.used, "Sorted_25" = sorted.25.used)
used <- data.table::rbindlist(sorted.used.all, use.names=TRUE, idcol="Signature", fill=TRUE)
write.csv(used, "Exp28_usedSortedProteins.csv", row.names=FALSE)

# check for completely missing
missing <- sorted.full[!(sorted.full$Gene %in% colnames(exp28notfull)[2:ncol(exp28notfull)]),] # 88
found <- sorted.full[(sorted.full$Gene %in% colnames(exp28notfull)[2:ncol(exp28notfull)]),] # 2144
missing25 <- sorted.25[!(sorted.25$Gene %in% colnames(exp28notfull)[2:ncol(exp28notfull)]),] # 0
found25 <- sorted.25[(sorted.25$Gene %in% colnames(exp28notfull)[2:ncol(exp28notfull)]),] # 25
