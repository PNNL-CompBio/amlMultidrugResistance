# does cell type correlate with drug sensitivity?
library(synapser)
library(DMEA)
library(tidyr)
library(tibble)
library(reshape2)
synapser::synLogin()

# load deconvolution results
camilo <- synapser::synTableQuery("select * from syn51788432 WHERE algorithm = 'xcell' AND matrix = 'AML' AND dataType = 'mrna'")$asDataFrame()
camilo <- reshape2::dcast(camilo, sample ~ `Cell type`, value.var = "cellPop")
wv <- read.table("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/decomprolute/compareSigs/aml/aml-tumor-prot-raw-wv-AML_vanGalen_100.tsv",fill=TRUE,sep = '\t')
colnames(wv) <- wv[1,]
wv <- wv[2:nrow(wv),]

# load Beat AML drug AUC
drug.BeatAML <- read.csv(synapser::synGet("syn51674470")$path)
drug.BeatAML$X <- NULL
colnames(drug.BeatAML)[1] <- "sample"
drug.BeatAML.wide <- reshape2::dcast(drug.BeatAML, sample ~ inhibitor, value.var = "auc")
colnames(drug.BeatAML)[1] <- "sample" # to match sample column name in camilo
ven.AUC <- drug.BeatAML[drug.BeatAML$inhibitor == "Venetoclax",]
AV.AUC <- drug.BeatAML[drug.BeatAML$inhibitor == "Azacytidine - Venetoclax",]
aza.AUC <- drug.BeatAML[drug.BeatAML$inhibitor == "Azacytidine",]

# for venetoclax or Aza+Ven, does drug AUC correlate with any cell type score?
camilo.ven.AUC <- merge(ven.AUC[,c("sample", "auc")], camilo, by="sample")
camilo.ven.corr <- DMEA::rank_corr(camilo.ven.AUC) # monocytes are top correlated cell type with Pearson r = 0.7470501, p = 4.751839e-21, q = 6.652574e-20, N = 111; similar to sorted sig r = 0.746, p = 7.99E-24 with 25 proteins
write.csv(camilo.ven.corr$result, "BeatAML_xcell_deconvolution_correlation_with_Ven_AUC.csv", row.names = FALSE)
ggsave("BeatAML_xcell_deconvolution_correlation_with_Ven_AUC_scatterPlots.pdf", camilo.ven.corr$scatter.plots)

camilo.AV.AUC <- merge(AV.AUC[,c("sample", "auc")], camilo, by="sample")
camilo.AV.corr <- DMEA::rank_corr(camilo.AV.AUC) # monocytes are top correlated cell type with Pearson r = 0.81483554, p = 1.214319e-05, q = 0.0001700046, N = 20; similar to sorted sig r = 0.804, p = 1.59E-5 with 25 proteins
write.csv(camilo.AV.corr$result, "BeatAML_xcell_deconvolution_correlation_with_AzaVen_AUC.csv", row.names = FALSE)
ggsave("BeatAML_xcell_deconvolution_correlation_with_AzaVen_AUC_scatterPlots.pdf", camilo.AV.corr$scatter.plots)


# other questions
# does ratio of mono/prog in Beat AML correlate with sorted cell CD14/CD34 signature?
globalFileDIA <- synapser::synGet("syn59429685") # filtered for max FDR of 0.05
global.sig.DIA <- read.csv(globalFileDIA$path) # 1842 proteins
global.sig.25 <- global.sig.DIA %>% slice_max(abs(Log2FC), n=25)

# load Beat AML global data
meta.BeatAML <- read.table(synapser::synGet("syn25807733")$path, 
                           sep = "\t", header = TRUE)
global.BeatAML <- read.table(synapser::synGet("syn25714248")$path,
                             sep = "\t", header = TRUE)

# change global.BeatAML column names from SampleID.abbrev to 
# Barcode.ID to match drug.BeatAML
global.ids <- names(global.BeatAML)

# remove X and any 0's from start of each column name and then
# replace SampleID.abbrev with Barcode.ID to match drug.BeatAML
for(i in seq_len(length(global.ids))){
  global.ids[i] <- substr(global.ids[i], 2, nchar(global.ids[i]))
  
  if(substring(global.ids[i], 1, 1) == 0){
    global.ids[i] <- substr(global.ids[i], 2, nchar(global.ids[i]))
  }
  
  if(global.ids[i] %in% meta.BeatAML$SampleID.abbrev){
    global.ids[i] <- meta.BeatAML[meta.BeatAML$SampleID.abbrev == global.ids[i], ]$Barcode.ID
  }
}

# replace global.BeatAML column names 
names(global.BeatAML) <- global.ids

# transpose global.BeatAML so that first column is Barcode.ID and 
# rest of columns are gene symbols
global.BeatAML <- as.data.frame(t(global.BeatAML))

# make first column sample
global.BeatAML[, "sample"] <- rownames(global.BeatAML)
global.BeatAML <- 
  global.BeatAML[ , c("sample", 
                      names(global.BeatAML[ , 1:(ncol(global.BeatAML)-1)]))]

# filter for genes with 100% coverage in Beat AML database
global100 <- global.BeatAML[,colSums(is.na(global.BeatAML)) == 0]

global100.wv <- DMEA::WV(global100, global.sig.25[,c("Gene", "Log2FC")]) # CPM was unused gene
write.csv(global100.wv$scores, "WV_BeatAML_top25DIAGlobal_by_absLog2FC.csv", row.names = FALSE)
global100.wv <- read.csv("WV_BeatAML_top25DIAGlobal_by_absLog2FC.csv")

# correlation between rna deconv of bulk samples and WV scores of sorted samples
all.scores <- merge(global100.wv$scores, camilo, by = "sample")
cor.bulk.sorted.scores <- DMEA::rank_corr(all.scores)
write.csv(cor.bulk.sorted.scores$result, "BeatAML_xcell_deconvolution_correlation_with_WV_sorted_CD14vsCD34_top_25_absLog2FC_DIA.csv", row.names = FALSE)
ggsave("BeatAML_xcell_deconvolution_correlation_with_WV_sorted_CD14vsCD34_top_25_absLog2FC_DIA_scatterPlots.pdf", cor.bulk.sorted.scores$scatter.plots)

# strongest correlation is again monocyte: Pearson r = 0.8603639, p = 8.192483e-48, q = 1.146948e-46
# does ratio of mono/prog in Beat AML better predict drug sensitivity than sorted cell CD14/CD34 signature? (similar) what about a combination of the two scores?
# try mono score of xcell rna van galen + wv sorted prot cd14 vs cd34
all.scores$sum_monocyte_WV <- all.scores$WV + all.scores$Monocyte
combo.score <- all.scores[,c("sample", "sum_monocyte_WV")]

# for venetoclax or Aza+Ven, does drug AUC correlate with any cell type score?
combo.drug.AUC <- merge(combo.score, drug.BeatAML.wide, by="sample")
combo.corr <- DMEA::rank_corr(combo.drug.AUC) # Ven, Aza+Ven are correlated
# Ven result is a bit worse Pearson r = 0.694086775, p = 3.011629e-17, q = 9.336048e-15 and rank 8 with N = 111
# Aza - Ven is 2nd top correlated drug treatment but still a bit worse with Pearson r = 0.782821788, p = 4.498523e-05, q = 7.339695e-04, N = 20
write.csv(combo.corr$result, "BeatAML_xcell_deconvolution_plus_sorted_WV_correlation_with_drug_AUC.csv", row.names = FALSE)
ggsave("BeatAML_xcell_deconvolution_plus_sorted_WV_correlation_with_drug_AUC_scatterPlots.pdf", combo.corr$scatter.plots)

# what if we just use the 10 patients studied for sorted proteomics?
sorted.patients <- c("16-01184", "17-01060","18-00103","18-00105","19-00074",
                     "21-00432","21-00839","22-00117","22-00251","22-00571")

# for venetoclax or Aza+Ven, does drug AUC correlate with any cell type score?
camilo.10 <- camilo[camilo$sample %in% sorted.patients,] # why only 3 patients?
camilo.ven.AUC.10 <- camilo.ven.AUC[camilo.ven.AUC$sample %in% sorted.patients,] # why only 2 patients?

camilo.AV.AUC.10 <- camilo.AV.AUC[camilo.AV.AUC$sample %in% sorted.patients,] # why only 3 patients?
camilo.AV.corr.10 <- DMEA::rank_corr(camilo.AV.AUC.10) # monocytes are top correlated cell type with Pearson r = 0.81483554, p = 1.214319e-05, q = 0.0001700046, N = 20; similar to sorted sig r = 0.804, p = 1.59E-5 with 25 proteins
write.csv(camilo.AV.corr.10$result, "BeatAML_xcell_deconvolution_correlation_with_AzaVen_AUC_3_sorted_patients.csv", row.names = FALSE)
ggsave("BeatAML_xcell_deconvolution_correlation_with_AzaVen_AUC_3_sorted_patients_scatterPlots.pdf", camilo.AV.corr.10$scatter.plots)
# Monocyte was 3rd ranked cell type with Pearson r = 0.992, p = 0.0813, q = 0.219

# correlation between rna deconv of bulk samples and WV scores of sorted samples
all.scores.10 <- all.scores[all.scores$sample %in% sorted.patients,]
cor.bulk.sorted.scores.10 <- DMEA::rank_corr(all.scores.10)
write.csv(cor.bulk.sorted.scores.10$result, "BeatAML_xcell_deconvolution_correlation_with_WV_sorted_CD14vsCD34_top_25_absLog2FC_DIA_3_sorted_patients.csv", row.names = FALSE)
ggsave("BeatAML_xcell_deconvolution_correlation_with_WV_sorted_CD14vsCD34_top_25_absLog2FC_DIA_3_sorted_patients_scatterPlots.pdf", cor.bulk.sorted.scores.10$scatter.plots)
# monocyte is 6th top cell type: Pearson r = 0.9361791, p = 0.22867240, q = 0.3389524
# does ratio of mono/prog in Beat AML better predict drug sensitivity than sorted cell CD14/CD34 signature? (similar) what about a combination of the two scores?
# try mono score of xcell rna van galen + wv sorted prot cd14 vs cd34
combo.drug.AUC.10 <- combo.drug.AUC[combo.drug.AUC$sample %in% sorted.patients,] # 3 patients
combo.corr.10 <- DMEA::rank_corr(combo.drug.AUC.10)
# ven didn't have 3+ points to be evaluated
# Aza - Ven is 45th top correlated drug treatment with Pearson r = 0.8840447, p = 0.3096201, q = 0.6890094, N = 3
write.csv(combo.corr.10$result, "BeatAML_xcell_deconvolution_plus_sorted_WV_correlation_with_drug_AUC_3_sorted_patients.csv", row.names = FALSE)
ggsave("BeatAML_xcell_deconvolution_plus_sorted_WV_correlation_with_drug_AUC_scatterPlots_3_sorted_patients.pdf", combo.corr.10$scatter.plots)

# try just WV 10:
WV.drug.AUC <- merge(all.scores[,c("sample", "WV")], drug.BeatAML.wide) # 159 patients
WV.drug.AUC.10 <- WV.drug.AUC[WV.drug.AUC$sample %in% sorted.patients,] # 3 patients
WV.corr.10 <- DMEA::rank_corr(WV.drug.AUC.10)
# ven didn't have 3+ points to be evaluated
# Aza - Ven is 45th top correlated drug treatment with Pearson r = 0.8838070, p = 0.3099436, q = 0.6896608, N = 3
write.csv(WV.corr.10$result, "BeatAML_sorted_WV_correlation_with_drug_AUC_3_sorted_patients.csv", row.names = FALSE)
ggsave("BeatAML_sorted_WV_correlation_with_drug_AUC_scatterPlots_3_sorted_patients.pdf", WV.corr.10$scatter.plots)

WV.corr <- DMEA::rank_corr(WV.drug.AUC) # Ven, Aza+Ven are correlated
# ven is 8th top correlated drug treatment with Pearson r = 0.69385670, p = 3.115443e-17, q = 9.657874e-15, N = 111
# Aza - Ven is 2nd top correlated drug treatment with Pearson r = 0.7826665695, p= 4.524763e-05, q = 7.382508e-04, N = 20
write.csv(WV.corr$result, "BeatAML_sorted_WV_correlation_with_drug_AUC.csv", row.names = FALSE)
ggsave("BeatAML_sorted_WV_correlation_with_drug_AUC_scatterPlots.pdf", WV.corr$scatter.plots)

# next steps:
# try with CD14 vs. all score instead of CD14 vs. CD34
sorted.AML <- read.table("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/decomprolute/compareSigs/aml/aml-tumor-prot-raw-wv-AML_sorted_100.tsv",fill=TRUE,sep = '\t')
colnames(sorted.AML) <- sorted.AML[1,]
sorted.AML <- sorted.AML[2:nrow(sorted.AML),]
rownames(sorted.AML) <- c("Monocyte", "Progenitor", "MSC") # remove -like from end of cell types

mono.sorted <- t(sorted.AML[1,])
# need to map V1, V2, etc. to sample IDs

###
# redo sig comparison with cd14 vs other sig
cd14.sig <- read.csv(synapser::synGet("syn59424916")$path) # DIA
cd14.sig.25 <- cd14.sig %>% slice_max(abs(Log2FC), n=25)

# Beat AML WV
global100.wv <- DMEA::WV(global100, cd14.sig.25[,c("Gene", "Log2FC")]) # CPM was unused gene
write.csv(global100.wv$scores, "WV_BeatAML_top25DIAGlobalCD14_by_absLog2FC.csv", row.names = FALSE)
global100.wv.cd14 <- read.csv("WV_BeatAML_top25DIAGlobalCD14_by_absLog2FC.csv")

# correlation between rna deconv of bulk samples and WV scores of sorted samples
all.scores <- merge(global100.wv$scores, camilo, by = "sample")
cor.bulk.sorted.scores <- DMEA::rank_corr(all.scores)
write.csv(cor.bulk.sorted.scores$result, "BeatAML_xcell_deconvolution_correlation_with_WV_sorted_CD14vsOther_top_25_absLog2FC_DIA.csv", row.names = FALSE)
ggsave("BeatAML_xcell_deconvolution_correlation_with_WV_sorted_CD14vsOther_top_25_absLog2FC_DIA_scatterPlots.pdf", cor.bulk.sorted.scores$scatter.plots)

# strongest correlation is again monocyte: Pearson r = 0.8603639, p = 8.192483e-48, q = 1.146948e-46
# does ratio of mono/prog in Beat AML better predict drug sensitivity than sorted cell CD14/CD34 signature? (similar) what about a combination of the two scores?
# try mono score of xcell rna van galen + wv sorted prot cd14 vs cd34
all.scores$sum_monocyte_WV <- all.scores$WV + all.scores$Monocyte
combo.score <- all.scores[,c("sample", "sum_monocyte_WV")]

# for venetoclax or Aza+Ven, does drug AUC correlate with any cell type score?
combo.drug.AUC <- merge(combo.score, drug.BeatAML.wide, by="sample")
combo.corr <- DMEA::rank_corr(combo.drug.AUC) # Ven, Aza+Ven are correlated
# Ven result is a bit worse Pearson r = 0.694086775, p = 3.011629e-17, q = 9.336048e-15 and rank 8 with N = 111
# Aza - Ven is 2nd top correlated drug treatment but still a bit worse with Pearson r = 0.782821788, p = 4.498523e-05, q = 7.339695e-04, N = 20
write.csv(combo.corr$result, "BeatAML_xcell_deconvolution_plus_sorted_WV_CD14vsOther_correlation_with_drug_AUC.csv", row.names = FALSE)
ggsave("BeatAML_xcell_deconvolution_plus_sorted_WV_CD14vsOther_correlation_with_drug_AUC_scatterPlots.pdf", combo.corr$scatter.plots)

# what about just CD14 vs. other
combo.drug.AUC <- merge(global100.wv$scores, drug.BeatAML.wide, by="sample")
combo.corr <- DMEA::rank_corr(combo.drug.AUC) # Ven, Aza+Ven are correlated
# Ven result is a bit worse Pearson r = 0.694086775, p = 3.011629e-17, q = 9.336048e-15 and rank 8 with N = 111
# Aza - Ven is 2nd top correlated drug treatment but still a bit worse with Pearson r = 0.782821788, p = 4.498523e-05, q = 7.339695e-04, N = 20
write.csv(combo.corr$result, "BeatAML_sorted_WV_CD14vsOther_correlation_with_drug_AUC.csv", row.names = FALSE)
ggsave("BeatAML_sorted_WV_CD14vsOther_correlation_with_drug_AUC_scatterPlots.pdf", combo.corr$scatter.plots)

# is this better than GBT
gbt.prot.ven.pearson <- c("0.972325815","0.55449849","-0.28147852",
                          "-0.14253104","0.328844856","0.711020363",
                          "-0.457582208","-0.099684927","0.607069431",
                          "0.494702538","0.512396604","0.968494829",
                          "0.930644642","-0.925853916")
gbt.rna.ven.pearson <- c("0.905819016","-0.194143678","0.658213465",
                         "-0.343268218","0.241723849","0.648055495",
                         "0.930182236","0.806121592","0.95229431",
                         "0.780614575","0.746470005","0.627643838",
                         "0.998350314","0.969552207","-0.190122135")
t.test(gbt.prot.ven.pearson, mu=0.7, alternative = "less")
t.test(gbt.rna.ven.pearson, mu=0.746, alternative = "less")

# try combining scores - scaling wv this time
all.scores <- merge(global100.wv, global100.wv.cd14, by = "sample", suffixes = c("_CD14vsCD34", "_CD14vsOther"))
all.scores <- merge(all.scores, camilo, by = "sample")
all.scores$wv_scaled <- scales::rescale(all.scores$WV_CD14vsCD34)
all.scores$wv_cd14_scaled <- scales::rescale(all.scores$WV_CD14vsOther)
all.scores$combo <- all.scores$wv_scaled + all.scores$Monocyte
all.scores$combo_cd14 <- all.scores$wv_cd14_scaled + all.scores$Monocyte
all.scores[,c("mean_combo", "mean_combo_cd14")] <- NA
for (i in 1:nrow(all.scores)) {
  all.scores$mean_combo[i] <- mean(c(all.scores$wv_scaled[i], all.scores$Monocyte[i]))
  all.scores$mean_combo_cd14[i] <- mean(c(all.scores$wv_cd14_scaled[i], all.scores$Monocyte[i])) 
}

# sum
# first try combo score using CD14 vs CD34 signature
combo.drug.AUC <- merge(all.scores[,c("sample", "combo")], drug.BeatAML.wide, by="sample")
combo.corr <- DMEA::rank_corr(combo.drug.AUC) 
write.csv(combo.corr$result, "BeatAML_xcell_deconvolution_plus_scaled_sorted_WV_CD14vsCD34_correlation_with_drug_AUC.csv", row.names = FALSE)
ggsave("BeatAML_xcell_deconvolution_plus_scaled_sorted_WV_CD14vsCD34_correlation_with_drug_AUC_scatterPlots.pdf", combo.corr$scatter.plots)

# next try combo score using CD14 vs other signature
combo.drug.AUC <- merge(all.scores[,c("sample", "combo_cd14")], drug.BeatAML.wide, by="sample")
combo.corr <- DMEA::rank_corr(combo.drug.AUC) 
write.csv(combo.corr$result, "BeatAML_xcell_deconvolution_plus_scaled_sorted_WV_CD14vsOther_correlation_with_drug_AUC.csv", row.names = FALSE)
ggsave("BeatAML_xcell_deconvolution_plus_scaled_sorted_WV_CD14vsOther_correlation_with_drug_AUC_scatterPlots.pdf", combo.corr$scatter.plots)

# mean
# first try combo score using CD14 vs CD34 signature
combo.drug.AUC <- merge(all.scores[,c("sample", "mean_combo")], drug.BeatAML.wide, by="sample")
combo.corr <- DMEA::rank_corr(combo.drug.AUC) 
write.csv(combo.corr$result, "BeatAML_mean_xcell_deconvolution_with_scaled_sorted_WV_CD14vsCD34_correlation_with_drug_AUC.csv", row.names = FALSE)
ggsave("BeatAML_mean_xcell_deconvolution_with_scaled_sorted_WV_CD14vsCD34_correlation_with_drug_AUC_scatterPlots.pdf", combo.corr$scatter.plots)

# next try combo score using CD14 vs other signature
combo.drug.AUC <- merge(all.scores[,c("sample", "mean_combo_cd14")], drug.BeatAML.wide, by="sample")
combo.corr <- DMEA::rank_corr(combo.drug.AUC) 
write.csv(combo.corr$result, "BeatAML_mean_xcell_deconvolution_with_scaled_sorted_WV_CD14vsOther_correlation_with_drug_AUC.csv", row.names = FALSE)
ggsave("BeatAML_mean_xcell_deconvolution_with_scaled_sorted_WV_CD14vsOther_correlation_with_drug_AUC_scatterPlots.pdf", combo.corr$scatter.plots)

# try to optimize signature (# of proteins, which proteins, Log2FC or adj.P.Value or both?)

# try with same deconvolution algorithm
# would like to do fair comparison by adding DMEA to decomprolute

# does measured cell fraction correlate with drug response?
base.path <- "~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/data/clinical_metadata/"
setwd(base.path)
patient.key <- readxl::read_xlsx("lab_dbgap_key.xlsx")
patient.meta <- readxl::read_xlsx("Table_S1.xlsx",sheet=1)
patient.meta <- as.data.frame(patient.meta)
#rownames(patient.meta) <- patient.meta[,1]

# does metadata correlate Ven response?
key.cols <- colnames(patient.meta)[colnames(patient.meta) %in% colnames(patient.key)]
num.cols <- colnames(dplyr::select_if(patient.meta, is.numeric))
num.meta <- dplyr::distinct(patient.meta[,c(key.cols, num.cols)])
num.meta <- merge(patient.key, num.meta, by = key.cols)
num.meta <- num.meta[,c("labId", num.cols[2:length(num.cols)])] # leave out dbgap_subject_id from numeric columns
colnames(num.meta)[1] <- "sample"

meta.ven <- merge(ven.AUC[,c(1,3)], num.meta, by="sample")
meta.ven.corr <- DMEA::rank_corr(meta.ven, variable="metadata", value="value", FDR = 1, xlab = "Venetoclax AUC") # set FDR 1 to get all plots
write.csv(meta.ven.corr$result, "BeatAML_Ven_AUC_correlation_with_metadata.csv", row.names = FALSE)
ggsave("BeatAML_Ven_AUC_correlation_with_metadata.pdf", meta.ven.corr$scatter.plots)

meta.av <- merge(AV.AUC[,c(1,3)], num.meta, by="sample")
meta.av.corr <- DMEA::rank_corr(meta.av, variable="metadata", value="value", FDR=1, xlab = "Aza + Ven AUC") # set FDR 1 to get all plots
write.csv(meta.av.corr$result, "BeatAML_AzaVen_AUC_correlation_with_metadata.csv", row.names = FALSE)
ggsave("BeatAML_AzaVen_AUC_correlation_with_metadata.pdf", meta.av.corr$scatter.plots)

meta.aza <- merge(aza.AUC[,c(1,3)], num.meta, by="sample")
meta.aza.corr <- DMEA::rank_corr(meta.aza, variable="metadata", value="value", FDR = 1, xlab = "Azacitidine AUC") # set FDR 1 to get all plots
write.csv(meta.aza.corr$result, "BeatAML_Aza_AUC_correlation_with_metadata.csv", row.names = FALSE)
ggsave("BeatAML_Aza_AUC_correlation_with_metadata.pdf", meta.aza.corr$scatter.plots)

# are NPM1-mutant AMLs more resistant?
npm1.mut <- na.omit(unique(patient.meta[patient.meta$NPM1 == "negative",]$dbgap_subject_id)) # 610
npm1.wt <- na.omit(unique(patient.meta[patient.meta$NPM1 == "positive",]$dbgap_subject_id)) # 202
overlap <- npm1.mut[npm1.mut %in% npm1.wt] # 9
npm1.mut.noOverlap <- npm1.mut[!(npm1.mut %in% overlap)] # 601
npm1.wt.noOverlap <- npm1.wt[!(npm1.wt %in% overlap)] # 193
id.npm1.mut <- unique(patient.key[patient.key$dbgap_subject_id %in% npm1.mut.noOverlap,]$labId) # 690
id.npm1.wt <- unique(patient.key[patient.key$dbgap_subject_id %in% npm1.wt.noOverlap,]$labId) # 227
npm1.ven.test <- stats::t.test(ven.AUC[ven.AUC$sample %in% id.npm1.mut,]$auc,
                               ven.AUC[ven.AUC$sample %in% id.npm1.wt,]$auc) # p = 0.91 from 68 vs. 39 samples
npm1.aza.test <- stats::t.test(aza.AUC[aza.AUC$sample %in% id.npm1.mut,]$auc,
                               aza.AUC[aza.AUC$sample %in% id.npm1.wt,]$auc) # p = 0.41 from 75 vs. 41 samples
npm1.av.test <- stats::t.test(AV.AUC[AV.AUC$sample %in% id.npm1.mut,]$auc,
                               AV.AUC[AV.AUC$sample %in% id.npm1.wt,]$auc) # p = 0.31 from 11 vs. 5 samples
