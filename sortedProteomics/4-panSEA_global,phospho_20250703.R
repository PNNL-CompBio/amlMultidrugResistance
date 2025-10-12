# Differential expression & enrichment analyses: global & phospho
# PTRC2: exp 24 TMT
# Author: Belinda B. Garana
# Created: 2023-12-06
# Last edit: 2024-02-09
remove(list=ls())
library(readxl); library(panSEA); library(synapser)
library(stringr); library(tidyr); library(dplyr); library(Biobase)

setwd("~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/")
#source("panSEA_helper_20240508_updated20241015.R")
#source("~/OneDrive - PNNL/Documents/GitHub/Exp21_NRAS_ASO_treated_patients/proteomics/panSEA_helper_20240508_updated20250131.R")
source("https://raw.githubusercontent.com/PNNL-CompBio/Exp21_NRAS_ASO_treated_patients/refs/heads/main/proteomics/panSEA_helper_20240508_updated20250131.R")

# overview
# 1. Import metadata & crosstabs
# 2. Import BeatAML data formatted for DMEA
# 3. Run panSEA across omics for each contrast

#### 1. Import metadata & crosstabs ####
### TMT
setwd("data")
meta.df <- readxl::read_excel("Exp24metadataTable_TMT.xlsx") 
meta.df$MeasurementName <- as.character(rownames(meta.df))
rownames(meta.df) <- paste0("X", rownames(meta.df))
meta.df$row.name <- rownames(meta.df)

# add other drug info & make sure sensitivity is correctly labeled
sens.info <- read.csv("Exp24_drug_sensitivity_20240209.csv") # also syn53627410
meta.df <- merge(meta.df, sens.info)
meta.df$X <- NULL
# add 'X' to MeasurementName to match colnames of data
meta.df$id <- paste0('X', meta.df$MeasurementName)

# add other metadata for contrasts
meta.df$CD14 <- "Neg"
#meta.df[grepl("CD14+", meta.df$SampleType),]$CD14_Pos <- "Pos"
meta.df[meta.df$SampleType == "CD14+" | meta.df$SampleType == "CD14+ Flow", ]$CD14 <- "Pos"

meta.df$CD34 <- "Neg"
#meta.df[grepl("CD34+", meta.df$SampleType),]$CD34_Pos <- "Pos"
meta.df[meta.df$SampleType == "CD34+" | meta.df$SampleType == "CD34+ Flow", ]$CD34 <- "Pos"

meta.df$MSC <- "Non_MSC"
meta.df[meta.df$SampleType == "MSC Flow",]$MSC <- "MSC"

meta.df$'Sort Type' <- "Bead"
meta.df[grepl("Flow", meta.df$SampleType),]$'Sort Type' <- "Flow"
rownames(meta.df) <- meta.df$row.name
meta.df$`Sample Type` <- meta.df$SampleType
global.df <- read.table(
  "global_data/Exp24_crosstab_global_gene_corrected.txt", 
  sep = "\t")
phospho.df <- read.table(
  "phospho_data/Exp24_crosstab_phospho_SiteID_corrected.txt", 
  sep = "\t")
phospho.pep.TMT <- read.table(
  "phospho_data/Exp24_crosstab_phospho_peptide_corrected_labelsMatchDIAPhospho.txt", 
  sep = "\t") # 831 peptides

# fix column names to match DIA
sample.names <- meta.df[colnames(phospho.pep.TMT)[1:(ncol(phospho.pep.TMT)-1)],]$DIA_id
colnames(phospho.pep.TMT)[1:(ncol(phospho.pep.TMT)-1)] <- sample.names

# add column for feature names and later make it the first column
global.df$Gene <- rownames(global.df)
phospho.df$SUB_SITE <- rownames(phospho.df)
sample.names <- meta.df[colnames(global.df)[1:(ncol(global.df)-1)],]$DIA_id
colnames(global.df)[1:(ncol(global.df)-1)] <- sample.names
sample.names <- meta.df[colnames(phospho.df)[1:(ncol(phospho.df)-1)],]$DIA_id
colnames(phospho.df)[1:(ncol(phospho.df)-1)] <- sample.names
tmt <- list("meta" = meta.df,
            "global" = global.df, # 5169 gene symbols
            "phospho" = phospho.df) # 794 phospho-sites
rownames(tmt$meta) <- tmt$meta$DIA_id 

# determine outliers based on CD14 expression should be higher than CD34 in CD14+ samples and vice-versa
val.df <- reshape2::melt(tmt$global[c("CD14","CD34"),])
val.df$`Cell Type` <- "CD14+"
val.df[grepl("cd34", val.df$variable, ignore.case=TRUE),]$`Cell Type` <- "CD34+"
val.df[grepl("msc", val.df$variable, ignore.case=TRUE),]$`Cell Type` <- "MSC"
val.df$`Sorting Method` <- "Bead"
val.df[grepl("flow",val.df$variable, ignore.case=TRUE) | grepl("_f_", val.df$variable),]$`Sorting Method` <- "Flow"
val.df$Patient <- stringr::str_split_i(val.df$variable, "_", 1)
ggplot2::ggplot(val.df, aes(x=`Cell Type`, y=value)) + 
  geom_violin(position=position_dodge(width=0.4), alpha=0) + 
  geom_point(aes(color=Patient, shape=`Sorting Method`)) + 
  geom_boxplot(width=0.1, position = position_dodge(width=0.4), alpha=0) + 
  facet_wrap(.~Gene) + theme_classic() + ylab("Normalized Expression")
ggsave("CD14_CD34_expression_TMT.pdf", width=6, height=5)

ggplot2::ggplot(val.df[val.df$`Sorting Method` == "Bead",], aes(x=`Cell Type`, y=value)) + 
  geom_violin(position=position_dodge(width=0.4), alpha=0) + 
  geom_point(aes(color=Patient)) + 
  geom_boxplot(width=0.1, position = position_dodge(width=0.4), alpha=0) + 
  facet_wrap(.~Gene) + theme_classic() + ylab("Normalized Expression")
ggsave("CD14_CD34_expression_TMT_bead.pdf", width=6, height=5)

# require validation values to be higher than other values (e.g., CD34 for CD14+ sample should be less than CD14)
library(plyr);library(dplyr)
ratio.df <- plyr::ddply(val.df, .(variable, `Cell Type`, Patient, `Sorting Method`), summarize,
                        marker = ifelse(value[Gene=="CD14"] > value[Gene=="CD34"],"CD14","CD34"))
ratio.df <- na.omit(ratio.df)
ratio.df$outlier <- FALSE
ratio.df[ratio.df$`Cell Type` == "CD14+" & ratio.df$marker == "CD34",]$outlier <- TRUE
ratio.df[ratio.df$`Cell Type` == "CD34+" & ratio.df$marker == "CD14",]$outlier <- TRUE
outliers <- ratio.df[ratio.df$outlier,]$variable
length(outliers) # 3
outliers # X00839_CD34plus     X00839_CD34plusFlow X00117_CD34plus

ggplot2::ggplot(val.df[!(val.df$variable %in% outliers),], aes(x=`Cell Type`, y=value)) + 
  geom_violin(position=position_dodge(width=0.4), alpha=0) + 
  geom_point(aes(color=Patient, shape=`Sorting Method`)) + 
  geom_boxplot(width=0.1, position = position_dodge(width=0.4), alpha=0) + 
  facet_wrap(.~Gene) + theme_classic() + ylab("Normalized Expression")
ggsave("CD14_CD34_expression_TMT_noOutliers.pdf", width=6, height=5)

ggplot2::ggplot(val.df[val.df$`Sorting Method` == "Bead" & !(val.df$variable %in% outliers),], aes(x=`Cell Type`, y=value)) + 
  geom_violin(position=position_dodge(width=0.4), alpha=0) + 
  geom_point(aes(color=Patient)) + 
  geom_boxplot(width=0.1, position = position_dodge(width=0.4), alpha=0) + 
  facet_wrap(.~Gene) + theme_classic() + ylab("Normalized Expression")
ggsave("CD14_CD34_expression_TMT_bead_noOutliers.pdf", width=6, height=5)

tmt.wo.out <- list("meta" = meta.df[!(meta.df$DIA_id %in% outliers),],
                   "global" = tmt$global[,colnames(tmt$global)[!(colnames(tmt$global) %in% outliers)]],
                   "phospho" = tmt$phospho[,colnames(tmt$phospho)[!(colnames(tmt$phospho) %in% outliers)]])

# redo coverage again
tmt.wo.out.cov <- tmt.wo.out$global[, which(colMeans(!is.na(tmt.wo.out$global)) >= 0.75)] # 25 samples out of 25
tmt.wo.out.cov <- tmt.wo.out.cov[which(rowSums(is.na(tmt.wo.out.cov)) < ncol(tmt.wo.out.cov)/2),] # 5169 out of 5169 proteins
"CD14" %in% rownames(tmt.wo.out.cov) # TRUE
"CD34" %in% rownames(tmt.wo.out.cov) # TRUE
tmt.wo.out.covP <- tmt.wo.out$phospho[, which(colMeans(!is.na(tmt.wo.out$phospho)) >= 0.75)] # 18 samples out of 25
tmt.wo.out.covP <- tmt.wo.out.covP[which(rowSums(is.na(tmt.wo.out.covP)) < ncol(tmt.wo.out.covP)/2),] # 790 out of 794 sites
rownames(tmt.wo.out$meta) <- tmt.wo.out$meta$DIA_id
tmt.wo.out$meta$cellType <- "CD14"
tmt.wo.out$meta[grepl("CD34",tmt.wo.out$meta$SampleType),]$cellType <- "CD34"
tmt.wo.out$meta[grepl("MSC",tmt.wo.out$meta$SampleType),]$cellType <- "MSC"
tmt.wo.out <- list("meta" = tmt.wo.out$meta,
                   "global" = tmt.wo.out.cov,
                   "phospho" = tmt.wo.out.covP)
saveRDS(tmt.wo.out,"TMT_noOutliers_2025-07-07.rds")
tmt.wo.out <- readRDS("data/TMT_noOutliers_2025-07-07.rds")

### DIA
meta.df <- readxl::read_excel("Exp24metadataTable_DIA.xlsx") 
meta.df$id <- stringr::str_split_i(meta.df$patient, "-", -1)
meta.df$id <- paste0("X", meta.df$id)
meta.df[meta.df$`sample type` == "CD14+", ]$id <- paste0(meta.df[meta.df$`sample type` == "CD14+", ]$id, "_CD14plus")
meta.df[meta.df$`sample type` == "CD34+", ]$id <- paste0(meta.df[meta.df$`sample type` == "CD34+", ]$id, "_CD34plus")
meta.df[meta.df$`sample type` == "CD14+ Flow", ]$id <- paste0(meta.df[meta.df$`sample type` == "CD14+ Flow", ]$id, "_CD14plusFlow")
meta.df[meta.df$`sample type` == "CD34+ Flow", ]$id <- paste0(meta.df[meta.df$`sample type` == "CD34+ Flow", ]$id, "_CD34plusFlow")
meta.df[meta.df$`sample type` == "MSC Flow", ]$id <- paste0(meta.df[meta.df$`sample type` == "MSC Flow", ]$id, "_MSCflow")
rownames(meta.df) <- meta.df$id

# add other drug info & make sure sensitivity is correctly labeled
meta.df <- merge(meta.df, sens.info)
meta.df$X <- NULL

# add other metadata for contrasts
meta.df$CD14 <- "Neg"
meta.df[grepl("CD14+", meta.df$`sample type`),]$CD14 <- "Pos"

meta.df$CD34 <- "Neg"
meta.df[grepl("CD34+", meta.df$`sample type`),]$CD34 <- "Pos"

meta.df$MSC <- "Non_MSC"
meta.df[meta.df$`sample type` == "MSC Flow",]$MSC <- "MSC"

meta.df$'Sort Type' <- "Bead"
meta.df[grepl("Flow", meta.df$`sample type`),]$'Sort Type' <- "Flow"

synapser::synLogin()
globalFile <- synapser::synGet("syn63608740") # Samantha's unprocessed version
# previously: syn55234888 with 8897 proteins -> 6115 after requiring proteins in 50+% of samples
# now with larger library: syn63608740 with 7025 proteins after requiring proteins in 50+% of samples, 
global.df <- read.table(
  globalFile$path, 
  sep = "\t") # 7025 proteins, 48 samples
hist(rowSums(is.na(global.df)))
hist(colMeans(!is.na(global.df)))

# require proteins to be quantified in at least half of samples
global.df <- global.df[which(rowSums(is.na(global.df)) < ncol(global.df)/2),] # 7025 proteins, 48 samples

# require samples to have at least 50% of proteins quantified
#global.df <- global.df[ , which(colSums(is.na(global.df)) < nrow(global.df)/2)] # 42 out of 48 samples are kept
global.df75 <- global.df[ , which(colMeans(!is.na(global.df)) >= 0.75)] # 37 out of 48 samples are kept
hist(unlist(global.df75))

# add column for feature names and later make it the first column
global.df$Gene <- rownames(global.df) # 7025 gene symbols
global.df75$Gene <- rownames(global.df75) # 7025 gene symbols
rownames(meta.df) <- meta.df$id
meta.df$`Sample Type` <- meta.df$`sample type`
#write.csv(global.df, "Exp24_DIA_crosstab_global_gene_corrected.csv", row.names = FALSE)
#write.csv(global.df75, "Exp24_DIA_75PercentCoverage_crosstab_global_gene_corrected.csv", row.names = FALSE)
dia <- list("meta" = meta.df[meta.df$id %in% colnames(global.df), ],
            "global" = global.df)
dia75 <- list("meta" = meta.df[meta.df$id %in% colnames(global.df75), ],
            "global" = global.df75)

phospho.pep.DIA <- read.table("phospho_data/Exp24_DIA_crosstab_phospho_peptide_corrected.txt",sep = "\t") # 2331 peptides
pep <- list("meta" = meta.df,
            "TMT" = phospho.pep.TMT,
            "DIA" = phospho.pep.DIA)


### combine DIA & TMT
# #dia$meta$method <- "DIA"
# dia75$meta$method <- "DIA"
# tmt$meta$method <- "TMT"
# 
# # make sure dia & tmt meta data have same columns
# #dia$meta[, colnames(dia$meta)[!(colnames(dia$meta) %in% colnames(tmt$meta))]] <- NULL
# tmt$meta[, colnames(tmt$meta)[!(colnames(tmt$meta) %in% colnames(dia75$meta))]] <- NULL
# dia75$meta[, colnames(dia75$meta)[!(colnames(dia75$meta) %in% colnames(tmt$meta))]] <- NULL
# 
# dia.tmt75 <- list("meta" = rbind(dia75$meta, tmt$meta),
#                 "global" = merge(dia75$global, tmt$global, all = TRUE, 
#                                  suffixes = c("_DIA", "_TMT")))

method.data <- list("TMT" = tmt,
                    "DIA" = dia,
                    "DIA_75PercentCoverage" = dia75)

### new batch of DIA ###
meta.df <- readxl::read_excel("CPTAC Samples 03_25_25_long.xlsx") 
meta.df <- meta.df[meta.df$`run or not run?`=="run",]
meta.df$patient <- gsub("-",".",meta.df$Sample)
meta.df$id <- paste0("PTRC_",meta.df$patient,"_",meta.df$`cell type`,"_",
                     substring(meta.df$`sorting method`,1,1),"_",round(as.numeric(meta.df$cells)/1000,0),"K")
rownames(meta.df) <- meta.df$id

dia.meta.df <- readxl::read_excel("Exp24metadataTable_DIA.xlsx") 
dia.meta.df$id <- stringr::str_split_i(dia.meta.df$patient, "-", -1)
dia.meta.df$id <- paste0("X", dia.meta.df$id)
dia.meta.df[dia.meta.df$`sample type` == "CD14+", ]$id <- paste0(dia.meta.df[dia.meta.df$`sample type` == "CD14+", ]$id, "_CD14plus")
dia.meta.df[dia.meta.df$`sample type` == "CD34+", ]$id <- paste0(dia.meta.df[dia.meta.df$`sample type` == "CD34+", ]$id, "_CD34plus")
dia.meta.df[dia.meta.df$`sample type` == "CD14+ Flow", ]$id <- paste0(dia.meta.df[dia.meta.df$`sample type` == "CD14+ Flow", ]$id, "_CD14plusFlow")
dia.meta.df[dia.meta.df$`sample type` == "CD34+ Flow", ]$id <- paste0(dia.meta.df[dia.meta.df$`sample type` == "CD34+ Flow", ]$id, "_CD34plusFlow")
dia.meta.df[dia.meta.df$`sample type` == "MSC Flow", ]$id <- paste0(dia.meta.df[dia.meta.df$`sample type` == "MSC Flow", ]$id, "_MSCflow")
rownames(dia.meta.df) <- dia.meta.df$id
dia.meta.df$patient2 <- gsub("-","[.]",dia.meta.df$patient)

any(dia.meta.df$patient2 %in% meta.df$patient) # no patient overlap

# add other drug info & make sure sensitivity is correctly labeled
sens.info2 <- readxl::read_excel(synapser::synGet("syn65472730")$path)
colnames(sens.info2)[1] <- "Sample"
meta.df <- merge(meta.df, sens.info2,by="Sample")

# add other metadata for contrasts
meta.df$CD14 <- "Neg"
meta.df[grepl("cd14", meta.df$`cell type`),]$CD14 <- "Pos"

meta.df$CD34 <- "Neg"
meta.df[grepl("cd34", meta.df$`cell type`),]$CD34 <- "Pos"

meta.df$MSC <- "Non_MSC"
meta.df[meta.df$`cell type` == "msc",]$MSC <- "MSC"

meta.df$'Sort Type' <- "Bead"
meta.df[grepl("flow", meta.df$`sorting method`),]$'Sort Type' <- "Flow"

synapser::synLogin()
globalFile <- synapser::synGet("syn66694759") # from Samantha processed using: https://github.com/PNNL-CompBio/AML_sorted_proteomics/blob/main/proteomics/DIA/final_processing_DIA_sorted_cells.Rmd 
# she originally required proteins in 24+ samples because last batch had 48 samples, but now updated to ncol(m)/2 = 15.5 samples
global.df <- read.table(
  globalFile$path, 
  sep = "\t") # 8404 rows, 31 variables
hist(rowSums(is.na(global.df)))
hist(colMeans(!is.na(global.df)))

# require proteins to be quantified in at least half of samples
global.df <- global.df[which(rowSums(is.na(global.df)) < ncol(global.df)/2),] # 8404 proteins, 31 samples

# require samples to have at least 50% of proteins quantified
global.df75 <- global.df[ , which(colMeans(!is.na(global.df)) >= 0.75)] # 31 out of 31 samples are kept
hist(unlist(global.df75))
newCols <- sub("K_.*","K", colnames(global.df75))
colnames(global.df75) <- newCols
global.df75$Gene <- rownames(global.df75)
dia2 <- list("meta" = meta.df,
             "global" = global.df75)

## combine DIA batches
meta.df$Aza.Ven <- NA
meta.df$`Aza-Ven AUC` <- as.numeric(meta.df$`Aza-Ven AUC`)
meta.df[!is.na(meta.df$`Aza-Ven AUC`) & 
          meta.df$`Aza-Ven AUC`<100,]$Aza.Ven <- "Sensitive"
meta.df[!is.na(meta.df$`Aza-Ven AUC`) & 
          meta.df$`Aza-Ven AUC`>100,]$Aza.Ven <- "Resistant"

meta.df$Aza <- NA
meta.df$`Aza AUC` <- as.numeric(meta.df$`Aza AUC`)
meta.df[!is.na(meta.df$`Aza AUC`) & 
          meta.df$`Aza AUC`<200,]$Aza <- "Sensitive"
meta.df[!is.na(meta.df$`Aza AUC`) & 
          meta.df$`Aza AUC`>200,]$Aza <- "Resistant"

meta.df$Ven <- NA
meta.df$`Ven AUC` <- as.numeric(meta.df$`Ven AUC`)
meta.df[!is.na(meta.df$`Ven AUC`) & 
          meta.df$`Ven AUC`<100,]$Ven <- "Sensitive"
meta.df[!is.na(meta.df$`Ven AUC`) & 
          meta.df$`Ven AUC`>100,]$Ven <- "Resistant"
meta.df <- meta.df[,c("Sample","id","CD14","CD34","MSC","Sort Type","Aza.Ven","Ven","Aza")]
colnames(meta.df)[1] <- "patient"
meta.df$batch <- "batch2"
rownames(meta.df) <- meta.df$id
oldMeta <- dia75$meta[,c("patient","id","CD14","CD34","MSC","Sort Type","Aza.Ven","Ven","Aza")]
oldMeta$batch <- "batch1"
mCombo <- rbind(meta.df, oldMeta)
gCombo <- merge(global.df75, dia75$global, by="Gene", all=TRUE) # 8538 rows, 68 samples
mCombo$Flow <- FALSE
mCombo[mCombo$`Sort Type`=="Flow",]$Flow <- TRUE
rownames(gCombo) <- gCombo$Gene
gCombo <- gCombo[, which(colMeans(!is.na(gCombo)) >= 0.75)] # 60 samples out of 68
gCombo <- gCombo[which(rowSums(is.na(gCombo)) < ncol(gCombo)/2),] # 7009 proteins out of 8538

mCombo$cellType <- "CD34"
mCombo[mCombo$CD14=="Pos",]$cellType <- "CD14"
mCombo[mCombo$MSC=="MSC",]$cellType <- "MSC"
mCombo <- mCombo[colnames(gCombo)[colnames(gCombo) != "Gene"],]
diaCombo <- list("meta" = mCombo,
                 "global" = gCombo)
saveRDS(diaCombo, "DIA_2batches.rds")

#### 2. Histograms and PCA ####
base.path <- "~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/analysis"
setwd(base.path)

phenos <- c("Plex", "id", "patient", "Sample Type", "Aza", "Ven", "Aza.Ven", "Flow","batch","cellType")
#phenos <- "cellType"
library(MSnSet.utils)
write.csv(mCombo,"Exp24-27_metadata.csv")

dir.create("histograms_and_PCA")
setwd("histograms_and_PCA")
dir.create("20250707")
setwd("20250707")
#dir.create("20241001")
#setwd("20241001")
#dir.create("rmOutliers")
#setwd("rmOutliers")
#dir.create("rmMoreOutliers")
#setwd("rmMoreOutliers")
method.data <- list("DIA_2batches" = diaCombo)
method.data <- list("DIA_2batches_noOutliers" = dia.wo.out)
method.data <- list("TMT_noOutliers" = tmt.wo.out, "TMT" = tmt)
saveRDS(dia.wo.out,"DIA_2batches_noOutliers.rds")
for (k in 1:length(method.data)) {
  if (grepl("DIA", names(method.data)[k])) {
    omics <- list("Global" = method.data[[k]]$global)
  } else {
    omics <- list("Global" = method.data[[k]]$global,
                  "Phospho" = method.data[[k]]$phospho)
  }
  temp.meta <- method.data[[k]]$meta
  temp.phenos <- phenos[phenos %in% colnames(temp.meta)]
  
  for (i in 1:length(omics)) {
    # histogram
    melted.df <- reshape2::melt(omics[[i]])
    xlab <- paste("Normalized", names(omics)[i], "Expression")
    title <- names(method.data)[k]
    pdf(file.path(paste0(names(method.data)[k], "_", names(omics)[i], "_histogram_", Sys.Date(), ".pdf")))
    hist(melted.df$value, xlab = xlab, main = title)
    dev.off()
    
    # pca
    if (names(method.data)[k] != "DIA_&_TMT") {
      #sample.names <- temp.meta$id[temp.meta$id %in% colnames(omics[[i]])]
      if ("DIA_id" %in% colnames(temp.meta)) {
        sample.names <- colnames(omics[[i]])[colnames(omics[[i]]) %in% temp.meta$DIA_id]
      } else {
        sample.names <- colnames(omics[[i]])[colnames(omics[[i]]) %in% temp.meta$id]
      }
      pca.data <- MSnSet(exprs = omics[[i]][, sample.names] %>% as.matrix(),
                         pData = temp.meta[sample.names,])
      # write.table(exprs(pca.data),
      #             file = "~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/data/global_data/Exp24-27_crosstab_global_DIA_gene.txt",
      #             quote=F, sep="\t")
      for (j in 1:length(temp.phenos)) {
        MSnSet.utils::plot_pca(pca.data, phenotype = temp.phenos[j], label = "patient") + ggtitle(paste(names(method.data)[k], names(omics)[i], "PCA"))
        ggsave(paste0(names(method.data)[k], "_", names(omics)[i], "_PCA_", temp.phenos[j], "_", Sys.Date(), ".pdf"),width=5,height=5) # 3018 DIA (49%), 3975 (77%) TMT; wo outliers1: 3261 DIA (53%), 3975 TMT w or wo outliers (77%), 415 TMT phospho (2%) wo outliers or 14 w outliers
        
        if ("batch" %in% temp.phenos) {
          MSnSet.utils::plot_pca(pca.data, phenotype = temp.phenos[j], label = "batch") + ggtitle(paste(names(method.data)[k], names(omics)[i], "PCA"))
          ggsave(paste0(names(method.data)[k], "_", names(omics)[i], "_PCA_", temp.phenos[j], "_", Sys.Date(), "_batchLabel.pdf"),width=5,height=5) # 3018 DIA (49%), 3975 (77%) TMT; wo outliers1: 3261 DIA (53%)
        }
      
        MSnSet.utils::plot_pca(pca.data, phenotype = temp.phenos[j], label = "Sort Type") + ggtitle(paste(names(method.data)[k], names(omics)[i], "PCA"))
        ggsave(paste0(names(method.data)[k], "_", names(omics)[i], "_PCA_", temp.phenos[j], "_", Sys.Date(), "_flowLabel.pdf"),width=5,height=5) # 3018 DIA (49%), 3975 (77%) TMT; wo outliers1: 3261 DIA (53%), 3975 TMT w or wo outliers (77%), 415 TMT phospho (2%) wo outliers or 14 w outliers
        
        if ("cellType" %in% temp.phenos) {
        MSnSet.utils::plot_pca(pca.data, phenotype = temp.phenos[j], label = "cellType") + ggtitle(paste(names(method.data)[k], names(omics)[i], "PCA"))
        ggsave(paste0(names(method.data)[k], "_", names(omics)[i], "_PCA_", temp.phenos[j], "_", Sys.Date(), "_cellLabel.pdf"),width=5,height=5) # 3018 DIA (49%), 3975 (77%) TMT; wo outliers1: 3261 DIA (53%)
        } else {
          MSnSet.utils::plot_pca(pca.data, phenotype = temp.phenos[j], label = "Sample Type") + ggtitle(paste(names(method.data)[k], names(omics)[i], "PCA"))
          ggsave(paste0(names(method.data)[k], "_", names(omics)[i], "_PCA_", temp.phenos[j], "_", Sys.Date(), "_cellLabel.pdf"),width=5,height=5) # 3018 DIA (49%), 3975 (77%) TMT; wo outliers1: 3261 DIA (53%), 3975 TMT w or wo outliers (77%), 415 TMT phospho (2%) wo outliers or 14 w outliers
          
        }
      } 
      
      # try to resolve batch effect
      if ("batch" %in% temp.phenos){
        pca.data <- correct_batch_effect_NA(pca.data, "batch", "cellType", par.prior = T)
        # write.table(exprs(pca.data),
        #             file = "~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/data/global_data/Exp24-27_crosstab_global_DIA_gene_correctedConsideringCellType.txt",
        #             quote=F, sep="\t")
        
        for (j in 1:length(temp.phenos)) {
          MSnSet.utils::plot_pca(pca.data, phenotype = temp.phenos[j], label = "patient") + ggtitle(paste(names(method.data)[k], names(omics)[i], "PCA"))
          ggsave(paste0(names(method.data)[k], "_", names(omics)[i], "_PCA_", temp.phenos[j], "_", Sys.Date(), "_batchCorrectedTypeCov.pdf"),width=5,height=5) # 3018 DIA (49%), 3975 (77%) TMT; wo outliers1: 3261 DIA (53%), 3986 TMT (77%), 17 TMT phospho (2%)
          
          MSnSet.utils::plot_pca(pca.data, phenotype = temp.phenos[j], label = "batch") + ggtitle(paste(names(method.data)[k], names(omics)[i], "PCA"))
          ggsave(paste0(names(method.data)[k], "_", names(omics)[i], "_PCA_", temp.phenos[j], "_", Sys.Date(), "_batchCorrectedTypeCov_batchLabel.pdf"),width=5,height=5) # 3018 DIA (49%), 3975 (77%) TMT; wo outliers1: 3261 DIA (53%), 3986 TMT (77%), 17 TMT phospho (2%)
          
          MSnSet.utils::plot_pca(pca.data, phenotype = temp.phenos[j], label = "Sort Type") + ggtitle(paste(names(method.data)[k], names(omics)[i], "PCA"))
          ggsave(paste0(names(method.data)[k], "_", names(omics)[i], "_PCA_", temp.phenos[j], "_", Sys.Date(), "_batchCorrectedTypeCov_flowLabel.pdf"),width=5,height=5) # 3018 DIA (49%), 3975 (77%) TMT; wo outliers1: 3261 DIA (53%), 3986 TMT (77%), 17 TMT phospho (2%)
          
        } 
      }
    }
  }
}
# wo outliers2: 3368 DIA (55%), 4793 TMT (93%), 85 TMT phospho (11%)
# wo outliers2a: 3408 DIA (56%), 4860 TMT (94%), 122 TMT phospho (15%)
# not sample centered DIA, no outliers removed: 3018 DIA (49%), 3975 (77%) TMT, 14 TMT phospho (2%)
# 1917 complete rows for DIA, 3975 complete rows for TMT
# 4249 complete rows for DIA 2 batches after redoing coverage filters after combining
# 4314 complete rows after removing outliers based on CD14, CD34 and then redoing coverage filters again

# why are there some outliers?
cov.df <- reshape2::melt(gCombo)
cov.df <- plyr::ddply(cov.df, .(variable), summarize,
                      coverage = length(value[!is.na(value)])/length(value))
# some of the lowest coverage samples appear to be outliers
val.df <- reshape2::melt(gCombo[c("CD14","CD34"),])
val.df$`Cell Type` <- "CD14+"
val.df[grepl("cd34", val.df$variable, ignore.case=TRUE),]$`Cell Type` <- "CD34+"
val.df[grepl("msc", val.df$variable, ignore.case=TRUE),]$`Cell Type` <- "MSC"
val.df$Batch <- "1"
val.df[grepl("PTRC_",val.df$variable),]$Batch <- "2"
val.df$`Sorting Method` <- "Bead"
val.df[grepl("flow",val.df$variable, ignore.case=TRUE) | grepl("_f_", val.df$variable),]$`Sorting Method` <- "Flow"
val.df$Patient <- stringr::str_split_i(val.df$variable, "_", 1)
val.df[val.df$Patient=="PTRC",]$Patient <- stringr::str_split_i(val.df[val.df$Patient=="PTRC",]$variable, "_", 2)
ggplot2::ggplot(val.df, aes(x=`Cell Type`, y=value)) + 
  geom_violin(position=position_dodge(width=0.4), alpha=0) + 
  geom_point(aes(color=Patient, shape=`Sorting Method`)) + 
  geom_boxplot(width=0.1, position = position_dodge(width=0.4), alpha=0) + 
  facet_wrap(.~Gene) + theme_classic() + ylab("Normalized Expression")
ggsave("CD14_CD34_expression_DIA_2batches.pdf", width=6, height=5)

ggplot2::ggplot(val.df[val.df$`Sorting Method` == "Bead",], aes(x=`Cell Type`, y=value)) + 
  geom_violin(position=position_dodge(width=0.4), alpha=0) + 
  geom_point(aes(color=Patient)) + 
  geom_boxplot(width=0.1, position = position_dodge(width=0.4), alpha=0) + 
  facet_wrap(.~Gene) + theme_classic() + ylab("Normalized Expression")
ggsave("CD14_CD34_expression_DIA_2batches_bead.pdf", width=6, height=5)

# # maybe require CD14 to be detected in CD14+ samples and CD34 in CD34+ samples
# cd14.fail <- val.df[val.df$Gene == "CD14" & val.df$`Cell Type` == "CD14+" & is.na(val.df$value),]$variable
# length(cd14.fail) # 0
# cd34.fail <- val.df[val.df$Gene == "CD34" & val.df$`Cell Type` == "CD34+" & is.na(val.df$value),]$variable
# length(cd34.fail) # 0
# 
# # maybe require validation values > 0
# cd14.fail <- val.df[val.df$Gene == "CD14" & val.df$`Cell Type` == "CD14+" & val.df$value <= 0,]$variable
# length(cd14.fail) # 2: PTRC_18.00390_cd14_b_20K PTRC_18.00290_cd14_b_20K
# cd34.fail <- val.df[val.df$Gene == "CD34" & val.df$`Cell Type` == "CD34+" & val.df$value <= 0,]$variable
# length(cd34.fail) # 1: PTRC_20.00450_cd34_b_20K
# 
# # also require other values be less than median of validation group (e.g., CD34 for CD14+ sample should be less than CD34 median of CD34+ samples)
# cd14.med <- median(val.df[val.df$Gene == "CD14" & val.df$`Cell Type` == "CD14+",]$value) # 3.33
# cd34.med <- median(val.df[val.df$Gene == "CD34" & val.df$`Cell Type` == "CD34+",]$value) # 3.37
# cd14.fail.med <- na.omit(val.df[val.df$Gene == "CD34" & val.df$`Cell Type` == "CD14+" & val.df$value >= cd34.med,]$variable)
# length(cd14.fail.med) # 4: PTRC_20.00083_cd14_b_20K PTRC_18.00290_cd14_b_20K X00103_CD14plus          X00432_CD14plus
# cd34.fail.med <- na.omit(val.df[val.df$Gene == "CD14" & val.df$`Cell Type` == "CD34+" & val.df$value >= cd14.med,]$variable)
# length(cd34.fail.med) # 0
# msc.fail.med <- na.omit(c(val.df[val.df$Gene == "CD34" & val.df$`Cell Type` == "MSC" & val.df$value >= cd34.med,]$variable,
#                   val.df[val.df$Gene == "CD14" & val.df$`Cell Type` == "MSC" & val.df$value >= cd14.med,]$variable))
# length(msc.fail.med) # 0

# require validation values to be higher than other values (e.g., CD34 for CD14+ sample should be less than CD14)
ratio.df <- plyr::ddply(val.df, .(variable, `Cell Type`, Patient, `Sorting Method`, Batch), summarize,
                        marker = ifelse(value[Gene=="CD14"] > value[Gene=="CD34"],"CD14","CD34"))
ratio.df <- na.omit(ratio.df)
ratio.df$outlier <- FALSE
ratio.df[ratio.df$`Cell Type` == "CD14+" & ratio.df$marker == "CD34",]$outlier <- TRUE
ratio.df[ratio.df$`Cell Type` == "CD34+" & ratio.df$marker == "CD14",]$outlier <- TRUE
outliers <- ratio.df[ratio.df$outlier,]$variable
length(outliers) # 9
outliers
# outliers are: PTRC_18.00390_cd14_b_20K PTRC_20.00450_cd34_b_20K 
#PTRC_20.00083_cd14_b_20K PTRC_18.00290_cd14_b_20K PTRC_21.00176_cd14_b_15K 
# PTRC_21.00034_cd34_b_20K X00103_CD14plus X00432_CD14plus X00839_CD34plusFlow

#outliers <- unique(c(cd14.fail, cd34.fail, cd14.fail.med))

ggplot2::ggplot(val.df[!(val.df$variable %in% outliers),], aes(x=`Cell Type`, y=value)) + 
  geom_violin(position=position_dodge(width=0.4), alpha=0) + 
  geom_point(aes(color=Patient, shape=`Sorting Method`)) + 
  geom_boxplot(width=0.1, position = position_dodge(width=0.4), alpha=0) + 
  facet_wrap(.~Gene) + theme_classic() + ylab("Normalized Expression")
ggsave("CD14_CD34_expression_DIA_2batches_noOutliers.pdf", width=6, height=5)

ggplot2::ggplot(val.df[val.df$`Sorting Method` == "Bead" & !(val.df$variable %in% outliers),], aes(x=`Cell Type`, y=value)) + 
  geom_violin(position=position_dodge(width=0.4), alpha=0) + 
  geom_point(aes(color=Patient)) + 
  geom_boxplot(width=0.1, position = position_dodge(width=0.4), alpha=0) + 
  facet_wrap(.~Gene) + theme_classic() + ylab("Normalized Expression")
ggsave("CD14_CD34_expression_DIA_2batches_bead_noOutliers.pdf", width=6, height=5)

# # determine outliers based on PCAs and iterate
# #cd14outliers <- c("18-00390", "18-00290", "20-00083","21-00176","22-00251","21-00432")
# cd14outliers <- c("20-00083", "18-00290", "18-00390","21-00176","21-00432"
#                   #,"22-00251" # don't see it after redoing coverage filter
#                   # could also remove 18-00103
#                   )
# cd14outlierIDs <- mCombo[mCombo$patient %in% cd14outliers & mCombo$CD14=="Pos" & mCombo$`Sort Type`=="Bead",]$id # need to specify bead because otherwise 2 samples from 21-00432 get removed
# cd34outliers <- c("21-00034","20-00450")
# cd34outlierIDs <- mCombo[mCombo$patient %in% cd34outliers & mCombo$CD34=="Pos" & mCombo$`Sort Type`=="Bead",]$id
# mscOutliers <- c("18-00190", "18-00390")
# mscOutlierIDs <- mCombo[mCombo$patient %in% mscOutliers & mCombo$MSC=="MSC",]$id
# outliers <- c(cd14outlierIDs, cd34outlierIDs, mscOutlierIDs) # 9
# #dia <- exprs(pca.data)
dia <- read.table(
  synapser::synGet("syn66695021")$path, 
  sep = "\t") # corrected accounting for cell type

dia.wo.out <- list("meta" = mCombo[!(mCombo$id %in% outliers),],
                   "global" = dia[,colnames(dia)[!(colnames(dia) %in% outliers)]]) # 7009 proteins (same), 51 samples down from 60

# redo coverage again
dia.wo.out.cov <- dia.wo.out$global[, which(colMeans(!is.na(dia.wo.out$global)) >= 0.75)] # 51 samples out of 51
dia.wo.out.cov <- dia.wo.out.cov[which(rowSums(is.na(dia.wo.out.cov)) < ncol(dia.wo.out.cov)/2),] # 6887 out of 7009 proteins
"CD14" %in% rownames(dia.wo.out.cov) # TRUE
"CD34" %in% rownames(dia.wo.out.cov) # TRUE
dia.wo.out <- list("meta" = dia.wo.out$meta,
                   "global" = dia.wo.out.cov)
saveRDS(dia.wo.out,"DIA_2batches_noOutliers.rds")
dia.wo.out <- readRDS("analysis/DIA_2batches_noOutliers.rds")

method.data <- list("DIA_2batches_noOutliers" = dia.wo.out)
for (k in 1:length(method.data)) {
  if (grepl("DIA", names(method.data)[k])) {
    omics <- list("Global" = method.data[[k]]$global)
  } else {
    omics <- list("Global" = method.data[[k]]$global,
                  "Phospho" = method.data[[k]]$phospho)
  }
  temp.meta <- method.data[[k]]$meta
  temp.phenos <- phenos[phenos %in% colnames(temp.meta)]
  
  for (i in 1:length(omics)) {
    # histogram
    melted.df <- reshape2::melt(omics[[i]])
    xlab <- paste("Normalized", names(omics)[i], "Expression")
    title <- names(method.data)[k]
    pdf(file.path(paste0(names(method.data)[k], "_", names(omics)[i], "_histogram_", Sys.Date(), ".pdf")))
    hist(melted.df$value, xlab = xlab, main = title)
    dev.off()
    
    # pca
    if (names(method.data)[k] != "DIA_&_TMT") {
      #sample.names <- temp.meta$id[temp.meta$id %in% colnames(omics[[i]])]
      sample.names <- colnames(omics[[i]])[colnames(omics[[i]]) %in% temp.meta$id]
      pca.data <- MSnSet(exprs = omics[[i]][, sample.names] %>% as.matrix(),
                         pData = temp.meta[sample.names,])
      for (j in 1:length(temp.phenos)) {
        MSnSet.utils::plot_pca(pca.data, phenotype = temp.phenos[j], label = "patient") + ggtitle(paste(names(method.data)[k], names(omics)[i], "PCA"))
        ggsave(paste0(names(method.data)[k], "_", names(omics)[i], "_PCA_", temp.phenos[j], "_", Sys.Date(), "_correctedHighCov_noOutliers.pdf"),width=5,height=5) # 4416 complete rows
        
        MSnSet.utils::plot_pca(pca.data, phenotype = temp.phenos[j], label = "batch") + ggtitle(paste(names(method.data)[k], names(omics)[i], "PCA"))
        ggsave(paste0(names(method.data)[k], "_", names(omics)[i], "_PCA_", temp.phenos[j], "_", Sys.Date(), "_correctedHighCov_noOutliers_batchLabel.pdf"),width=5,height=5)
        
        MSnSet.utils::plot_pca(pca.data, phenotype = temp.phenos[j], label = "Sort Type") + ggtitle(paste(names(method.data)[k], names(omics)[i], "PCA"))
        ggsave(paste0(names(method.data)[k], "_", names(omics)[i], "_PCA_", temp.phenos[j], "_", Sys.Date(), "_correctedHighCov_noOutliers_flowLabel.pdf"),width=5,height=5)
        
        MSnSet.utils::plot_pca(pca.data, phenotype = temp.phenos[j], label = "cellType") + ggtitle(paste(names(method.data)[k], names(omics)[i], "PCA"))
        ggsave(paste0(names(method.data)[k], "_", names(omics)[i], "_PCA_", temp.phenos[j], "_", Sys.Date(), "_correctedHighCov_noOutliers_cellLabel.pdf"),width=5,height=5)
      } 
    }
  }
}

outliers <- c("X00839_CD34plusFlow", "X00117_CD34plus", 
              "X00432_CD14plus", "X00251_CD14plus", "X00105_CD14plusFlow")
outlier.ids <- c("X11", "X17", "X7", "X23", "X26")
tmt.wo.out <- list("meta" = tmt$meta[!(tmt$meta$id %in% outlier.ids),],
                   "global" = tmt$global[,colnames(tmt$global)[!(colnames(tmt$global) %in% outlier.ids)]],
                   "phospho" = tmt$phospho[,colnames(tmt$phospho)[!(colnames(tmt$phospho) %in% outlier.ids)]]) #5169 proteins (same), 794 phospho sites (same), 24 samples down from 28

method.data <- list("DIA" = dia.wo.out,
                    "TMT" = tmt.wo.out)
mCombo <- read.csv("Exp24-27_metadata.csv")
dia.wo.out <- readRDS("DIA_2batches_noOutliers.rds")
method.data <- list("DIA_2batches_noOutliers" = dia.wo.out)
length(na.omit(unique(dia.wo.out$meta$patient))) # 23
length(na.omit(unique(dia.wo.out$meta[dia.wo.out$meta$batch=="batch1",]$patient))) # 10
length(na.omit(unique(dia.wo.out$meta[dia.wo.out$meta$batch=="batch2",]$patient))) # 13

length(na.omit(unique(dia.wo.out$meta[dia.wo.out$meta$`Sort Type`=="Bead",]$patient))) # 20
length(na.omit(unique(dia.wo.out$meta[dia.wo.out$meta$batch=="batch1" &
                                        dia.wo.out$meta$`Sort Type`=="Bead",]$patient))) # 10
length(na.omit(unique(dia.wo.out$meta[dia.wo.out$meta$batch=="batch2" &
                                        dia.wo.out$meta$`Sort Type`=="Bead",]$patient))) # 10

phenos <- c("Plex", "id", "patient", "Sample Type", "Aza", "Ven", "Aza.Ven", "Flow","batch","cellType")
# also generate PCA just for DIA bead sorted (no flow, and therefore no MSC)
for (k in 1) {
  if (grepl("DIA", names(method.data)[k])) {
    omics <- list("Global" = method.data[[k]]$global)
  } else {
    omics <- list("Global" = method.data[[k]]$global,
                  "Phospho" = method.data[[k]]$phospho)
  }
  temp.meta <- method.data[[k]]$meta
  temp.meta <- temp.meta[temp.meta$`Sort Type` == "Bead",]
  #temp.phenos <- c("Sample Type")
  temp.phenos <- phenos[phenos %in% colnames(temp.meta)]
  
  for (i in 1:length(omics)) {
    # histogram
    melted.df <- reshape2::melt(omics[[i]])
    xlab <- paste("Normalized", names(omics)[i], "Expression")
    title <- names(method.data)[k]
    pdf(file.path(paste0(names(method.data)[k], "_beadSorted_", names(omics)[i], "_histogram_", Sys.Date(), ".pdf")))
    hist(melted.df$value, xlab = xlab, main = title)
    dev.off()
    
    # pca
    if (names(method.data)[k] != "DIA_&_TMT") {
      #sample.names <- temp.meta$id[temp.meta$id %in% colnames(omics[[i]])]
      sample.names <- colnames(omics[[i]])[colnames(omics[[i]]) %in% temp.meta$id]
      pca.data <- MSnSet(exprs = omics[[i]][, sample.names] %>% as.matrix(),
                         pData = temp.meta[sample.names,])
      for (j in 1:length(temp.phenos)) {
        MSnSet.utils::plot_pca(pca.data, phenotype = temp.phenos[j], label = "patient") + ggtitle(paste(names(method.data)[k], names(omics)[i], "PCA")) # 4612 complete rows for PCA
        ggsave(paste0(names(method.data)[k], "_beadSorted_", names(omics)[i], "_PCA_", temp.phenos[j], "_", Sys.Date(), "_correctedHighCov_noOutliers.pdf"),
               width=5, height=5) # 4670 2 DIA batches no outliers, high coverage; 3018 DIA (49%), 3975 (77%) TMT; wo outliers1: 3261 DIA (53%), 3986 TMT (77%), 17 TMT phospho (2%)
      } 
    }
  }
}

# also run DIA phospho peptide vs. TMT
#pep.wo.out <- list("meta" = pep$meta[!(pep$meta$id %in% outliers),],"DIA" = )
base.path <- "~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/analysis"
setwd(base.path)

phenos <- c("Plex", "id", "patient", "Sample Type", "Aza", "Ven", "Aza.Ven", "Flow")
library(MSnSet.utils)

dir.create("histograms_and_PCA")
setwd("histograms_and_PCA")
dir.create("phosphoPeptide")
setwd("phosphoPeptide")
omics <- pep[c("DIA","TMT")]
temp.meta <- pep$meta
temp.phenos <- phenos[phenos %in% colnames(temp.meta)]

for (i in 1:length(omics)) {
  # histogram
  melted.df <- reshape2::melt(omics[[i]])
  xlab <- paste("Normalized", "Phospho Peptide", "Expression")
  title <- names(omics)[i]
  pdf(file.path(paste0(title, "_histogram_", Sys.Date(), ".pdf")))
  hist(melted.df$value, xlab = xlab, main = title)
  dev.off()
  
  # pca
  #sample.names <- temp.meta$id[temp.meta$id %in% colnames(omics[[i]])]
  sample.names <- unique(colnames(omics[[i]])[colnames(omics[[i]]) %in% temp.meta$id])
  pca.data <- MSnSet(exprs = omics[[i]][, sample.names] %>% as.matrix(),
                     pData = temp.meta[sample.names,])
  complete.df <- omics[[i]][, sample.names][which(rowMeans(!is.na(omics[[i]][, sample.names])) == 1),]
  nCompleteRows <- nrow(complete.df)
  title <- paste0(names(omics)[i], " Phospho PCA (", nCompleteRows, "/", nrow(omics[[i]][, sample.names]), " peptides represented)")
  for (j in 1:length(temp.phenos)) {
    MSnSet.utils::plot_pca(pca.data, phenotype = temp.phenos[j], label = "patient") + ggtitle(title)
    ggsave(paste0(names(omics)[i], "_PCA_", temp.phenos[j], "_", Sys.Date(), ".pdf")) # 3018 DIA (49%), 3975 (77%) TMT
  } 
  
  non.outliers <- sample.names[!(sample.names %in% outliers)]
  pca.data <- MSnSet(exprs = omics[[i]][, non.outliers] %>% as.matrix(),
                     pData = temp.meta[non.outliers,])
  nonOutlier.df <- omics[[i]][, non.outliers]
  complete.df <- nonOutlier.df[which(rowMeans(!is.na(nonOutlier.df)) == 1),]
  nCompleteRows <- nrow(complete.df)
  title <- paste0(names(omics)[i], " Phospho PCA (", nCompleteRows, "/", nrow(nonOutlier.df), " peptides represented)")
  for (j in 1:length(temp.phenos)) {
    MSnSet.utils::plot_pca(pca.data, phenotype = temp.phenos[j], label = "patient") + ggtitle(title)
    ggsave(paste0(names(omics)[i], "_PCA_", temp.phenos[j], "_woOutliers_", Sys.Date(), ".pdf")) # 3018 DIA (49%), 3975 (77%) TMT
  } 
  
  # try filtering for samples with 75+% of peptides quantified
  highCoverage.df <- omics[[i]][,which(colMeans(!is.na(omics[[i]])) >= 0.75)] # 34 out of 48 samples kept
  highCoverage.samples <- sample.names[sample.names %in% colnames(highCoverage.df)]
  pca.data <- MSnSet(exprs = highCoverage.df[, highCoverage.samples] %>% as.matrix(),
                     pData = temp.meta[highCoverage.samples,])
  complete.df <- highCoverage.df[which(rowMeans(!is.na(highCoverage.df)) == 1),]
  nCompleteRows <- nrow(complete.df)
  title <- paste0(names(omics)[i], " Phospho PCA (", nCompleteRows, "/", nrow(highCoverage.df), " peptides represented)")
  for (j in 1:length(temp.phenos)) {
    MSnSet.utils::plot_pca(pca.data, phenotype = temp.phenos[j], label = "patient") + ggtitle(title)
    ggsave(paste0(names(omics)[i], "_PCA_", temp.phenos[j], "_min75PercentCoverage_", Sys.Date(), ".pdf")) # 1789 DIA (49%), 3975 (77%) TMT
  }
  
  # also try removing outliers based on global PCAs
  highCoverage.nonOutliers <- highCoverage.samples[!(highCoverage.samples %in% outliers)]
  pca.data <- MSnSet(exprs = highCoverage.df[, highCoverage.nonOutliers] %>% as.matrix(),
                     pData = temp.meta[highCoverage.nonOutliers,])
  complete.df <- highCoverage.df[, highCoverage.nonOutliers][which(rowMeans(!is.na(highCoverage.df[, highCoverage.nonOutliers])) == 1),]
  nCompleteRows <- nrow(complete.df)
  title <- paste0(names(omics)[i], " Phospho PCA (", nCompleteRows, "/", nrow(highCoverage.df[, highCoverage.nonOutliers]), " peptides represented)")
  for (j in 1:length(temp.phenos)) {
    MSnSet.utils::plot_pca(pca.data, phenotype = temp.phenos[j], label = "patient") + ggtitle(title)
    ggsave(paste0(names(omics)[i], "_PCA_", temp.phenos[j], "_min75PercentCoverage_woOutliers_", Sys.Date(), ".pdf")) # 1789 DIA (49%), 3975 (77%) TMT
  }
}
sample.names <- rownames(pep$meta)
non.outliers <- sample.names[!(sample.names %in% outliers)] # 44

# filter for high coverage
highCov.pep <- list("meta" = pep$meta,
                    "DIA" = pep$DIA[,which(colMeans(!is.na(pep$DIA)) >= 0.75)],
                    "TMT" = pep$TMT[,which(colMeans(!is.na(pep$TMT)) >= 0.75)])

non.outliers.DIA <- colnames(highCov.pep$DIA)[colnames(highCov.pep$DIA) %in% non.outliers] # 29
non.outliers.TMT <- colnames(highCov.pep$TMT)[colnames(highCov.pep$TMT) %in% non.outliers] # 18
non.outliers.overlapping <- non.outliers.DIA[non.outliers.DIA %in% non.outliers.TMT] # 17
highCov.pep.wo.out <- list("meta" = pep$meta[non.outliers.overlapping,],
                   "DIA" = highCov.pep$DIA[,c("SUB_SITE_MOD",non.outliers.overlapping)],
                   "TMT" = highCov.pep$TMT[,c("SUB_SITE_MOD",non.outliers.overlapping)])
saveRDS(highCov.pep.wo.out, "Phospho_peptide_DIA_TMT_overlappingSamples_min75PercentCoverage_noGlobalOutliers.rds")

#### 3. Run panSEA for each contrast & combination ####
setwd(base.path)
if (file.exists("gmt_BeatAML_drug_MOA.rds")) {
  gmt.drug <- readRDS("gmt_BeatAML_drug_MOA.rds")
} else {
  gmt.drug <- DMEA::as_gmt(moa.BeatAML, sep = ", ")
  saveRDS(gmt.drug, "gmt_BeatAML_drug_MOA.rds")
}

synapse_id <- "syn66698496"
all.degs <- data.frame()
all.degs.noMSC <- data.frame()
contrasts <- c("CD14", "CD34", #"Aza", "Ven", "Aza.Ven", # Aza + Ven has NA- need to debug why comparison isn't sensitive vs. resistant
               "Ven", "Sort Type", "MSC")
contrasts.noMSC <- c("CD14", "Sort Type", #"Aza", "Ven", "Aza.Ven"
                     "Ven")
contrasts <- c("CD14", "CD34", "Aza", "Ven", "Aza.Ven", 
               "Sort Type", "MSC")
contrasts.noMSC <- c("CD14", "Sort Type", "Aza", "Ven", "Aza.Ven")
#contrasts.noMSC <- c("CD14", "Sort Type")
#BeatAML.data <- load_not_norm_BeatAML_for_DMEA2()
#sorted.patients <- unique(mCombo$patient)

gmt1 <- readRDS("gmt1_more.rds")
dir.create("combined24-27")
setwd("combined24-27")
# gmt1 <- get_gmt1_v2()
# gmt1[[11]] <- NULL
# gmt[[13]] <- NULL
# names(gmt1) <- c("TFT_GTRD", "MIR_MIRDB", "GO_BP", "GO_CC", "GO_MF", 
#                  "Oncogenic_signatures", "BioCarta", "KEGG", "PID", "Reactome", 
#                  "Hallmark", "Positional")
# temp.gmt <- msigdbr::msigdbr(species = "Homo sapiens", category="C2", subcategory="CP:WIKIPATHWAYS")
# temp.gmt <- DMEA::as_gmt(as.data.frame(temp.gmt), element.names = "gene_symbol", set.names = "gs_name", descriptions = "gs_description")
# gmt1[["WikiPathways"]] <- temp.gmt
# gmt1 <- gmt1[c(1:12,14)]
# saveRDS(gmt1, "gmt1_more.rds")

# gmt.names <- c("Hallmark","PID", "Oncogenic_signatures", "KEGG")
# gmt1 <- gmt1[gmt.names]
# saveRDS(gmt1, "gmt1_more.rds")
#devtools::install_github("cstawitz/roomba")
# library(roomba)
# packageurl <- "http://cran.r-project.org/src/contrib/Archive/ggplot2/ggplot2_3.4.4.tar.gz"
# install.packages(packageurl, repos=NULL, type="source")
# library(ggplot2)
base.path <- getwd()
dia.wo.out$global <- as.data.frame(dia.wo.out$global)
dia.wo.out$global$Gene <- rownames(dia.wo.out$global)
unmatched <- rownames(dia.wo.out$meta)[!(rownames(dia.wo.out$meta) %in% colnames(dia.wo.out$global))] # none

method.data <- list("DIA_2batches_noOutliers"=dia.wo.out)
method.data <- list("TMT_noOutliers"=tmt.wo.out)
synapser::synLogin()
for (k in 1:length(method.data)) {
  setwd(base.path)
  method.path <- file.path(base.path, names(method.data)[k])
  dir.create(names(method.data)[k])
  setwd(names(method.data)[k])
  methodFolder <- 
    synapser::synStore(synapser::Folder(names(method.data)[k],
                                        parent = synapse_id))
  
  meta.df <- method.data[[k]]$meta
  if (names(method.data[k]) != "DIA") {
    meta.df$id <- meta.df$DIA_id
  }
  meta.df$Patient <- meta.df$patient
  meta.df$`Sample Type` <- meta.df$cellType
  
  sorted.patients <- unique(meta.df$patient)
  BeatAML.data <- load_not_norm_BeatAML_for_DMEA3(exclude.samples = sorted.patients)
  
  # run contrast combos
  if (grepl("DIA", names(method.data)[k])) {
    omics <- list("global" = method.data[[k]]$global)
    feature.names <- "Gene"
    temp.expr <- list(BeatAML.data$global)
  } else {
    omics <- list("global" = method.data[[k]]$global)
    feature.names <- "Gene"
    temp.expr <- list(BeatAML.data$global)
    
    # # don't use phospho since there were low cell counts
    # omics <- list("global" = method.data[[k]]$global,
    #               "phospho" = method.data[[k]]$phospho)
    # feature.names <- c("Gene", "SUB_SITE")
    # temp.expr <- list(BeatAML.data$global, BeatAML.data$phospho)
  }
  names(temp.expr) <- names(omics)
  panSEA2_combos2(contrasts, meta.df = meta.df,
                      omics = omics, 
                      expr = temp.expr,
                      gmt.drug = gmt.drug, drug.sens = BeatAML.data$drug,
                      base.path = base.path,
                      temp.path = method.path,
                      synapse_id = methodFolder)
  # Calculating Weighted Voting scores...
  # Error in 2:ncol(filtered.expr) : argument of length 0

  # get compiled DEGs
  methodDEGs <- as.list(synapser::synGetChildren(methodFolder, list("file"), sortBy = 'NAME'))
  if (length(methodDEGs) > 0) {
    if (methodDEGs[[1]]$name == "Differential_expression_results.csv") {
      methodFile <- synapser::synGet(methodDEGs[[1]]$id)
      methodDEG <- read.csv(methodFile$path)
      methodDEG$method <- names(method.data)[k]
      all.degs <- rbind(all.degs, methodDEG)
    }
  }
  
  # also filter for non-MSC
  setwd(base.path)
  method.path.noMSC <- file.path(base.path, paste0(names(method.data)[k],"_noMSC"))
  dir.create(paste0(names(method.data)[k],"_noMSC"))
  setwd(paste0(names(method.data)[k],"_noMSC"))
  methodFolder.noMSC <- 
    synapser::synStore(synapser::Folder(paste0(names(method.data)[k],"_noMSC"),
                                        parent = synapse_id))
  meta.df.noMSC <- meta.df[meta.df$MSC == "Non_MSC",]
  panSEA2_combos2(contrasts.noMSC, meta.df = meta.df.noMSC, 
                  omics = omics,
                  expr = temp.expr,
                  gmt.drug = gmt.drug, drug.sens = BeatAML.data$drug,
                  base.path = base.path,
                  temp.path = method.path.noMSC,
                  synapse_id = methodFolder.noMSC)
  # TMT no MSC: (just ran again from line above after)
  # [1] "Running Aza_Sensitive_vs_Resistant with Sort Type == Flow"
  # Error in .ebayes(fit = fit, proportion = proportion, stdev.coef.lim = stdev.coef.lim,  : 
  #                    No residual degrees of freedom in linear model fits
  #                  In addition: There were 50 or more warnings (use warnings() to see the first 50)
  
  # get compiled DEGs for analyses without MSCs
  methodDEGs.noMSC <- as.list(synapser::synGetChildren(methodFolder.noMSC, list("file"), sortBy = 'NAME'))
  if (length(methodDEGs.noMSC) > 0) {
    if (methodDEGs.noMSC[[1]]$name == "Differential_expression_results.csv") {
      methodFile <- synapser::synGet(methodDEGs.noMSC[[1]]$id)
      methodDEG <- read.csv(methodFile$path)
      methodDEG$method <- names(method.data)[k]
      all.degs.noMSC <- rbind(all.degs.noMSC, methodDEG)
    }
  }
}
setwd(base.path)
all.DEG.files <- list("Differential_expression_results.csv" = 
                        all.degs,
                      "Differential_expression_results_max_5_percent_FDR.csv" = 
                        all.degs[all.degs$adj.P.Val <= 0.05, ],
                      "Differential_expression_results_noMSC.csv" = 
                        all.degs.noMSC,
                      "Differential_expression_results_max_5_percent_FDR_noMSC.csv" = 
                        all.degs.noMSC[all.degs.noMSC$adj.P.Val <= 0.05, ])
all.DEG.files <- list("TMT_Differential_expression_results.csv" = 
                        all.degs,
                      "TMT_Differential_expression_results_max_5_percent_FDR.csv" = 
                        all.degs[all.degs$adj.P.Val <= 0.05, ],
                      "TMT_Differential_expression_results_noMSC.csv" = 
                        all.degs.noMSC,
                      "TMT_Differential_expression_results_max_5_percent_FDR_noMSC.csv" = 
                        all.degs.noMSC[all.degs.noMSC$adj.P.Val <= 0.05, ])
save_to_synapse(all.DEG.files, synapse_id)

#### panSEA with cell type as a factor ####
setwd(base.path)
if (file.exists("gmt_BeatAML_drug_MOA.rds")) {
  gmt.drug <- readRDS("gmt_BeatAML_drug_MOA.rds")
} else {
  gmt.drug <- DMEA::as_gmt(moa.BeatAML, sep = ", ")
  saveRDS(gmt.drug, "gmt_BeatAML_drug_MOA.rds")
}

synapse_id <- "syn68888753"
all.degs <- data.frame()
all.degs.noMSC <- data.frame()
all.degs.noMSC.bead <- data.frame()
contrasts <- c("CD14", "CD34", #"Aza", "Ven", "Aza.Ven", # Aza + Ven has NA- need to debug why comparison isn't sensitive vs. resistant
               "Ven", "Sort Type", "MSC")
contrasts.noMSC <- c("CD14", "Sort Type", #"Aza", "Ven", "Aza.Ven"
                     "Ven")
contrasts <- c("CD14", "CD34", "Aza", "Ven", "Aza.Ven", 
               "Sort Type", "MSC")
contrasts.noMSC <- c("CD14", "Sort Type", "Aza", "Ven", "Aza.Ven")
contrasts.noMSC.sort <- c("CD14", "Aza", "Ven", "Aza.Ven")
#contrasts.noMSC <- c("CD14", "Sort Type")
#BeatAML.data <- load_not_norm_BeatAML_for_DMEA2()
#sorted.patients <- unique(mCombo$patient)

#gmt1 <- readRDS("gmt1_more.rds")
gmt1 <- msigdbr::msigdbr(collection="H")
gmt1 <- DMEA::as_gmt(as.data.frame(gmt1),"gene_symbol","gs_name")
gmt1 <- list("Hallmark" = gmt1)
saveRDS(gmt1,"gmt1_more.rds")
setwd("analysis")
dir.create("combined24-27")
setwd("combined24-27")
dir.create("using_cellType-sortType-patient_factors")
setwd("using_cellType-sortType-patient_factors")
# gmt1 <- get_gmt1_v2()
# gmt1[[11]] <- NULL
# gmt[[13]] <- NULL
# names(gmt1) <- c("TFT_GTRD", "MIR_MIRDB", "GO_BP", "GO_CC", "GO_MF", 
#                  "Oncogenic_signatures", "BioCarta", "KEGG", "PID", "Reactome", 
#                  "Hallmark", "Positional")
# temp.gmt <- msigdbr::msigdbr(species = "Homo sapiens", category="C2", subcategory="CP:WIKIPATHWAYS")
# temp.gmt <- DMEA::as_gmt(as.data.frame(temp.gmt), element.names = "gene_symbol", set.names = "gs_name", descriptions = "gs_description")
# gmt1[["WikiPathways"]] <- temp.gmt
# gmt1 <- gmt1[c(1:12,14)]
# saveRDS(gmt1, "gmt1_more.rds")

# gmt.names <- c("Hallmark","PID", "Oncogenic_signatures", "KEGG")
# gmt1 <- gmt1[gmt.names]
# saveRDS(gmt1, "gmt1_more.rds")
#devtools::install_github("cstawitz/roomba")
# library(roomba)
# packageurl <- "http://cran.r-project.org/src/contrib/Archive/ggplot2/ggplot2_3.4.4.tar.gz"
# install.packages(packageurl, repos=NULL, type="source")
# library(ggplot2)
base.path <- getwd()
dia.wo.out$global <- as.data.frame(dia.wo.out$global)
dia.wo.out$global$Gene <- rownames(dia.wo.out$global)
unmatched <- rownames(dia.wo.out$meta)[!(rownames(dia.wo.out$meta) %in% colnames(dia.wo.out$global))] # none

method.data <- list("DIA_2batches_noOutliers"=dia.wo.out, "TMT_noOutliers"=tmt.wo.out)
synapser::synLogin()
for (k in 1:length(method.data)) {
  setwd(base.path)
  meta.df <- method.data[[k]]$meta
  if (names(method.data[k]) != "DIA") {
    meta.df$id <- meta.df$DIA_id
  }
  meta.df$Patient <- meta.df$patient
  meta.df$`Sample Type` <- meta.df$cellType
  meta.df$patientID <- make.names(meta.df$patient)
  
  sorted.patients <- unique(meta.df$patient)
  BeatAML.data <- load_not_norm_BeatAML_for_DMEA3(exclude.samples = sorted.patients)
  
  # run contrast combos
  if (grepl("DIA", names(method.data)[k])) {
    omics <- list("global" = method.data[[k]]$global)
    feature.names <- "Gene"
    temp.expr <- list(BeatAML.data$global)
  } else {
    omics <- list("global" = method.data[[k]]$global)
    feature.names <- "Gene"
    temp.expr <- list(BeatAML.data$global)
    
    # # don't use phospho since there were low cell counts
    # omics <- list("global" = method.data[[k]]$global,
    #               "phospho" = method.data[[k]]$phospho)
    # feature.names <- c("Gene", "SUB_SITE")
    # temp.expr <- list(BeatAML.data$global, BeatAML.data$phospho)
  }
  names(temp.expr) <- names(omics)
  
  # # remove MSCs
  # method.path.noMSC <- file.path(base.path, paste0(names(method.data)[k],"_noMSC"))
  # dir.create(paste0(names(method.data)[k],"_noMSC"))
  # setwd(paste0(names(method.data)[k],"_noMSC"))
  # methodFolder.noMSC <- 
  #   synapser::synStore(synapser::Folder(paste0(names(method.data)[k],"_noMSC"),
  #                                       parent = synapse_id))
  # meta.df.noMSC <- meta.df[meta.df$MSC == "Non_MSC",]
  # panSEA2_combos2(contrasts.noMSC, contrast2=c("cellType","Sort Type","patientID"), 
  #                 meta.df = meta.df.noMSC, omics = omics, expr = temp.expr,
  #                 gmt.drug = gmt.drug, drug.sens = BeatAML.data$drug,
  #                 base.path = base.path, temp.path = method.path.noMSC,
  #                 synapse_id = methodFolder.noMSC, n.net=0, DMEA=FALSE)
  # 
  # # get compiled DEGs for analyses without MSCs
  # methodDEGs.noMSC <- as.list(synapser::synGetChildren(methodFolder.noMSC, list("file"), sortBy = 'NAME'))
  # if (length(methodDEGs.noMSC) > 0) {
  #   if (methodDEGs.noMSC[[1]]$name == "Differential_expression_results.csv") {
  #     methodFile <- synapser::synGet(methodDEGs.noMSC[[1]]$id)
  #     methodDEG <- read.csv(methodFile$path)
  #     methodDEG$method <- names(method.data)[k]
  #     all.degs.noMSC <- rbind(all.degs.noMSC, methodDEG)
  #   }
  # }
  # 
  # filter for sorting method
  for (temp.sort in unique(meta.df$`Sort Type`)) {
    method.path.noMSC.bead <- file.path(base.path, paste0(names(method.data)[k],"_noMSC_",temp.sort))
    setwd(base.path)
    dir.create(paste0(names(method.data)[k],"_noMSC_",temp.sort))
    setwd(paste0(names(method.data)[k],"_noMSC_",temp.sort))
    methodFolder.noMSC.bead <- 
      synapser::synStore(synapser::Folder(paste0(names(method.data)[k],"_noMSC_",temp.sort),
                                          parent = synapse_id))
    meta.df.noMSC.bead <- meta.df[meta.df$MSC == "Non_MSC" & meta.df$`Sort Type`==temp.sort,]
    panSEA2_combos2(contrasts.noMSC.sort, contrast2=c("cellType","patientID"), 
                    meta.df = meta.df.noMSC.bead, omics = omics, expr = temp.expr,
                    gmt.drug = gmt.drug, drug.sens = BeatAML.data$drug,
                    base.path = base.path, temp.path = method.path.noMSC.bead,
                    synapse_id = methodFolder.noMSC.bead, n.net=0, DMEA=FALSE)
    # re-run line above to pick up where you left off if there is an error
    
    # get compiled DEGs
    methodDEGs.noMSC.bead <- as.list(synapser::synGetChildren(methodFolder.noMSC.bead, list("file"), sortBy = 'NAME'))
    if (length(methodDEGs.noMSC.bead) > 0) {
      if (methodDEGs.noMSC.bead[[1]]$name == "Differential_expression_results.csv") {
        methodFile <- synapser::synGet(methodDEGs.noMSC.bead[[1]]$id)
        methodDEG <- read.csv(methodFile$path)
        methodDEG$method <- names(method.data)[k]
        methodDEG$sortType <- temp.sort
        all.degs.noMSC.bead <- rbind(all.degs.noMSC.bead, methodDEG)
      }
    }
  }
  
}
setwd(base.path)
all.DEG.files <- list("Differential_expression_results_noMSC.csv" = 
                        all.degs.noMSC,
                      "Differential_expression_results_max_5_percent_FDR_noMSC.csv" = 
                        all.degs.noMSC[all.degs.noMSC$adj.P.Val <= 0.05, ])
all.DEG.files <- list("Differential_expression_results_noMSC_sortFiltered.csv" = 
                        all.degs.noMSC.bead,
                      "Differential_expression_results_max_5_percent_FDR_noMSC_sortFiltered.csv" = 
                        all.degs.noMSC.bead[all.degs.noMSC.bead$adj.P.Val <= 0.05, ])
save_to_synapse(all.DEG.files, synapse_id)


#### look at STRING network for DIA bead: CD14+ vs. CD34+ ####
library(PCSF)
data("STRINGv12")
de <- read.csv("analysis/combined24-27/using_cellType-sortType-patient_factors/DIA_2batches_noOutliers_noMSC_Bead/no_filter/CD14_Pos_vs_Neg/global/Differential_expression/Differential_expression_results.csv")
de <- de[de$adj.P.Val<=0.05,] # 2597 gene symbols
pos.de <- de[de$Log2FC>0,] # 1136
neg.de <- de[de$Log2FC<0,] # 1461
hist(pos.de$Log2FC) # maybe look at log2FC>1
hist(neg.de$Log2FC) # maybe look at log2FC<-1

pos.de1 <- de[de$Log2FC > 1,] # 478
neg.de1 <- de[de$Log2FC < -1,] # 455

pos.de2 <- de[de$Log2FC > 2,] # 117
neg.de2 <- de[de$Log2FC < -2,] # 24

pos.de3 <- de[de$Log2FC > 3,] # 30
neg.de3 <- de[de$Log2FC < -3,] # 1: CD34

pos.de4 <- de[de$Log2FC > 4,] # 7
neg.de4 <- de[de$Log2FC < -4,] # 0

pos.hall.edges <- STRINGv12[STRINGv12$from %in% pos.de3$Gene |
                          STRINGv12$to %in% pos.de3$Gene,]
pos.hall.vert <- data.frame(unique(c(pos.hall.edges$from, pos.hall.edges$to)), type="Inferred", absLog2FC=NA)
colnames(pos.hall.vert)[1] <- "name"
pos.hall.vert[pos.hall.vert$name %in% pos.de3$Gene,]$type <- "Input"
used.input <- unique(pos.hall.vert[pos.hall.vert$type=="Input",]$name) # 115
for (i in used.input) {
  pos.hall.vert[pos.hall.vert$name == i,]$absLog2FC <- abs(pos.de3[pos.de3$Gene==i,]$Log2FC)
}
topGraph <- igraph::graph_from_data_frame(pos.hall.edges, directed=FALSE, vertices = pos.hall.vert)
centrality.df <- data.frame(name = V(topGraph)$name,
                            degree = igraph::degree(topGraph, mode="all"),
                            closeness = igraph::closeness(topGraph, mode="all"),
                            betweenness = igraph::betweenness(topGraph, directed = FALSE),
                            eigen_centrality = igraph::eigen_centrality(topGraph, directed = FALSE)$vector,
                            hub_score = igraph::hub_score(topGraph)$vector,
                            authority_score = igraph::authority_score(topGraph)$vector)
RCy3::createNetworkFromIgraph(topGraph, title="pos.hall_exp24")

library(plyr);library(dplyr)
top5pos <- de %>% slice_max(Log2FC, n=5)
top5pos.edges <- STRINGv12[STRINGv12$from %in% top5pos$Gene |
                          STRINGv12$to %in% top5pos$Gene,]
top5pos.vert <- data.frame(unique(c(top5pos.edges$from, top5pos.edges$to)), type="Inferred", absLog2FC=NA)
colnames(top5pos.vert)[1] <- "name"
top5pos.vert[top5pos.vert$name %in% top5pos$Gene,]$type <- "Input"
used.input <- unique(top5pos.vert[top5pos.vert$type=="Input",]$name) # 5
for (i in used.input) {
  top5pos.vert[top5pos.vert$name == i,]$absLog2FC <- abs(top5pos[top5pos$Gene==i,]$Log2FC)
}
topGraph <- igraph::graph_from_data_frame(top5pos.edges, directed=FALSE, vertices = top5pos.vert)
centrality.df <- data.frame(name = V(topGraph)$name,
                            degree = igraph::degree(topGraph, mode="all"),
                            closeness = igraph::closeness(topGraph, mode="all"),
                            betweenness = igraph::betweenness(topGraph, directed = FALSE),
                            eigen_centrality = igraph::eigen_centrality(topGraph, directed = FALSE)$vector,
                            hub_score = igraph::hub_score(topGraph)$vector,
                            authority_score = igraph::authority_score(topGraph)$vector)
RCy3::createNetworkFromIgraph(topGraph, title="top5pos_exp24")

# filter for what connects these 5 DEGs
top5pos.vert$directCon <- NA
for (i in used.input) {
  directConnections <- unique(c(top5pos.edges[top5pos.edges$from == i,]$to,
                                top5pos.edges[top5pos.edges$to == i,]$from))
  top5pos.vert[top5pos.vert$name == i,]$directCon <- paste0(directConnections, collapse=", ")
}
allDirCon <- na.omit(unique(unlist(strsplit(top5pos.vert$directCon, ", ")))) # 516

# top5pos.vert$inputCon <- NA
# for (i in 1:nrow(top5pos.vert)) {
#   j <- top5pos.vert$name[i]
#   if (!(j %in% used.input)) {
#     # find what inputs the protein is connected to
#     tempInputCon <- 
#   }
# }

nInputCon <- plyr::ddply(top5pos.edges[!(top5pos.edges$from %in% used.input),], 
                         .(from), summarize,
                         inputCon = paste0(unique(to[to %in% used.input]), collapse=", "),
                         n = length(unique(to[to %in% used.input])))
top5posCon.edges <- top5pos.edges[top5pos.edges$from %in% nInputCon[nInputCon$n>1,]$from,] # 163
top5posCon.vert <- top5pos.vert[top5pos.vert$name %in% c(top5posCon.edges$from, top5posCon.edges$to),] # 80
topGraph <- igraph::graph_from_data_frame(top5posCon.edges, directed=FALSE, vertices = top5posCon.vert)
centrality.df <- data.frame(name = V(topGraph)$name,
                            degree = igraph::degree(topGraph, mode="all"),
                            closeness = igraph::closeness(topGraph, mode="all"),
                            betweenness = igraph::betweenness(topGraph, directed = FALSE),
                            eigen_centrality = igraph::eigen_centrality(topGraph, directed = FALSE)$vector,
                            hub_score = igraph::hub_score(topGraph)$vector,
                            authority_score = igraph::authority_score(topGraph)$vector)
RCy3::createNetworkFromIgraph(topGraph, title="top5posCon_exp24")

hall <- msigdbr::msigdbr(collection="H")
hall.nfkb <- hall[grepl("NFKB",hall$gs_name),] # 200 genes
pos3.edges <- STRINGv12[STRINGv12$from %in% pos.de3$Gene |
                          STRINGv12$to %in% pos.de3$Gene,]
pos.hall.edges <- pos3.edges[pos3.edges$from %in% c(pos.de3$Gene, hall.nfkb$gene_symbol) &
                               pos3.edges$to %in% c(pos.de3$Gene, hall.nfkb$gene_symbol),] # 198
pos.hall.vert <- data.frame(unique(c(pos.hall.edges$from, pos.hall.edges$to)), type="Inferred", absLog2FC=NA) # 62
colnames(pos.hall.vert)[1] <- "name"
pos.hall.vert[pos.hall.vert$name %in% pos.de3$Gene,]$type <- "Input"
used.input <- unique(pos.hall.vert[pos.hall.vert$type=="Input",]$name) # 22 / 30
for (i in used.input) {
  pos.hall.vert[pos.hall.vert$name == i,]$absLog2FC <- abs(pos.de3[pos.de3$Gene==i,]$Log2FC)
}
pos.hall.vert$NFKB <- FALSE
pos.hall.vert[pos.hall.vert$name %in% hall.nfkb$gene_symbol,]$NFKB <- TRUE
topGraph <- igraph::graph_from_data_frame(pos.hall.edges, directed=FALSE, vertices = pos.hall.vert)
RCy3::createNetworkFromIgraph(topGraph, title="posHallNFKB3_exp24")

hall.nfkb <- hall[grepl("NFKB",hall$gs_name),] # 200 genes
hdac.genes <- unique(c(STRINGv12$from[startsWith(STRINGv12$from, "HDAC")],STRINGv12$to[startsWith(STRINGv12$to, "HDAC")])) # 11: HDAC1-11
pos.hall.edges <- STRINGv12[STRINGv12$from %in% c(pos.de3$Gene, hall.nfkb$gene_symbol, "MAPK14", hdac.genes) & # MAPK14: p38 MAPK
                               STRINGv12$to %in% c(pos.de3$Gene, hall.nfkb$gene_symbol, "MAPK14", hdac.genes),] # 1996
pos.hall.vert <- data.frame(unique(c(pos.hall.edges$from, pos.hall.edges$to)), type="Inferred", absLog2FC=NA) # 200
colnames(pos.hall.vert)[1] <- "name"
pos.hall.vert[pos.hall.vert$name %in% pos.de3$Gene,]$type <- "Input"
used.input <- unique(pos.hall.vert[pos.hall.vert$type=="Input",]$name) # 22 / 30
for (i in used.input) {
  pos.hall.vert[pos.hall.vert$name == i,]$absLog2FC <- abs(pos.de3[pos.de3$Gene==i,]$Log2FC)
}
pos.hall.vert$NFKB <- FALSE
pos.hall.vert[pos.hall.vert$name %in% hall.nfkb$gene_symbol,]$NFKB <- TRUE
pos.hall.vert$HDAC <- FALSE
pos.hall.vert[pos.hall.vert$name %in% hdac.genes,]$HDAC <- TRUE
pos.hall.vert$p38 <- FALSE
pos.hall.vert[pos.hall.vert$name == "MAPK14",]$p38 <- TRUE
topGraph <- igraph::graph_from_data_frame(pos.hall.edges, directed=FALSE, vertices = pos.hall.vert)
RCy3::createNetworkFromIgraph(topGraph, title="posHallNFKB3-HDAC-p38_exp24")

pos.hall.edges <- STRINGv12[STRINGv12$from %in% c(pos.de3$Gene, "MAPK14", hdac.genes) & # MAPK14: p38 MAPK
                              STRINGv12$to %in% c(pos.de3$Gene, "MAPK14", hdac.genes),] # 122
pos.hall.vert <- data.frame(unique(c(pos.hall.edges$from, pos.hall.edges$to)), type="Inferred", absLog2FC=NA) # 31
colnames(pos.hall.vert)[1] <- "name"
pos.hall.vert[pos.hall.vert$name %in% pos.de3$Gene,]$type <- "Input"
used.input <- unique(pos.hall.vert[pos.hall.vert$type=="Input",]$name) # 19 / 30
for (i in used.input) {
  pos.hall.vert[pos.hall.vert$name == i,]$absLog2FC <- abs(pos.de3[pos.de3$Gene==i,]$Log2FC)
}
pos.hall.vert$NFKB <- FALSE
pos.hall.vert[pos.hall.vert$name %in% hall.nfkb$gene_symbol,]$NFKB <- TRUE # none
pos.hall.vert$HDAC <- FALSE
pos.hall.vert[pos.hall.vert$name %in% hdac.genes,]$HDAC <- TRUE
pos.hall.vert$p38 <- FALSE
pos.hall.vert[pos.hall.vert$name == "MAPK14",]$p38 <- TRUE
topGraph <- igraph::graph_from_data_frame(pos.hall.edges, directed=FALSE, vertices = pos.hall.vert)
RCy3::createNetworkFromIgraph(topGraph, title="posHDAC-p38_exp24")

pos.hall.edges <- STRINGv12[STRINGv12$from %in% c(hall.nfkb$gene_symbol, "MAPK14", hdac.genes) & # MAPK14: p38 MAPK
                              STRINGv12$to %in% c(hall.nfkb$gene_symbol, "MAPK14", hdac.genes),] # 1794
pos.hall.vert <- data.frame(unique(c(pos.hall.edges$from, pos.hall.edges$to)), type="Inferred", absLog2FC=NA) # 177
colnames(pos.hall.vert)[1] <- "name"
pos.hall.vert[pos.hall.vert$name %in% pos.de$Gene,]$type <- "Input"
used.input <- unique(pos.hall.vert[pos.hall.vert$type=="Input",]$name) # 24 / 1136
any(used.input %in% pos.de3$Gene) # FALSE
for (i in used.input) {
  pos.hall.vert[pos.hall.vert$name == i,]$absLog2FC <- abs(pos.de[pos.de$Gene==i,]$Log2FC)
}
pos.hall.vert$NFKB <- FALSE
pos.hall.vert[pos.hall.vert$name %in% hall.nfkb$gene_symbol,]$NFKB <- TRUE # none
pos.hall.vert$HDAC <- FALSE
pos.hall.vert[pos.hall.vert$name %in% hdac.genes,]$HDAC <- TRUE
pos.hall.vert$p38 <- FALSE
pos.hall.vert[pos.hall.vert$name == "MAPK14",]$p38 <- TRUE
pos.hall.vert$Pathway <- ""
pos.hall.vert[pos.hall.vert$name %in% hall.nfkb$gene_symbol,]$Pathway <- "NFKB"
pos.hall.vert[pos.hall.vert$name %in% hdac.genes,]$Pathway <- "HDAC"
pos.hall.vert[pos.hall.vert$name == "MAPK14",]$Pathway <- "p38 MAPK"
topGraph <- igraph::graph_from_data_frame(pos.hall.edges, directed=FALSE, vertices = pos.hall.vert)
RCy3::createNetworkFromIgraph(topGraph, title="posHallNFKB-HDAC-p38_exp24")

pos.hall.edges <- STRINGv12[STRINGv12$from %in% c("MAPK14", hdac.genes) & # MAPK14: p38 MAPK
                              STRINGv12$to %in% c("MAPK14", hdac.genes),] # 66
pos.hall.vert <- data.frame(unique(c(pos.hall.edges$from, pos.hall.edges$to)), type="Inferred", absLog2FC=NA) # 12
colnames(pos.hall.vert)[1] <- "name"
pos.hall.vert[pos.hall.vert$name %in% pos.de$Gene,]$type <- "Input"
used.input <- unique(pos.hall.vert[pos.hall.vert$type=="Input",]$name) # 24 / 1136
any(used.input %in% pos.de3$Gene) # FALSE
for (i in used.input) {
  pos.hall.vert[pos.hall.vert$name == i,]$absLog2FC <- abs(pos.de[pos.de$Gene==i,]$Log2FC)
}
pos.hall.vert$NFKB <- FALSE
pos.hall.vert[pos.hall.vert$name %in% hall.nfkb$gene_symbol,]$NFKB <- TRUE # none
pos.hall.vert$HDAC <- FALSE
pos.hall.vert[pos.hall.vert$name %in% hdac.genes,]$HDAC <- TRUE
pos.hall.vert$p38 <- FALSE
pos.hall.vert[pos.hall.vert$name == "MAPK14",]$p38 <- TRUE
pos.hall.vert$Pathway <- ""
pos.hall.vert[pos.hall.vert$name %in% hall.nfkb$gene_symbol,]$Pathway <- "NFKB"
pos.hall.vert[pos.hall.vert$name %in% hdac.genes,]$Pathway <- "HDAC"
pos.hall.vert[pos.hall.vert$name == "MAPK14",]$Pathway <- "p38 MAPK"
topGraph <- igraph::graph_from_data_frame(pos.hall.edges, directed=FALSE, vertices = pos.hall.vert)
RCy3::createNetworkFromIgraph(topGraph, title="posHDAC-p38_stricter_exp24")

# only HDAC4 was upregulated (in addition to MAPK14) - were the other hdacs quantified?
quant.hdac <- hdac.genes[hdac.genes %in% de$Gene] # HDAC2 and HDAC4 were detected

# maybe just look at things upregulated
pos.nfkb.p38.hdac <- pos.de[pos.de$Gene %in% c(hall.nfkb$gene_symbol, "MAPK14", hdac.genes),] # 28 / 1136
pos.hall.edges <- STRINGv12[STRINGv12$from %in% pos.nfkb.p38.hdac$Gene & # MAPK14: p38 MAPK
                              STRINGv12$to %in% pos.nfkb.p38.hdac$Gene,] # 54
pos.hall.vert <- data.frame(unique(c(pos.hall.edges$from, pos.hall.edges$to)), type="Input", absLog2FC=NA) # 15
colnames(pos.hall.vert)[1] <- "name"
used.input <- unique(pos.hall.vert[pos.hall.vert$type=="Input",]$name) # 15 / 28
any(used.input %in% pos.de3$Gene) # FALSE
for (i in used.input) {
  pos.hall.vert[pos.hall.vert$name == i,]$absLog2FC <- abs(pos.de[pos.de$Gene==i,]$Log2FC)
}
pos.hall.vert$NFKB <- FALSE
pos.hall.vert[pos.hall.vert$name %in% hall.nfkb$gene_symbol,]$NFKB <- TRUE # none
pos.hall.vert$HDAC <- FALSE
pos.hall.vert[pos.hall.vert$name %in% hdac.genes,]$HDAC <- TRUE
pos.hall.vert$p38 <- FALSE
pos.hall.vert[pos.hall.vert$name == "MAPK14",]$p38 <- TRUE
pos.hall.vert$Pathway <- ""
pos.hall.vert[pos.hall.vert$name %in% hall.nfkb$gene_symbol,]$Pathway <- "NFKB"
pos.hall.vert[pos.hall.vert$name %in% hdac.genes,]$Pathway <- "HDAC"
pos.hall.vert[pos.hall.vert$name == "MAPK14",]$Pathway <- "p38 MAPK"
topGraph <- igraph::graph_from_data_frame(pos.hall.edges, directed=FALSE, vertices = pos.hall.vert)
RCy3::createNetworkFromIgraph(topGraph, title="posDEHallNFKB-HDAC-p38_exp24")

pos.p38.hdac <- pos.de[pos.de$Gene %in% c("MAPK14", hdac.genes),] # 2 / 1136
pos.hall.edges <- STRINGv12[STRINGv12$from %in% pos.p38.hdac$Gene | # MAPK14: p38 MAPK
                              STRINGv12$to %in% pos.p38.hdac$Gene,] # 0 if AND logic, 1388 if OR logic
# pos.hall.vert <- data.frame(unique(c(pos.hall.edges$from, pos.hall.edges$to)), type="Input", absLog2FC=NA) # 15
# colnames(pos.hall.vert)[1] <- "name"
# used.input <- unique(pos.hall.vert[pos.hall.vert$type=="Input",]$name) # 15 / 28
# any(used.input %in% pos.de3$Gene) # FALSE
# for (i in used.input) {
#   pos.hall.vert[pos.hall.vert$name == i,]$absLog2FC <- abs(pos.de[pos.de$Gene==i,]$Log2FC)
# }
# pos.hall.vert$NFKB <- FALSE
# pos.hall.vert[pos.hall.vert$name %in% hall.nfkb$gene_symbol,]$NFKB <- TRUE # none
# pos.hall.vert$HDAC <- FALSE
# pos.hall.vert[pos.hall.vert$name %in% hdac.genes,]$HDAC <- TRUE
# pos.hall.vert$p38 <- FALSE
# pos.hall.vert[pos.hall.vert$name == "MAPK14",]$p38 <- TRUE
# pos.hall.vert$Pathway <- ""
# pos.hall.vert[pos.hall.vert$name %in% hall.nfkb$gene_symbol,]$Pathway <- "NFKB"
# pos.hall.vert[pos.hall.vert$name %in% hdac.genes,]$Pathway <- "HDAC"
# pos.hall.vert[pos.hall.vert$name == "MAPK14",]$Pathway <- "p38 MAPK"
# topGraph <- igraph::graph_from_data_frame(pos.hall.edges, directed=FALSE, vertices = pos.hall.vert)
# RCy3::createNetworkFromIgraph(topGraph, title="posDEHDAC-p38_exp24")

# are NFKB proteins more likely to be diffexp than not?
diffexp <- read.csv("analysis/combined24-27/using_cellType-sortType-patient_factors/DIA_2batches_noOutliers_noMSC_Bead/no_filter/CD14_Pos_vs_Neg/global/Differential_expression/Differential_expression_results.csv")
quant.nfkb <- diffexp[diffexp$Gene %in% hall.nfkb$gene_symbol,] # 57
pos.nfkb <- pos.de[pos.de$Gene %in% quant.nfkb$Gene,] # 24
neg.nfkb <- neg.de[neg.de$Gene %in% quant.nfkb$Gene,] # 2

# is MAPK14 or HDAC4 more likely to interact with NFKB signaling proteins than others?

# before updating panSEA package
# [1] "Running CD14_Pos_vs_Neg with no filter"
# Loading required namespace: snow
# Loading required namespace: doSNOW
# Running ssGSEA using global data
# Running enrichment analysis...
# Error in `colnames<-`(`*tmp*`, value = `*vtmp*`) : 
#   attempt to set 'colnames' on an object with less than two dimensions
# In addition: Warning message:
#   In GSEA_custom(input, gmt, num.permutations, stat.type, min.per.set,  :
#                    Removing drug sets with less than 6 drugs observed in data set

# update panSEA & then debug why not 2 dimensions in input to GSEA?

method.data <- list("phosphoPeptide" = highCov.pep.wo.out)
# there are duplicates in "SUB_SITE_MOD" resulting in duplicate rownames
# so we need to switch out the feature names
dia.features <- highCov.pep.wo.out$DIA$SUB_SITE_MOD
saveRDS(dia.features, paste0("DIA_SUB_SITE_MOD_",Sys.Date(),".rds"))
tmt.features <- highCov.pep.wo.out$TMT$SUB_SITE_MOD
saveRDS(tmt.features, paste0("TMT_SUB_SITE_MOD_",Sys.Date(),".rds"))
highCov.pep.wo.out$DIA$SUB_SITE_MOD <- NULL
highCov.pep.wo.out$TMT$SUB_SITE_MOD <- NULL
highCov.pep.wo.out$DIA$DIA <- seq(1,nrow(highCov.pep.wo.out$DIA))
highCov.pep.wo.out$TMT$TMT <- seq(1,nrow(highCov.pep.wo.out$TMT))
highCov.pep.wo.out$DIA <- highCov.pep.wo.out$DIA[,c("DIA",colnames(highCov.pep.wo.out$DIA)[1:(ncol(highCov.pep.wo.out$DIA)-1)])]
highCov.pep.wo.out$TMT <- highCov.pep.wo.out$TMT[,c("TMT",colnames(highCov.pep.wo.out$TMT)[1:(ncol(highCov.pep.wo.out$TMT)-1)])]
highCov.pep.wo.out$DIA$DIA <- as.character(highCov.pep.wo.out$DIA$DIA)
highCov.pep.wo.out$TMT$TMT <- as.character(highCov.pep.wo.out$TMT$TMT)
method.data <- list("phosphoPeptide" = highCov.pep.wo.out)

for (k in 1:length(method.data)) {
  setwd(base.path)
  method.path <- file.path(base.path, names(method.data)[k])
  dir.create(names(method.data)[k])
  setwd(names(method.data)[k])
  methodFolder <- 
    synapser::synStore(synapser::Folder(names(method.data)[k],
                                        parent = synapse_id))
  
  meta.df <- method.data[[k]]$meta
  meta.df$Patient <- meta.df$patient
  
  # run contrast combos
  omics <- method.data[[k]][c("DIA","TMT")]
  feature.names <- c("DIA", "TMT")
  #temp.expr <- list(BeatAML.data$phospho, BeatAML.data$phospho)
  temp.expr <- list(data.frame(),data.frame())
  names(temp.expr) <- names(omics)
  # panSEA2_combos2(contrasts, meta.df = meta.df,
  #                 omics = omics,
  #                 expr = temp.expr,
  #                 gmt.drug = gmt.drug, drug.sens = BeatAML.data$drug,
  #                 base.path = base.path,
  #                 temp.path = method.path,
  #                 synapse_id = methodFolder)
  # # there is an error at Aza.Ven_Sensitive_vs_Resistant with MSC == MSC but can just re-run loop or from line 340 to continue
  # 
  # # get compiled DEGs
  # methodDEGs <- as.list(synapser::synGetChildren(methodFolder, list("file"), sortBy = 'NAME'))
  # if (length(methodDEGs) > 0) {
  #   if (methodDEGs[[1]]$name == "Differential_expression_results.csv") {
  #     methodFile <- synapser::synGet(methodDEGs[[1]]$id)
  #     methodDEG <- read.csv(methodFile$path)
  #     methodDEG$method <- names(method.data)[k]
  #     all.degs <- rbind(all.degs, methodDEG)
  #   }
  # }
  
  # also filter for non-MSC
  setwd(base.path)
  method.path.noMSC <- file.path(base.path, paste0(names(method.data)[k],"_noMSC"))
  dir.create(paste0(names(method.data)[k],"_noMSC"))
  setwd(paste0(names(method.data)[k],"_noMSC"))
  methodFolder.noMSC <- 
    synapser::synStore(synapser::Folder(paste0(names(method.data)[k],"_noMSC"),
                                        parent = synapse_id))
  meta.df.noMSC <- meta.df[meta.df$MSC == "Non_MSC",]
  panSEA2_combos2(contrasts.noMSC, meta.df = meta.df.noMSC, 
                  omics = omics, gmt.list2=c(),
                  expr = temp.expr,
                  gmt.drug = gmt.drug, drug.sens = BeatAML.data$drug,
                  base.path = base.path,
                  temp.path = method.path.noMSC,
                  synapse_id = methodFolder.noMSC, GSEA = FALSE, DMEA = FALSE)
  # TMT no MSC: (just ran again from line above after)
  # [1] "Running Aza_Sensitive_vs_Resistant with Sort Type == Flow"
  # Error in .ebayes(fit = fit, proportion = proportion, stdev.coef.lim = stdev.coef.lim,  : 
  #                    No residual degrees of freedom in linear model fits
  #                  In addition: There were 50 or more warnings (use warnings() to see the first 50)
  
  # get compiled DEGs for analyses without MSCs
  methodDEGs.noMSC <- as.list(synapser::synGetChildren(methodFolder.noMSC, list("file"), sortBy = 'NAME'))
  if (length(methodDEGs.noMSC) > 0) {
    if (methodDEGs.noMSC[[1]]$name == "Differential_expression_results.csv") {
      methodFile <- synapser::synGet(methodDEGs.noMSC[[1]]$id)
      methodDEG <- read.csv(methodFile$path)
      methodDEG$method <- names(method.data)[k]
      all.degs.noMSC <- rbind(all.degs.noMSC, methodDEG)
    }
  }
}
setwd(base.path)

# switch feature names back
DIA.degs.noMSC <- all.degs.noMSC[all.degs.noMSC$Feature_type == "DIA",]
TMT.degs.noMSC <- all.degs.noMSC[all.degs.noMSC$Feature_type == "TMT",]
replaced.degs.noMSC <- data.frame()
types <- c("DIA", "TMT")
for (i in types) {
  temp.degs <- all.degs.noMSC[all.degs.noMSC$Feature_type == i,]
  if (i == "DIA") {
    conversion <- dia.features
  } else {
    conversion <- tmt.features
  }
  dict <- data.frame(conversion)
  dict$id <- rownames(dict)
  write.csv(dict, paste0("Peptide_labels_", i, ".csv"), row.names = FALSE)
  
  temp.degs$Feature <- with(dict, conversion[match(temp.degs$Feature, id)])
  replaced.degs.noMSC <- rbind(replaced.degs.noMSC, temp.degs)
}

all.DEG.files <- list(
  # "Differential_expression_results.csv" = 
  #                       all.degs,
  #                     "Differential_expression_results_max_5_percent_FDR.csv" = 
  #                       all.degs[all.degs$adj.P.Val <= 0.05, ],
                      "Differential_expression_results_noMSC_replacedPeptideLabels.csv" = 
                        replaced.degs.noMSC,
                      "Differential_expression_results_max_5_percent_FDR_noMSC_replacedPeptideLabels.csv" = 
                        replaced.degs.noMSC[replaced.degs.noMSC$adj.P.Val <= 0.05, ])
save_to_synapse(all.DEG.files, methodFolder.noMSC)

setwd("analysis/phosphoPeptide_noMSC/")
replaced.degs.noMSC <- read.csv("Differential_expression_results_noMSC_replacedPeptideLabels.csv")

# does the diffexp correlate between DIA & TMT?
temp.degs <- replaced.degs.noMSC[is.na(replaced.degs.noMSC$Filter) &
                                   replaced.degs.noMSC$Contrast == "CD14_Pos_vs_Neg",]
compareSigs <- function(sigs, venn.names = names(sigs), 
                        fname = "Monocyte_vs_progenitor_signatures") {
  og.path <- getwd()
  library(plyr)
  dir.create(fname)
  # import signatures and generate venn data
  venn.data <- list()
  sig.venn.data <- list()
  for (i in names(sigs)) {
    venn.data[[i]] <- sigs[[i]]$Feature
    sig.venn.data[[i]] <- sigs[[i]][sigs[[i]]$adj.P.Val <= 0.05,]$Feature
  }
  
  # Venn diagram
  setwd(fname)
  my.venn <- ggvenn::ggvenn(venn.data)
  ggplot2::ggsave("venn_diagram.pdf",my.venn, 
                  width = 7, height = 7)
  
  my.venn <- ggvenn::ggvenn(sig.venn.data)
  ggplot2::ggsave("venn_diagram_significant.pdf",my.venn, 
                  width = 7, height = 7)
  
  # Correlation matrix
  ## combine signatures into one data frame
  # start with first signature
  sig.df <- sigs[[1]]
  sig.df <- sig.df[,1:2]
  colnames(sig.df)[2] <- names(sigs)[1]
  
  # merge in other signatures one at a time
  for (i in 2:length(sigs)) {
    temp.sig <- sigs[[i]]
    colnames(temp.sig)[2] <- names(sigs)[i]
    sig.df <- merge(sig.df, temp.sig[,1:2], by="Feature")
  }
  write.csv(sig.df, "signature_Log2FC.csv", row.names = FALSE)
  #rownames(sig.df) <- sig.df$Feature
  
  if (length(sigs) == 2) {
    comparison <- paste0(names(sigs), collapse = "_and_")
    corr.df <- data.frame(comparison)
    corr.df[,c("Pearson.est", "Pearson.p", "Spearman.est", "Spearman.p")] <- NA
    pearson <- stats::cor.test(sig.df[,2], sig.df[,3], method="pearson")
    corr.df$Pearson.est <- pearson$estimate
    corr.df$Pearson.p <- pearson$p.value
    spearman <- stats::cor.test(sig.df[,2], sig.df[,3], method="spearman")
    corr.df$Spearman.est <- spearman$estimate
    corr.df$Spearman.p <- spearman$p.value
    write.csv(corr.df, "correlation_2signatures.csv", row.names = FALSE)
  }
  
  # create correlation matrix
  corr.mat <- stats::cor(as.matrix(sig.df[,2:ncol(sig.df)]))
  write.csv(corr.mat, "correlations.csv")
  
  # plot correlation matrix
  corr.mat.plot <- ggcorrplot::ggcorrplot(corr.mat)
  ggplot2::ggsave("correlation_matrix.pdf", 
                  corr.mat.plot, width = 7, height = 7)
  
  #try another approach
  DEG.df <- data.table::rbindlist(sigs, use.names = TRUE, idcol = "Signature", fill = TRUE)
  DEG.df$minusLogP <- -log(DEG.df$P.Value, base = 10)
  DEG.df$minusLogFDR <- -log(DEG.df$adj.P.Val, base = 10)
  DEG.df$sig <- FALSE
  DEG.df[DEG.df$adj.P.Val < 0.05, ]$sig <- TRUE
  
  # bar plot
  DEG.df$Direction <- "Upregulated"
  if (nrow(DEG.df[!is.na(DEG.df$Log2FC) & 
                  DEG.df$Log2FC < 0,]) > 0) {
    DEG.df[!is.na(DEG.df$Log2FC) &
             DEG.df$Log2FC < 0,]$Direction <- "Downregulated"
  }
  DEG.df$Significance <- "Not Significant"
  if (any(DEG.df$sig)) {
    DEG.df[DEG.df$sig,]$Significance <- DEG.df[DEG.df$sig,]$Direction
  }
  library(ggplot2)
  bar.plot3log <- ggplot2::ggplot(DEG.df, aes(fill = Significance, x=forcats::fct_infreq(Signature))) + 
    geom_bar(position = "dodge", stat="count") + ggplot2::theme_classic() + 
    scale_y_continuous(trans='log10') +
    ggplot2::xlab("Signature Source") +
    ggplot2::ylab("Number of Quantified Features")
  ggplot2::ggsave("barPlot_logScale.pdf", bar.plot3log, width = 7, height = 7, device = "pdf")
  
  bar.plot4log <- ggplot2::ggplot(DEG.df[DEG.df$Significance != "Not Significant",], 
                                  aes(fill = Direction, x=forcats::fct_infreq(Signature))) + 
    geom_bar(position = "dodge", stat="count") + ggplot2::theme_classic() + 
    scale_y_continuous(trans='log10') +
    ggplot2::xlab("Signature Source") +
    ggplot2::ylab("Number of Differentially Expressed Features")
  ggplot2::ggsave("barPlot_significant_logScale.pdf", bar.plot4log, width = 7, height = 7, device = "pdf")
  
  mean.DEG.df <- plyr::ddply(DEG.df, .(Feature), summarize,
                             mean_Log2FC = mean(Log2FC),
                             sd_Log2FC = sd(Log2FC),
                             Fisher_p = metap::sumlog(na.omit(P.Value))$p,
                             types = paste0(Signature, collapse = ", "),
                             N_types = length(unique(Signature)),
                             N_sig = length(sig[sig]),
                             sig_types = paste0(Signature[sig], collapse = ", "))
  if (length(unique(mean.DEG.df$Fisher_p)) > 1) {
    mean.DEG.df$adj_Fisher_p <- 
      qvalue::qvalue(mean.DEG.df$Fisher_p, pi0=1)$qvalues
  } else {
    mean.DEG.df$adj_Fisher_p <- NA
  }
  write.csv(mean.DEG.df, "compiled_signatures.csv", row.names = FALSE)
  write.csv(DEG.df, "signatures.csv", row.names = FALSE)
  setwd(og.path)
}
setwd(base.path)
sig.list <- list("DIA" = temp.degs[temp.degs$Feature_type=="DIA",], 
                 "TMT" = temp.degs[temp.degs$Feature_type=="TMT",])
compareSigs(sig.list, fname = "CD14_vs_CD34_peptide")

temp.degs2 <- na.omit(replaced.degs.noMSC[replaced.degs.noMSC$Filter=="Sort Type_Flow" &
                                   replaced.degs.noMSC$Contrast == "CD14_Pos_vs_Neg",]) # 0

temp.degs2 <- na.omit(replaced.degs.noMSC[replaced.degs.noMSC$Filter=="Sort Type_Bead" &
                                    replaced.degs.noMSC$Contrast == "CD14_Pos_vs_Neg",])
sig.list <- list("DIA" = temp.degs2[temp.degs2$Feature_type=="DIA",], 
                 "TMT" = temp.degs2[temp.degs2$Feature_type=="TMT",])
compareSigs(sig.list, fname = "CD14_vs_CD34_peptide_beadSorted")
#Error in cor.test.default(sig.df[, 2], sig.df[, 3], method = "pearson") : 
#  not enough finite observations

# also try MIT kinase enrichment: https://kinase-library.mit.edu/ea?a=ps
setwd("analysis/phosphoPeptide_noMSC/")
replaced.degs.noMSC <- read.csv("Differential_expression_results_noMSC_replacedPeptideLabels.csv")

# does the diffexp correlate between DIA & TMT?
temp.degs <- replaced.degs.noMSC[is.na(replaced.degs.noMSC$Filter) &
                                   replaced.degs.noMSC$Contrast == "CD14_Pos_vs_Neg",]
sig.list <- list("DIA" = temp.degs[temp.degs$Feature_type=="DIA",], 
                 "TMT" = temp.degs[temp.degs$Feature_type=="TMT",])
temp.degs2 <- na.omit(replaced.degs.noMSC[replaced.degs.noMSC$Filter=="Sort Type_Bead" &
                                            replaced.degs.noMSC$Contrast == "CD14_Pos_vs_Neg",])
sig.list2 <- list("DIA" = temp.degs2[temp.degs2$Feature_type=="DIA",], 
                 "TMT" = temp.degs2[temp.degs2$Feature_type=="TMT",])
all.sigs <- list("CD14vsCD34" = sig.list,
                 "CD14vsCD34_beadSorted" = sig.list2)
for (j in 1:length(all.sigs)) {
  fname <- names(all.sigs)[j]
  sig.list <- all.sigs[[j]]
  for (i in 1:length(sig.list)) {
    sig.name <- names(sig.list)[i]
    temp.sig <- sig.list[[i]]
    
    # only extract peptide sequence itself for MIT kinase enrichment
    temp.sig$Feature <- sub(".*@","",temp.sig$Feature)
    
    # only keep peptides with asterisk denoting phosphorylation
    temp.sig <- temp.sig[grepl("[*]",temp.sig$Feature),]
    
    foreground <- as.data.frame(temp.sig[temp.sig$adj.P.Val <= 0.05,]$Feature)
    background <- as.data.frame(temp.sig$Feature) # background needs to include foreground
    write.table(foreground, paste0(sig.name,"_foreground_peptides_",fname,".txt"), 
                sep="\t", row.names = FALSE, col.names = FALSE)
    write.table(background, paste0(sig.name,"_background_peptides_",fname,".txt"), 
                sep="\t", row.names = FALSE, col.names = FALSE)
  } 
}
# Error with each DIA and TMT:
# An error occurred during enrichment analysis. Invalid format?
#   Error message: At least 95% of sites must have at least 5 amino acids before and after phosphorylation site

# try ksdb instead
# ksdb <- read.csv(paste0("https://raw.githubusercontent.com/PNNL-CompBio/",
#                         "panSEA/shiny-app/data/ksdb_20231101.csv"))
# organism <- "human"
# if (organism %in% unique(na.omit(ksdb$KIN_ORGANISM)) &
#     organism %in% unique(na.omit(ksdb$SUB_ORGANISM))) {
#   ksdb <- ksdb[ksdb$KIN_ORGANISM == organism &
#                  ksdb$SUB_ORGANISM == organism, ]
# }
# # looks like the phospho sites are indicated by lower case letters instead of asterisk
# # need to replace s with S* and so on
# ksdb$sequence <- ksdb$SITE_...7_AA
# ksdb$phosphoAA <- stringr::str_extract_all(ksdb$sequence, "[[:lower:]]+")
# # extract all lowercase letters
# phosphoAA <- unique(unlist(stringr::str_extract_all(ksdb$sequence, "[[:lower:]]+")))
# # split up doublets (e.g., ys) and so on into individual letters
# phosphoAA <- unique(unlist(stringr::str_split(phosphoAA, "")))
# # replace lowercase letters with matching uppercase letter and asterisk
# for (i in phosphoAA) {
#   replacement <- paste0(toupper(i),"*")
#   ksdb$sequence <- gsub(i,replacement, ksdb$sequence)
# }
# ksdb$SUB_SITE_MOD <- paste0(ksdb$SUB_GENE,"@",ksdb$sequence)
# gmt <- DMEA::as_gmt(ksdb, "SUB_SITE_MOD", "KINASE",
#                     descriptions = "KIN_ACC_ID")
# saveRDS(gmt,"gmt_kdsb_human_peptide.rds") # 391 kinases with 6+ substrates
gmtv2 <- DMEA::as_gmt(ksdb, "SUB_SITE_MOD", "KINASE", descriptions = "KIN_ACC_ID", min.per.set=5)
saveRDS(gmtv2,"gmt_kdsb_human_peptide_min5perset.rds")
gmt <- readRDS("gmt_kdsb_human_peptide.rds")
gmt2 <- list(gmtv2, gmtv2)
names(gmt2) <- rep("ksdb", 2)
n.net <- 5
# base.path <- "~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/analysis"
# setwd(base.path)
# highCov.pep.wo.out <- readRDS("Phospho_peptide_DIA_TMT_overlappingSamples_min75PercentCoverage_noGlobalOutliers.rds")
# cc.df <- highCov.pep.wo.out$meta
synapse_id <- "syn64501548"
synapser::synLogin()
methodFolder.noMSC <- 
  synapser::synStore(synapser::Folder("phosphoPeptide_noMSC",
                                      parent = synapse_id))
library(plyr);library(dplyr)
all.ksea.files <- list()
for (j in 2:length(all.sigs)) {
  fname <- names(all.sigs)[j]
  sig.list <- all.sigs[[j]]
  
  # average over any duplicates
  for (i in 1:length(sig.list)) {
    temp.name <- names(sig.list)[i]
    temp.df <- sig.list[[temp.name]]
    temp.df <- plyr::ddply(temp.df, .(Feature), summarize,
                           Log2FC = mean(Log2FC, na.rm = TRUE))
    sig.list[[temp.name]] <- temp.df
  }
  
  gsea2 <- panSEA::mGSEA(sig.list, gmt2, types=names(sig.list),
                               feature.names=rep("Feature",2), min.per.set=5)
  # error: each feature must have only 1 rank.metric value (i.e., Log2FC)
  ksea.files <- list("KSEA_results_compiled.csv" = gsea2$compiled.results$mean.results,
                     "KSEA_results.csv" = gsea2$compiled.results$results,
                     "KSEA_venn_diagram.pdf" = gsea2$compiled.results$venn.diagram,
                     "KSEA_dot_plot.pdf" = gsea2$compiled.results$dot.plot,
                     "KSEA_correlations.csv" = gsea2$compiled.results$corr,
                     "KSEA_correlation_matrix.pdf" = gsea2$compiled.results$corr.matrix)
  for (i in 1:length(sig.list)) {
    temp.name <- names(sig.list)[i]
    PDFname <- paste0("KSEA_volcano_plot_",temp.name,".pdf")
    ksea.files[[PDFname]] <- gsea2$all.results[[temp.name]]$volcano.plot
  }
  all.ksea.files[[fname]] <- ksea.files
  
  #saveRDS(gsea2, paste0("KSEA_peptide_",fname,".rds"))
}
all.files <- list("KSEA" = all.ksea.files)
save_to_synapse(all.files, methodFolder.noMSC)
# error: insufficient coverage (need 2+ kinases with 6+ substrate sites each)

# 
# #### 4. Cell type predictions using weighted voting of CD14+ vs. CD34+ ####
# #### CD14 vs CD34: same as MSC_Non_MSC: CD14_Pos_vs_Neg
# #### Weighted voting: can CD14 vs CD34 signature rank data when not normalized across samples?
# ### prep raw DIA
# synapser::synLogin()
# globalFileDIA <- synapser::synGet("syn55234888") # Samantha's unprocessed version
# global.df.DIA <- read.table(
#   globalFileDIA$path, 
#   sep = "\t") # 8897 proteins, 48 samples
# 
# # require proteins to be quantified in at least half of samples
# global.df.DIA <- global.df.DIA[which(rowSums(is.na(global.df.DIA)) < ncol(global.df.DIA)/2),] # 6115 proteins, 48 samples
# 
# # require samples to have at least 75% of proteins quantified
# global.df.DIA <- global.df.DIA[ , which(colMeans(!is.na(global.df.DIA)) >= 0.75)] # 36 out of 48 samples are kept
# 
# # remove outliers
# global.df.DIA <- global.df.DIA[,colnames(global.df.DIA)[!(colnames(global.df.DIA) %in% c(outliers, outliers2a))]]
# 
# # log2-transform DIA data
# global.df.DIA <- log(global.df.DIA[,colnames(global.df.DIA)!="Gene"],2)
# 
# # subtract column (sample) medians
# global_sample_coef.DIA <- apply(global.df.DIA, 2, median, na.rm = T)
# global.df.DIA <- sweep(global.df.DIA, 2, global_sample_coef.DIA, FUN = '-')
# 
# # transform to prep for weighted voting
# global.df.DIA.trans <- t(global.df.DIA)
# 
# ### prep raw TMT
# globalFileTMT <- synapser::synGet("syn53493077") # unprocessed version
# global.df.TMT <- read.table(
#   globalFileTMT$path, 
#   sep = "\t") # 8897 proteins, 48 samples
# 
# # require proteins to be quantified in at least half of samples
# global.df.TMT <- global.df.TMT[which(rowSums(is.na(global.df.TMT)) < ncol(global.df.TMT)/2),] # 6115 proteins, 48 samples
# 
# # require samples to have at least 75% of proteins quantified
# global.df.TMT <- global.df.TMT[ , which(colMeans(!is.na(global.df.TMT)) >= 0.75)] # 36 out of 48 samples are kept
# 
# # remove outliers
# global.df.TMT <- global.df.TMT[,colnames(global.df.TMT)[!(colnames(global.df.TMT) %in% c(outlier.ids, outlier.ids2a))]]
# 
# # subtract column (sample) meTMTns
# global_sample_coef.TMT <- apply(global.df.TMT, 2, median, na.rm = T)
# global.df.TMT <- sweep(global.df.TMT, 2, global_sample_coef.TMT, FUN = '-')
# 
# # transform to prep for weighted voting
# global.df.TMT.trans <- t(global.df.TMT)
# 
# ### load signature of CD14+ vs. CD34+ cells
# synapser::synLogin()
# globalFileDIA <- synapser::synGet("syn59429685") # filtered for max FDR of 0.05
# global.sig.DIA <- read.csv(globalFileDIA$path) # 1842 proteins
# 
# globalFileTMT <- synapser::synGet("syn59438946") # filtered for max FDR of 0.05
# global.sig.TMT <- read.csv(globalFileTMT$path) # 703 proteins
# 
# ### add sample names
# global.df.DIA.trans <- as.data.frame(global.df.DIA.trans)
# global.df.DIA.trans$id <- rownames(global.df.DIA.trans)
# global.df.TMT.trans <- as.data.frame(global.df.TMT.trans)
# global.df.TMT.trans$id <- rownames(global.df.TMT.trans)
# global.df.DIA.trans <- global.df.DIA.trans[,c("id", global.sig.DIA$Gene[global.sig.DIA$Gene %in% colnames(global.df.DIA.trans)])]
# global.df.TMT.trans <- global.df.TMT.trans[,c("id", global.sig.TMT$Gene[global.sig.TMT$Gene %in% colnames(global.df.TMT.trans)])]
# 
# ### only consider proteins which have 100% coverage
# global.df.DIA.trans100 <- global.df.DIA.trans[,colSums(is.na(global.df.DIA.trans)) == 0] # 1120 proteins
# global.df.TMT.trans100 <- global.df.TMT.trans[,colSums(is.na(global.df.TMT.trans)) == 0] # 679 proteins
# global.df.DIA.trans100$id <- rownames(global.df.DIA.trans100)
# global.df.DIA.trans100 <- global.df.DIA.trans100[,c("id", colnames(global.df.DIA.trans100)[colnames(global.df.DIA.trans100)!="id"])]
# global.df.TMT.trans100$id <- rownames(global.df.TMT.trans100)
# global.df.TMT.trans100 <- global.df.TMT.trans100[,c("id", colnames(global.df.TMT.trans100)[colnames(global.df.TMT.trans100)!="id"])]
# 
# ### run WV
# DIA.WV <- DMEA::WV(global.df.DIA.trans100, global.sig.DIA, sample.names = "id")
# TMT.WV <- DMEA::WV(global.df.TMT.trans100, global.sig.TMT, sample.names = "id")
# 
# # merge with meta.df
# WV.results <- rbind(DIA.WV$scores, TMT.WV$scores)
# WV.df <- merge(WV.results, dia.tmt.wo.out2a$meta)
# WV.df$Pooled <- "MSC"
# WV.df[WV.df$Pooled_CD14_Pos, ]$Pooled <- "CD14+"
# WV.df[WV.df$Pooled_CD34_Pos, ]$Pooled <- "CD34+"
# WV.df$SampleType <- WV.df$Pooled
# WV.df[WV.df$Flow, ]$SampleType <- paste(WV.df[WV.df$Flow, ]$SampleType, "Flow")
# WV.df <- na.omit(WV.df)
# 
# # DIA p-values
# WV.p.CD14.DIA <- t.test(WV.df[WV.df$Pooled_CD14_Pos & WV.df$method == "DIA",]$WV, WV.df[!WV.df$Pooled_CD14_Pos & WV.df$method == "DIA",]$WV)$p.value
# WV.p.CD34.DIA <- t.test(WV.df[WV.df$Pooled_CD34_Pos & WV.df$method == "DIA",]$WV, WV.df[!WV.df$Pooled_CD34_Pos & WV.df$method == "DIA",]$WV)$p.value
# WV.p.MSC.DIA <- t.test(WV.df[WV.df$MSC_Flow & WV.df$method == "DIA",]$WV, WV.df[!WV.df$MSC_Flow & WV.df$method == "DIA",]$WV)$p.value
# 
# # TMT p-values
# WV.p.CD14.TMT <- t.test(WV.df[WV.df$Pooled_CD14_Pos & WV.df$method == "TMT",]$WV, WV.df[!WV.df$Pooled_CD14_Pos & WV.df$method == "TMT",]$WV)$p.value
# WV.p.CD34.TMT <- t.test(WV.df[WV.df$Pooled_CD34_Pos & WV.df$method == "TMT",]$WV, WV.df[!WV.df$Pooled_CD34_Pos & WV.df$method == "TMT",]$WV)$p.value
# WV.p.MSC.TMT <- t.test(WV.df[WV.df$MSC_Flow & WV.df$method == "TMT",]$WV, WV.df[!WV.df$MSC_Flow & WV.df$method == "TMT",]$WV)$p.value
# 
# ### violin plot: WV score vs. Sample Type
# setwd(base.path)
# library(ggplot2)
# marker.violin <- ggplot2::ggplot(WV.df, 
#                                  aes(fill = method, x=Pooled, y=WV)) + 
#   geom_violin(position=position_dodge(width=0.4), alpha=0.5) + 
#   geom_boxplot(width=0.1, position = position_dodge(width=0.4), alpha=0.5) + 
#   bg.theme3 + xlab("Sample Type") + ylab("CD14+ vs. CD34+ Score")
# ggsave(paste0("CD14_vs_CD34_WV_100PercentCoverage_by_pooled_sample_type_", Sys.Date(), ".pdf"), marker.violin)
# marker.violin <- ggplot2::ggplot(WV.df, 
#                                  aes(fill = method, x=SampleType, y=WV)) + 
#   geom_violin(position=position_dodge(width=0.4), alpha=0.5) + 
#   geom_boxplot(width=0.1, position = position_dodge(width=0.4), alpha=0.5) + 
#   bg.theme3 + xlab("Sample Type") + ylab("CD14+ vs. CD34+ Score")
# ggsave(paste0("CD14_vs_CD34_WV_100PercentCoverage_by_sample_type_", Sys.Date(), ".pdf"), marker.violin)
# 
# #### 5. Correlations between DIA & TMT ####
# corr.scatter <- function(df, rank.var, value, xlab, ylab, title, Pearson.p, Pearson.est, 
#                          se = TRUE, position.x = "min", position.y = "max", 
#                          shape = NULL, color = NULL, symmetrical = TRUE) {
#   # load themes for plots
#   ng.theme <- ggplot2::theme(
#     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#     panel.border = element_rect(fill = NA),
#     panel.background = element_blank(),
#     axis.line = element_line(colour = "black"),
#     axis.text.x = element_text(colour = "black"),
#     axis.text.y = element_text(colour = "black"),
#     axis.ticks.x = element_line(colour = "black"),
#     axis.ticks.y = element_line(colour = "black"),
#     #legend.title = element_blank(),
#     axis.title.y = element_text(size = 8, colour = "black")
#   )
#   
#   bg.theme <- ggplot2::theme(
#     legend.background = element_rect(), legend.position = "right",
#     legend.text = element_text(size = 14), 
#     legend.title = element_text(size=16),
#     legend.key = element_blank(),
#     axis.title.x = element_text(size = 20),
#     axis.text.x = element_text(size = 16),
#     axis.title.y = element_text(size = 20),
#     axis.text.y = element_text(size = 16),
#     plot.title = element_text(lineheight = .8, face = "bold", size = 36)
#   )
#   
#   # set plot parameters
#   min.x <- min(df[, c(rank.var)])
#   max.x <- max(df[, c(rank.var)])
#   mid.x <- 0.5 * (min.x + max.x)
#   min.y <- min(df[, c(value)])
#   max.y <- max(df[, c(value)])
#   mid.y <- 0.5 * (min.y + max.y)
#   if (position.x == "min") {
#     pos.x <- min.x
#   } else if (position.x == "mid") {
#     pos.x <- mid.x
#   } else if (position.x == "max") {
#     pos.x <- max.x
#   } else if (is.numeric(position.x)) {
#     pos.x <- position.x
#   }
#   if (position.y == "min") {
#     pos.y <- min.y
#   } else if (position.y == "mid") {
#     pos.y <- mid.y
#   } else if (position.y == "max") {
#     pos.y <- max.y
#   } else if (is.numeric(position.y)) {
#     pos.y <- position.y
#   }
#   
#   stats_pearson <- substitute(
#     r == est * "," ~ ~"p" ~ "=" ~ p,
#     list(
#       est = format(Pearson.est, digits = 3),
#       p = format(Pearson.p, digits = 3)
#     )
#   )
#   
#   scatter.plot <- ggplot2::ggplot(data = df,
#                                   aes_string(x = rank.var, y = value)) +
#     ggplot2::geom_point(aes_string(
#       shape = shape, color = color)) +
#     ggplot2::labs(x = xlab, y = ylab) +
#     ggplot2::ggtitle(title) +
#     ggplot2::geom_smooth(method = "lm", size = 1.5,
#                          linetype = "solid", color = "blue",
#                          se = se, na.rm = TRUE) +
#     ggplot2::geom_text(
#       x = pos.x, y = pos.y, vjust = "inward", hjust = "inward",
#       colour = "blue", parse = TRUE,
#       label = as.character(as.expression(stats_pearson)), size = 8
#     ) +
#     ng.theme +
#     bg.theme + theme_minimal()
#   
#   if (symmetrical) {
#     abs.max <- max(abs(c(min.x, max.x, min.y, max.y)))
#     
#     if (position.y == "min") {
#       pos.y <- -abs.max
#     } else if (position.y == "mid") {
#       pos.y <- 0
#     } else {pos.y <- abs.max}
#     
#     if (position.x == "min") {
#       pos.x <- -abs.max
#     } else if (position.x == "mid") {
#       pos.x <- 0
#     } else {pos.x <- abs.max}
#     
#     scatter.plot <- ggplot2::ggplot(data = df,
#                                     aes_string(x = rank.var, y = value)) +
#       ggplot2::geom_point(aes_string(
#         shape = shape, color = color)) +
#       ggplot2::labs(x = xlab, y = ylab) +
#       ggplot2::ggtitle(title) +
#       ggplot2::geom_smooth(method = "lm", size = 1.5,
#                            linetype = "solid", color = "blue",
#                            se = se, na.rm = TRUE) +
#       ggplot2::geom_text(
#         x = pos.x, y = pos.y, vjust = "inward", hjust = "inward",
#         colour = "blue", parse = TRUE,
#         label = as.character(as.expression(stats_pearson)), size = 8
#       ) +
#       ng.theme +
#       bg.theme + theme_minimal() + xlim(c(-abs.max, abs.max)) + 
#       ylim(c(-abs.max, abs.max)) + geom_hline(yintercept = 0) +
#       geom_vline(xintercept = 0) + theme(axis.title.x = element_text(size = 20),
#                                          axis.text.x = element_text(size = 16),
#                                          axis.title.y = element_text(size = 20),
#                                          axis.text.y = element_text(size = 16),
#                                          plot.title = element_text(lineheight = .8, face = "bold", size = 36))
#   }
#   
#   return(plot = scatter.plot)
# }
# 
# #### raw DIA vs. TMT paired for each sample in both
# # try giving TMT same sample names as DIA to correlate them
# rownames(meta.df) <- meta.df$id
# meta.wo.out2a <- meta.df[meta.df$id %in% colnames(global.df.TMT),]
# meta.wo.out2a <- meta.wo.out2a[colnames(global.df.TMT),]
# dia.names <- meta.wo.out2a$DIA_id
# colnames(global.df.TMT) <- dia.names
# 
# # make sure we are comparing the same genes in the same orders
# gene.names <- rownames(global.df.TMT)
# gene.names <- gene.names[gene.names %in% rownames(global.df.DIA)]
# matching.DIA <- global.df.DIA[gene.names,]
# matching.TMT <- global.df.TMT[gene.names,]
# 
# # make sure we have the same samples in the same order in DIA & TMT
# shared.dia.names <- dia.names[dia.names %in% colnames(matching.DIA)]
# matching.DIA <- matching.DIA[,shared.dia.names]
# matching.TMT <- matching.TMT[,shared.dia.names]
# DIA.TMT.correlations <- corrr::correlate(matching.DIA, matching.TMT, use="pairwise.complete.obs", diagonal = 1)
# rownames(DIA.TMT.correlations) <- DIA.TMT.correlations$term
# heatmap(as.matrix(dplyr::select_if(DIA.TMT.correlations, is.numeric)))
# pheatmap::pheatmap(as.matrix(dplyr::select_if(DIA.TMT.correlations, is.numeric)), cluster_rows = FALSE, cluster_cols = FALSE)
# 
# #### WV scores: DIA vs. TMT
# library(ggplot2)
# # prepare data frame
# DIA.WV.results <- merge(DIA.WV$scores, meta.df, by.x="id", by.y = "DIA_id", suffixes = c("_DIA", "_TMT"))
# WV.df <- merge(TMT.WV$scores, DIA.WV.results, by.x = "id", by.y = "id_TMT", suffixes = c("_TMT", "_DIA"))
# WV.df$'Sample' <- WV.df$SampleType
# WV.df$'Patient' <- WV.df$patient
# 
# # run correlation
# corr.DIA.TMT <- cor.test(WV.df$WV_DIA, WV.df$WV_TMT)
# p <- corr.DIA.TMT$p.value
# est <- as.numeric(corr.DIA.TMT$estimate)
# 
# # create scatter plot
# WV.DIA.TMT.plot <- corr.scatter(WV.df, "WV_DIA", "WV_TMT", "DIA", "TMT", 
#                                 "CD14+ vs. CD34+ Score", p, est, 
#                                 shape = "Sample", color = "Patient")
# ggsave("DIA_vs_TMT_CD14_vs_CD34_score_v3.pdf", WV.DIA.TMT.plot, height=7, width=7)
# 
# #### differentially expressed proteins based on adjusted p-values <= 0.05: DIA vs. TMT
# # prepare data frame
# global.sig <- merge(global.sig.DIA, global.sig.TMT, by="Gene", suffixes=c("_DIA", "_TMT")) # 525
# 
# # run correlation
# corr.DIA.TMT <- cor.test(global.sig$Log2FC_DIA, global.sig$Log2FC_TMT)
# p <- corr.DIA.TMT$p.value
# est <- as.numeric(corr.DIA.TMT$estimate)
# 
# # create scatter plot
# sig.DIA.TMT.plot <- corr.scatter(global.sig, "Log2FC_DIA", "Log2FC_TMT", 
#                                  "DIA", "TMT", "CD14+ vs. CD34+ Log2FC", p, est)
# ggsave("DIA_vs_TMT_CD14_vs_CD34_Log2FC_v3.pdf", sig.DIA.TMT.plot, height=7, width=7)
# 
# #### Venn diagrams: DIA vs. TMT
# ### differentially expressed proteins based on adjusted p-values <= 0.05
# venn.data <- list("TMT" = global.sig.TMT$Gene,
#                   "DIA" = global.sig.DIA$Gene)
# sig.venn <- ggvenn::ggvenn(venn.data, set_name_size = 8, text_size = 8)
# ggsave("Venn_diagram_DIA_vs_TMT_CD14_vs_CD34_differential_expression.pdf", sig.venn, height=8, width=8)
# 
# ### all proteins quantified after filters
# all.venn.data <- list("TMT" = rownames(global.df.TMT), 
#                       "DIA" = rownames(global.df.DIA))
# all.venn <- ggvenn::ggvenn(all.venn.data, set_name_size = 8, text_size = 8)
# ggsave("Venn_diagram_DIA_vs_TMT_proteins_quantified_after_filters.pdf", all.venn, height=8, width=8)
# 
# ##### 6. Plot cell markers and DEGS #####
# # plot markers
# # MSCs (Mesenchymal Stem Cells): should be high in CD73, CD90, CD105, CD106, CD146, STRO-1 and low in CD14, CD34, CD45, HLA-DR
# # AML monocytes should be high in CD4, CD11c, CD14, and CD64; source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5907644/#:~:text=Surface%20antigens%20that%20indicate%20leukemic,mostly%20expressed%20on%20mature%20monocytes.
# # Leukemia stem cells (LSCs):  CD34, CD38, CD123, TIM3, CD25, CD32 and CD96; source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5388677/
# # LSCs are emphasized in recent PTRC paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10950354/
# 
# markers <- unique(c("CD73", "CD90", "CD105", "CD106", "CD146", "STRO-1", "CD14", 
#                     "CD34", "CD45", "HLA-DR", "CD4", "CD11c", "CD14", "CD64", 
#                     "CD34", "CD38", "CD123", "TIM3", "CD25", "CD32", "CD96"))
# markers <- unique(c("NT5E", "THY1", "ENG", "VCAM1", "MCAM", "CD14", "CD34", 
#                     "PTPRC", "CD4", "ITGAX", "ITGAX", "FCGR1A", "CD38", "IL3RA", 
#                     "HAVCR2", "IL2RA", "ISG20", "FCGR2A", "FCGR2B", "CD96"))
# markers.from.Anupriya <- c("CD3", "HLA-DR", "CD1A", "CD4", "CD5", "ITGAL", 
#                            "ITGAM", "CD14", "FUT4", "ITGB2", "CD19", "IL2RA",
#                            "ISG20", "CD38", "NCAM1", "SELE", "SELP", "IL3RA",
#                            "CDH5", "FASLG", "CD9", "ITGB1", "CD44", "CD46",
#                            "CD47", "ITGA1", "ITGA2", "ITGA5", "CD58", "CD59",
#                            "ITGB3", "CD63", "NT5E", "CD81", "THY1", "SLC3A2",
#                            "SLC7A5", "BSG", "CD151", "CD200", "HLA-A", "HLA-B",
#                            "HLA-C", "ITGA3", "ITGAV", "FAS", "ENG", "CD13", 
#                            "ANPEP", "CD33")
# markers.from.Anupriya <- c(markers.from.Anupriya, c("CD31", "PECAM1", "CD90", "THY1", "CD105", "ENG"))
# all.markers <- unique(c(markers, markers.from.Anupriya))
# dia.tmt.markers <- dia.tmt.wo.out$global[dia.tmt.wo.out$global$Gene %in% markers,
#                                     c("Gene", colnames(dia.tmt.wo.out$global)[colnames(dia.tmt.wo.out$global) != "Gene"])]
# write.csv(dia.tmt.markers, paste0("Cell_markers_global_DIA_TMT_75percentCoverage_", Sys.Date(), ".csv"), row.names = FALSE)
# 
# ### violin plots of cell markers
# library(ggplot2)
# # load theme for plots
# bg.theme3 <- ggplot2::theme(
#   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#   panel.border = element_rect(fill = NA),
#   panel.background = element_blank(),
#   axis.line = element_line(colour = "black"), 
#   axis.title.x = element_text(size = 16, colour = "black"),
#   axis.title.y = element_text(size = 16, colour = "black"),
#   axis.text.x = element_text(colour = "black", size =16), 
#   axis.text.y = element_text(colour = "black", size = 16),
#   axis.ticks.x = element_line(colour = "black"), 
#   axis.ticks.y = element_line(colour = "black"),
#   legend.title = element_blank(), legend.background = element_rect(), 
#   legend.position = "top",
#   legend.text = element_text(size = 14), legend.key = element_blank(),
#   plot.title = element_text(lineheight = .8, face = "bold", size = 36)
# )
# 
# # prepare long data frame
# long.global <- reshape2::melt(dia.tmt.wo.out$global, variable.name = "id")
# long.global <- merge(long.global, dia.tmt.wo.out$meta, by = "id")
# long.global$Pooled <- gsub(" .*$", "", long.global$`Sample Type`)
# write.csv(long.global[long.global$Gene %in% markers,], paste0("Cell_markers_global_DIA_TMT_75percentCoverage_with_metadata_", Sys.Date(), ".csv"), row.names = FALSE)
# write.csv(long.global[long.global$Gene == "CD14",], paste0("CD14_global_DIA_TMT_75percentCoverage_with_metadata_", Sys.Date(), ".csv"), row.names = FALSE)
# write.csv(long.global[long.global$Gene == "CD34",], paste0("CD34_global_DIA_TMT_75percentCoverage_with_metadata_", Sys.Date(), ".csv"), row.names = FALSE)
# long.CD14 <- long.global[long.global$Gene == "CD14",]
# p.CD14.DIA <- t.test(long.CD14[long.CD14$method == "DIA" & long.CD14$CD14 == "Pos",]$value,
#                      long.CD14[long.CD14$method == "DIA" & long.CD14$CD14 == "Neg",]$value, alternative = "greater")$p.value # 2.69E-8
# p.CD14.TMT <- t.test(long.CD14[long.CD14$method == "TMT" & long.CD14$CD14 == "Pos",]$value,
#                      long.CD14[long.CD14$method == "TMT" & long.CD14$CD14 == "Neg",]$value, alternative = "greater")$p.value # 4.22E-6
# 
# long.CD34 <- long.global[long.global$Gene == "CD34",]
# p.CD34.DIA <- t.test(long.CD34[long.CD34$method == "DIA" & long.CD34$CD34 == "Pos",]$value,
#                      long.CD34[long.CD34$method == "DIA" & long.CD34$CD34 == "Neg",]$value, alternative = "greater")$p.value # 2.83E-6
# p.CD34.TMT <- t.test(long.CD34[long.CD34$method == "TMT" & long.CD34$CD34 == "Pos",]$value,
#                      long.CD34[long.CD34$method == "TMT" & long.CD34$CD34 == "Neg",]$value, alternative = "greater")$p.value # 2.51E-4
# 
# 
# setwd(base.path)
# dir.create("markers")
# setwd("markers")
# dir.create("20241003")
# setwd("20241003")
# for (i in 1:length(markers)) {
#   marker.df <- long.global[long.global$Gene == markers[i],]
#   if (nrow(marker.df) > 0) {
#     marker.violin <- ggplot2::ggplot(marker.df, 
#                                      aes(fill = method, x=Pooled, y=value)) + 
#       geom_violin(position=position_dodge(width=0.4), alpha=0.5) + 
#       geom_boxplot(width=0.1, position = position_dodge(width=0.4), alpha=0.5) + 
#       bg.theme3 + xlab("Sample Type") + ylab("Normalized Protein Expression")
#     ggsave(paste0(markers[i],"_75PercentCoverage_by_pooled_sample_type_", Sys.Date(), ".pdf"), marker.violin,
#            height = 7, width = 7)
#     marker.violin <- ggplot2::ggplot(marker.df, 
#                                      aes(fill = method, x=`Sample Type`, y=value)) + 
#       geom_violin(position=position_dodge(width=0.4), alpha=0.5) + 
#       geom_boxplot(width=0.1, position = position_dodge(width=0.4), alpha=0.5) + 
#       bg.theme3 + xlab("Sample Type") + ylab("Normalized Protein Expression")
#     ggsave(paste0(markers[i],"_75PercentCoverage_by_sample_type_", Sys.Date(), ".pdf"), marker.violin,
#            height = 7, width = 7)
#   }
# }
# 
# ### pathway data for heatmaps
# # repeat for KRAS pathway hits based on Jeff's paper
# jeff.markers <- c("KRAS", "PTPN11", "BCLXL", "MCL1", "CD40", "CD14", "CLEC7A", 
#              "TRAF2", "IRAK1", "NFKB", "BCL2A1", "BCL2", "BAK", "BAX")
# 
# # Myc targets
# msigdb.info <- msigdbr::msigdbr("Homo sapiens", "H")
# myc.targets <- unique(msigdb.info[grepl("MYC_TARGETS", msigdb.info$gs_name, 
#                                         ignore.case = TRUE), ]$gene_symbol) # 240
# myc.targetsv1 <- unique(msigdb.info[grepl("MYC_TARGETS_V1", msigdb.info$gs_name, 
#                                           ignore.case = TRUE), ]$gene_symbol) # 200
# myc.targetsv2 <- unique(msigdb.info[grepl("MYC_TARGETS_V2", msigdb.info$gs_name, 
#                                           ignore.case = TRUE), ]$gene_symbol) # 58
# 
# # JAK/STAT
# jak.stat <- unique(msigdb.info[grepl("JAK_STAT", msigdb.info$gs_name, 
#                                      ignore.case = TRUE), ]$gene_symbol) 
# 
# # KRAS
# kras <- unique(msigdb.info[grepl("KRAS", msigdb.info$gs_name, 
#                                  ignore.case = TRUE), ]$gene_symbol) 
# 
# # TGF beta
# msigdb.info <- msigdbr::msigdbr("Homo sapiens", "C2", "CP:KEGG")
# tgf <- unique(msigdb.info[grepl("TGF_BETA", msigdb.info$gs_name, 
#                                 ignore.case = TRUE), ]$gene_symbol) 
# 
# genesets <- list('Myc_targets' = myc.targets,
#                  'Myc_targets_v1' = myc.targetsv1,
#                  'Myc_targets_v2' = myc.targetsv2,
#                  "JAK-STAT" = jak.stat,
#                  "KRAS" = kras,
#                  "TGF_Beta" = tgf,
#                  "Jeff" = jeff.markers)
# 
# # prep annotations for heatmaps
# meta.df <- dia.tmt.wo.out$meta
# cc.df <- meta.df[,c("Sort Type", "Sample Type", "Patient")]
# 
# # create heatmaps
# DIA.heatmaps <- make_heatmaps(dia.wo.out$global, cc.df, top.gmt = genesets, fontsize=6)
# TMT.heatmaps <- make_heatmaps(tmt.wo.out$global, cc.df, top.gmt = genesets, fontsize=6)
# heatmaps <- list("Global_DIA" = DIA.heatmaps,
#                  "Global_TMT" = TMT.heatmaps)
# 
# # save heatmaps
# dir.create("heatmaps_of_interest")
# setwd("heatmaps_of_interest")
# dir.create("fontsize_6")
# setwd("fontsize_6")
# save_to_synapse(heatmaps)
# 
# # instead just look at individual data points
# dir.create("markers")
# setwd("markers")
# for (i in 1:length(markers)) {
#   marker.df <- na.omit(long.global[long.global$Gene == markers[i],])
#   if (nrow(marker.df) > 0) {
#     marker.violin <- ggplot2::ggplot(marker.df, 
#                                      aes(color = patient, x=Pooled, y=value, shape = method)) + 
#       geom_point(position=position_dodge(width=0.4)) +  
#       bg.theme3 + xlab("Sample Type") + ylab("Normalized Protein Expression")
#     ggsave(paste0(markers[i],"_75PercentCoverage_by_pooled_sample_type_dotplot_", Sys.Date(), ".pdf"), marker.violin, width = 11)
#     marker.violin <- ggplot2::ggplot(marker.df, 
#                                      aes(color = patient, x=`Sample Type`, y=value, shape = method)) + 
#       geom_point(position=position_dodge(width=0.4)) +  
#       bg.theme3 + xlab("Sample Type") + ylab("Normalized Protein Expression")
#     ggsave(paste0(markers[i],"_75PercentCoverage_by_sample_type_dotplot_", Sys.Date(), ".pdf"), marker.violin, width = 11)
#   }
# }
# 
