# how should we define our sorted proteomics signature for Monocyte?
# assuming trends with WV deconvolution will hold true with other deconvolution methods
# 1- what is the best filter for accurately scoring monocytes (based on t-test, correct category per sample)
# 2- what is the best filter for predicting drug sensitivity (based on Ven, Aza+Ven AUC)

library(synapser)
library(DMEA)
library(tidyr)
library(tibble)
library(reshape2)
library(plyr)
library(dplyr)
library(ggplot2)
synapser::synLogin()
source("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/helperScripts/circBar.R")
source("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/helperScripts/panSEA_helper_20240913.R")

evalOneMonoSig <- function(temp.sig, BeatAML, type="rna", gmt) {
  # calculate WV scores
  temp.sig <- na.omit(temp.sig)
  global.df100 <- BeatAML$global[,colSums(is.na(BeatAML$global)) == 0]
  if (nrow(temp.sig) > 0) {
    # perform weighted voting on input omics
    temp.sig <- temp.sig[,c("Gene", colnames(temp.sig)[1])]
    temp.sig2 <- temp.sig[temp.sig$Gene %in% colnames(global.df100)[2:ncol(global.df100)],]
    if (nrow(temp.sig2) > 0) {
      global.df100 <- global.df100[,c("Barcode.ID",temp.sig2$Gene)]
      if (nrow(global.df100) > 0) {
        sorted.wv <- panSEA::WV(global.df100, temp.sig2)
        temp.wv <- sorted.wv$scores }
      } else {
        temp.wv <- data.frame()
      }
    } else {
      temp.wv <- data.frame()
    }
    
  AUC.df <- merge(temp.wv, BeatAML$drug, by="Barcode.ID")
  Drug <- colnames(AUC.df)[3:ncol(AUC.df)]
  Barcode.ID <- AUC.df$Barcode.ID
  pt.corr <- data.frame(Barcode.ID, SSE = NA, Pearson.est = NA, Pearson.p = NA,
                        Spearman.est = NA, Spearman.p = NA, N = NA)
  all.corr <- data.frame()
  for (i in 1:nrow(AUC.df)) {
    # leave out one patient at a time
    temp.AUC <- AUC.df[-i,]
    corr <- data.frame(Drug, Slope = NA, Intercept = NA, R.squared = NA, N = NA)
    # for each drug:
    for (j in 3:ncol(temp.AUC)) {
      # determine line of best fit 
      x <- as.numeric(temp.AUC$WV)
      y <- as.numeric(temp.AUC[,j])
      Regression <- stats::lm(y ~ x)
      corr$Slope[j-2] <- Regression$coeff[[2]]
      corr$Intercept[j-2] <- Regression$coeff[[1]]
      corr$R.squared[j-2] <- summary(Regression)$r.squared
      corr$N[j-2] <- length(x)
    }
    corr$Barcode.ID <- AUC.df$Barcode.ID[i]
    
    # predict drug sensitivity using line of best fit
    corr$WV <- temp.wv[temp.wv$Barcode.ID == AUC.df$Barcode.ID[i],]$WV
    corr$AUC_predicted <- corr$Intercept + corr$Slope * corr$WV
    
    # calculate accuracy
    corr$AUC_measured <- as.numeric(as.vector(AUC.df[i,Drug]))
    corr$delta_AUC_squared <- ifelse(is.na(corr$AUC_predicted) | is.na(corr$AUC_measured),
                                     NA, (corr$AUC_predicted - corr$AUC_measured)^2)
    all.corr <- rbind(all.corr, corr)
    pt.corr[pt.corr$Barcode.ID == AUC.df$Barcode.ID[i],]$SSE <- sum(corr$delta_AUC_squared)
    
    corr.df <- na.omit(corr)
    pt.corr[pt.corr$Barcode.ID == AUC.df$Barcode.ID[i],]$N <- nrow(corr.df)
    if (nrow(corr.df) > 2) {
      Pearson <- stats::cor.test(corr.df$AUC_measured, corr.df$AUC_predicted, method = "pearson")
      pt.corr[pt.corr$Barcode.ID == AUC.df$Barcode.ID[i],]$Pearson.est <- Pearson$estimate
      pt.corr[pt.corr$Barcode.ID == AUC.df$Barcode.ID[i],]$Pearson.p <- Pearson$p.value
      
      Spearman <- stats::cor.test(corr.df$AUC_measured, corr.df$AUC_predicted, method = "spearman")
      pt.corr[pt.corr$Barcode.ID == AUC.df$Barcode.ID[i],]$Spearman.est <- Spearman$estimate
      pt.corr[pt.corr$Barcode.ID == AUC.df$Barcode.ID[i],]$Spearman.p <- Spearman$p.value 
    }
  }
  
  # for each drug, evaluate accuracy
  drug.corr <- data.frame(Drug, SSE = NA, Pearson.est = NA, Pearson.p = NA,
                        Spearman.est = NA, Spearman.p = NA, N = NA)
  for (j in Drug) {
    drug.df <- na.omit(all.corr[all.corr$Drug == j,])
    drug.corr[drug.corr$Drug == j,]$N <- nrow(drug.df)
    drug.corr$SSE <- sum(drug.df$delta_AUC_squared)
    
    if (nrow(drug.df) > 2) {
      Pearson <- stats::cor.test(drug.df$AUC_measured, drug.df$AUC_predicted, method = "pearson")
      drug.corr[drug.corr$Drug == j,]$Pearson.est <- Pearson$estimate
      drug.corr[drug.corr$Drug == j,]$Pearson.p <- Pearson$p.value
      
      Spearman <- stats::cor.test(drug.df$AUC_measured, drug.df$AUC_predicted, method = "spearman")
      drug.corr[drug.corr$Drug == j,]$Spearman.est <- Spearman$estimate
      drug.corr[drug.corr$Drug == j,]$Spearman.p <- Spearman$p.value 
    }
  }
  
  # run DMEA to see which MOAs have the best accuracy based on corr estimate
  dmeaInput <- na.omit(drug.corr[,c("Drug", "Pearson.est", "Spearman.est")])
  if (nrow(dmeaInput) > 6) {
    dmeaPearson <- panSEA::drugSEA_ties(drug.corr, gmt)
    dmeaSpearman <- panSEA::drugSEA_ties(drug.corr, gmt, rank.metric="Spearman.est") 
  } else {
    dmeaPearson <- list()
    dmeaSpearman <- list()
  }
  
  return(list(wv = temp.wv, drug = all.corr, 
              pt.corr = pt.corr, drug.corr = drug.corr, 
              DMEA = dmeaPearson, DMEA.Spearman = dmeaSpearman))
}

evalMonoSig <- function(sig.matrix, BeatAML, types=rep("rna",ncol(sig.matrix)), gmt) {
  # evaluate each signature
  wv.df <- data.frame()
  drug.df <- data.frame()
  pt.corr.df <- data.frame()
  drug.corr.df <- data.frame()
  sig.matrix$Gene <- rownames(sig.matrix)
  DMEA.results <- list()
  DMEA.results.Spearman <- list()
  for (i in 1:(ncol(sig.matrix)-1)) {
    cat("evaluating",names(sig.matrix)[i],"as",types[i],"\n")
    temp.result <- evalOneMonoSig(sig.matrix[,c(i,ncol(sig.matrix))],
                                  BeatAML, types[i], gmt)
    DMEA.results[[names(sig.matrix)[i]]] <- temp.result$DMEA
    DMEA.results.Spearman[[names(sig.matrix)[i]]] <- temp.result$DMEA.Spearman
    temp.wv.df <- temp.result$wv
    temp.drug.df <- temp.result$drug
    temp.pt.corr.df <- temp.result$pt.corr
    temp.drug.corr.df <- temp.result$drug.corr
    temp.wv.df$Signature <- names(sig.matrix)[i]
    temp.drug.df$Signature <- names(sig.matrix)[i]
    temp.pt.corr.df$Signature <- names(sig.matrix)[i]
    temp.drug.corr.df$Signature <- names(sig.matrix)[i]
    
    wv.df <- rbind(wv.df, temp.wv.df)
    drug.df <- rbind(drug.df, temp.drug.df)
    pt.corr.df <- rbind(pt.corr.df, temp.pt.corr.df)
    drug.corr.df <- rbind(drug.corr.df, temp.drug.corr.df)
  }
  
  return(list(wv = wv.df, drug = drug.df, pt.corr = pt.corr.df, 
              drug.corr = drug.corr.df, DMEA = DMEA.results,
              DMEA.Spearman = DMEA.results.Spearman))
}

compareSigs <- function(sigs, value.var = "Log2FC", BeatAML, 
                        types=rep("rna",length(sigs)), gmt, 
                        fillVals = RColorBrewer::brewer.pal(length(sigs), "Set2")) {
  # combine signatures into matrix
  filtered.sigs.df <- data.table::rbindlist(sigs, use.names = TRUE, idcol = "Signature")
  sig.matrix <- reshape2::dcast(filtered.sigs.df, Gene ~ Signature, mean,
                                value.var = value.var)
  rownames(sig.matrix) <- sig.matrix$Gene
  sig.matrix$Gene <- NULL
  sig.matrix <- sig.matrix[,names(sigs)]
  
  # test signature matrix
  sigResults <- evalMonoSig(sig.matrix, BeatAML,types,gmt)
  wv.df <- sigResults$wv
  drug.df <- sigResults$drug
  pt.corr.df <- sigResults$pt.corr
  drug.corr.df <- sigResults$drug.corr
  
  # compare accuracy for Aza, Ven, Aza+Ven across signatures
  rank.metrics <- c("Pearson.est", "Spearman.est", "SSE")
  for (i in rank.metrics) {
    descr <- stringr::str_split_1(i, "[.]")[1]
    if ("Drug" %in% colnames(drug.corr.df)) {
      doi <- c("Azacytidine", "Venetoclax", "Azacytidine - Venetoclax")
      drug.corr.df$`Drug Treatment` <- NA
      drug.corr.df[drug.corr.df$Drug == "Azacytidine",]$`Drug Treatment` <- "Aza"
      drug.corr.df[drug.corr.df$Drug == "Azacytidine - Venetoclax",]$`Drug Treatment` <- "Aza + Ven"
      drug.corr.df[drug.corr.df$Drug == "Venetoclax",]$`Drug Treatment` <- "Ven"
      doi.names <- c("Aza", "Aza + Ven", "Ven")
      
      # aza, aza + ven, ven correlations
      plot.df <- drug.corr.df[drug.corr.df$Drug %in% doi,]
      plot.df$rank <- plot.df[,i]
      sigOrder <- na.omit(unique(plot.df[order(plot.df$rank, decreasing=TRUE),]$Signature))
      ylab <- paste0(descr," Correlation Estimate")
      if (grepl("SSE", i)) {
        ylab <- "Sum Squared Error"
      }
      ggplot(na.omit(plot.df), aes(x=Signature, y=rank, fill = Signature)) + 
        geom_col(alpha=0.5) + theme_minimal(base_size = 12) + ylab(ylab) + 
        facet_wrap(~ `Drug Treatment`) +
        ggplot2::scale_x_discrete(limits = sigOrder) +
        theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
        scale_fill_manual(values=fillVals, 
                          breaks=c("Sorted","Lasry","Triana","van Galen"))+
        ggtitle("Monocytic signatures predict drug sensitivity")
      ggsave(paste0("AzaVen_DIA_WV_signatureFill_",descr,".pdf"), width = 5, height = 5)
      
      for (j in doi.names) {
        plot.df <- drug.corr.df[drug.corr.df$`Drug Treatment` == j,]
        plot.df$rank <- plot.df[,i]
        sigOrder <- na.omit(unique(plot.df[order(plot.df$rank, decreasing=TRUE),]$Signature))
        plot.df$alpha <- 0.5
        circBar(na.omit(plot.df), x="Signature", y = "rank", fill = "Signature", 
                alpha = "alpha", ymin = 0, ymax = 1, alpha_range=0.5, 
                ytick_yScale = 2/3, ytick_yShift = 0, fillVals=fillVals,
                title=paste("Monocytic signatures predict", j, "sensitivity"),
                fname=paste0(j,"_", descr, "_DIA_WV_signatureFill_circBarPlot.pdf"))
        ggplot(plot.df, aes(x=Signature, y=rank, fill = Signature)) + 
          geom_col(alpha=0.5) + theme_classic(base_size = 12) + 
          ylab(ylab) + 
          ggplot2::scale_x_discrete(limits = sigOrder) +
          scale_fill_manual(values=fillVals, 
                            breaks=c("Sorted","Lasry","Triana","van Galen"))+
          ggtitle(paste("Monocytic signatures predict", j, "sensitivity"))
        ggsave(paste0(j,"_", descr, "_DIA_WV_signatureFill_barPlot.pdf"), width = 5, height = 5)
      }
    }
  }
  
  # compare accuracy across all drugs for each signature
  rank.metrics <- c("Pearson.est", "Spearman.est", "SSE")
  for (i in rank.metrics) {
    descr <- stringr::str_split_1(i, "[.]")[1]
    if ("Drug" %in% colnames(drug.corr.df)) {
      plot.df <- drug.corr.df
      plot.df$rank <- plot.df[,i]
      sigOrder <- na.omit(unique(plot.df[order(plot.df$rank, decreasing=TRUE),]$Signature))
      ylab <- paste0(descr," Correlation Estimate")
      if (grepl("SSE", i)) {
        ylab <- "Sum Squared Error"
      }
      ggplot(na.omit(plot.df), aes(x=Signature, y=rank, fill = Signature)) + 
        geom_col(alpha=0.5) + theme_minimal(base_size = 12) + ylab(ylab) + 
        facet_wrap(~ Drug) +
        ggplot2::scale_x_discrete(limits = sigOrder) +
        theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
        scale_fill_manual(values=fillVals, 
                          breaks=c("Sorted","Lasry","Triana","van Galen"))+
        ggtitle("Monocytic signatures predict drug sensitivity")
      ggsave(paste0("Drug_DIA_WV_signatureFill_",descr,".pdf"), width = 5, height = 5)
    }
  }
  
  return(list(wv = wv.df, drug = drug.df, pt.corr = pt.corr.df,
              drug.corr = drug.corr.df, DMEA = sigResults$DMEA,
              DMEA.Spearman = sigResults$DMEA.Spearman))
}

load_not_norm_BeatAML_for_DMEA3 <- function(BeatAML.path = "BeatAML_DMEA_inputs_not_normalized",
                                            exclude.samples = c()) {
  message("Loading Beat AML data for DMEA")
  BeatAML_synapse_id <- list("drug_response.csv" = "syn51674470", 
                             "Ex10_metadata.txt" = "syn25807733",
                             "ptrc_ex10_crosstab_global_gene_original.txt" = "syn25714254",
                             "ptrc_ex10_crosstab_phospho_siteID_original.txt" = "syn25714936")
  
  ### download files if any not already downloaded
  if (!file.exists(BeatAML.path)) {
    lapply(BeatAML_synapse_id, synapser::synGet, downloadLocation = BeatAML.path)
  } else if (!any(FALSE %in% lapply(names(BeatAML_synapse_id), file.exists))) {
    lapply(BeatAML_synapse_id, synapser::synGet, downloadLocation = BeatAML.path)
  }
  
  ### load files
  drug.BeatAML <- read.csv(file.path(BeatAML.path, names(BeatAML_synapse_id)[1]))
  meta.BeatAML <- read.table(file.path(BeatAML.path, names(BeatAML_synapse_id)[2]), 
                             sep = "\t", header = TRUE)
  global.BeatAML <- read.table(file.path(BeatAML.path, names(BeatAML_synapse_id)[3]),
                               sep = "\t", header = TRUE)
  phospho.BeatAML <- read.table(file.path(BeatAML.path, names(BeatAML_synapse_id)[4]),
                                sep = "\t", header = TRUE)
  rna.BeatAML <- synapser::synTableQuery("select * from syn26545877")$asDataFrame()
  
  ### format BeatAML data for DMEA
  sample.names <- "Barcode.ID"
  
  ## format drug sensitivity data frame
  # format drug.BeatAML wide (samples in first column, drug names for rest of columns)
  drug.BeatAML <- reshape2::dcast(drug.BeatAML, sample_id ~ inhibitor, 
                                  value.var = "auc", fill = NA)
  
  # change sample column name to match expression data
  names(drug.BeatAML)[1] <- sample.names
  
  ## format global proteomics data frame
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
  
  # subtract sample medians
  sample.names <- colnames(dplyr::select_if(global.BeatAML, is.numeric))
  #global.BeatAML[,sample.names] <- log(global.BeatAML[,sample.names], 2)
  global_sample_coef <- apply(global.BeatAML[,sample.names], 2, median, na.rm = T)
  global.BeatAML[,sample.names] <- sweep(global.BeatAML[,sample.names], 2, global_sample_coef, FUN = '-')
  
  # transpose global.BeatAML so that first column is Barcode.ID and 
  # rest of columns are gene symbols
  global.BeatAML <- as.data.frame(t(global.BeatAML))
  
  # make first column Barcode.ID
  global.BeatAML[,"Barcode.ID"] <- rownames(global.BeatAML)
  global.BeatAML <- 
    global.BeatAML[ , c("Barcode.ID", 
                        names(global.BeatAML[ , 1:(ncol(global.BeatAML)-1)]))]
  
  ## format phospho-proteomics data frame
  # change global.BeatAML column names from SampleID.abbrev to Barcode.ID to match drug.BeatAML
  phospho.ids <- names(phospho.BeatAML)
  
  # remove X and any 0's from start of each column name and then
  # replace SampleID.abbrev with Barcode.ID to match drug.BeatAML
  for(i in seq_len(length(phospho.ids))){
    phospho.ids[i] <- substr(phospho.ids[i], 2, nchar(phospho.ids[i]))
    
    if(substring(phospho.ids[i], 1, 1) == 0){
      phospho.ids[i] <- substr(phospho.ids[i], 2, nchar(phospho.ids[i]))
    }
    
    if(phospho.ids[i] %in% meta.BeatAML$SampleID.abbrev){
      phospho.ids[i] <- meta.BeatAML[
        meta.BeatAML$SampleID.abbrev == phospho.ids[i], ]$Barcode.ID
    }
  }
  
  # replace phospho.BeatAML column names
  names(phospho.BeatAML) <- phospho.ids
  
  # subtract sample medians
  sample.names <- colnames(dplyr::select_if(phospho.BeatAML, is.numeric))
  #phospho.BeatAML[,sample.names] <- log(phospho.BeatAML[,sample.names], 2)
  phospho_sample_coef <- apply(phospho.BeatAML[,sample.names], 2, median, na.rm = T)
  phospho.BeatAML[,sample.names] <- sweep(phospho.BeatAML[,sample.names], 2, phospho_sample_coef, FUN = '-')
  
  # transpose phospho.BeatAML so that first column is Barcode.ID and rest of columns are gene symbols
  phospho.BeatAML <- as.data.frame(t(phospho.BeatAML))
  
  # make first column Barcode.ID
  phospho.BeatAML[, "Barcode.ID"] <- rownames(phospho.BeatAML)
  phospho.BeatAML <- phospho.BeatAML[ , c("Barcode.ID", names(phospho.BeatAML[ , 1:(ncol(phospho.BeatAML)-1)]))]
  
  ## format rnaSeq data frame
  rna.BeatAML$Barcode.ID <- rna.BeatAML$labId
  rna.BeatAML <- reshape2::dcast(rna.BeatAML, Barcode.ID ~ display_label, mean,
                                 value.var = "RNA counts")
  rownames(rna.BeatAML) <- rna.BeatAML$Barcode.ID
  
  return(list(meta = meta.BeatAML[!(meta.BeatAML$Barcode.ID %in% exclude.samples),], 
              drug = drug.BeatAML[!(drug.BeatAML$Barcode.ID %in% exclude.samples),], 
              rna = rna.BeatAML[!(rna.BeatAML$Barcode.ID %in% exclude.samples),],
              global = global.BeatAML[!(global.BeatAML$Barcode.ID %in% exclude.samples),],
              phospho = phospho.BeatAML[!(phospho.BeatAML$Barcode.ID %in% exclude.samples),]))
}

#### predict using full signatures ####
setwd("~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/")
setwd("data")

# drug.info <- read.csv("~/OneDrive - PNNL/Documents/PTRC2/BeatAML_single_drug_moa_2025-01-20.csv",
#                       stringsAsFactors=FALSE, fileEncoding="latin1")
# gmt.drug <- DMEA::as_gmt(drug.info, sep=", ")
# saveRDS(gmt.drug, "gmt_BeatAML_drug_MOA_2025-01-20.rds")
gmt.drug <- readRDS("gmt_BeatAML_drug_MOA_2025-01-20.rds")

# load sorted proteomics signature
# sig.paths <- list("Sorted" = "analysis/combined24-27/DIA_2batches_noOutliers_noMSC/Sort Type_Bead/CD14_Pos_vs_Neg/global/Differential_expression/Differential_expression_results.csv",
#                   "van Galen" = "data/externalSignatures/formatted/notFilteredForMalignant/Differential_expression_van_Galen_AML_D0_Mono-like_vs_Prog-like_noNA.csv",
#                   "Triana" = "data/externalSignatures/formatted/Triana_RNA_AML_100PercentCells_Classical-Monocytes_vs_HSCs-and-MPPs_differentialExpression.csv",
#                   "Lasry" = "data/externalSignatures/formatted/notFilteredForMalignant/Differential_expression_Lasry_AML_CD14PosMonocyte_vs_HSC_protein-coding.csv")

sig.paths <- list("Sorted" = "analysis/combined24-27/using_cellType-sortType-patient_factors/DIA_2batches_noOutliers_noMSC_Bead/no_filter/CD14_Pos_vs_Neg/global/Differential_expression/Differential_expression_results.csv",
                  "van Galen" = "data/externalSignatures/formatted/notFilteredForMalignant/Differential_expression_van_Galen_AML_D0_Mono-like_vs_Prog-like.csv",
                  "Triana" = "data/externalSignatures/formatted/Triana_RNA_AML_100PercentCells_Classical-Monocytes_vs_HSCs-and-MPPs_differentialExpression.csv",
                  "Lasry" = "data/externalSignatures/formatted/notFilteredForMalignant/Differential_expression_Lasry_AML_CD14PosMonocyte_vs_MPP.csv")

# import signatures and filter
sigs <- list()
setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/")
for (i in names(sig.paths)) {
  sigs[[i]] <- read.csv(sig.paths[[i]])
  sigs[[i]] <- na.omit(sigs[[i]][sigs[[i]]$adj.P.Val <= 0.05,c("Gene","Log2FC")])
}
sigs <- readRDS("mono_vs_prog_sigs_ProteinCoding.rds")
setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/")
dir.create("Monocyte_vs_progenitor_signatures_beadOnly_LOO_2025-09-11")
setwd("Monocyte_vs_progenitor_signatures_beadOnly_LOO_2025-09-11")

dia.wo.out <- readRDS("~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/analysis/DIA_2batches_noOutliers.rds")
# sorted.patients <- c("18-00105", "21-00839", "22-00571", "22-00117", "16-01184",
#                      "19-00074", "18-00103", "21-00432", "17-01060", "22-00251")
sorted.patients <- unique(dia.wo.out$meta$patient)
BeatAML <- load_not_norm_BeatAML_for_DMEA3(exclude.samples=sorted.patients)

# evaluate signature
evalResults <- compareSigs(sigs, BeatAML = BeatAML, types=c("global", "rna", "rna", "rna"), gmt = gmt.drug) 
write.csv(evalResults$wv, "wv.csv", row.names = FALSE)
write.csv(evalResults$drug, "predictions.csv", row.names = FALSE)
write.csv(evalResults$drug.corr, "drugAccuracy.csv", row.names = FALSE)
write.csv(evalResults$pt.corr, "patientAccuracy.csv", row.names = FALSE)
saveRDS(evalResults$DMEA, "DMEA.rds")
saveRDS(evalResults$DMEA.Spearman, "DMEA_Spearman.rds")
all.DMEA.files <- list()
for (i in names(sigs)) {
  DMEA.files <- list("DMEA_results.csv" =
                       evalResults$DMEA[[i]]$result,
                     "DMEA_results_Spearman.csv" =
                       evalResults$DMEA.Spearman[[i]]$result,
                     "DMEA_volcano_plot.pdf" =
                       evalResults$DMEA[[i]]$volcano.plot,
                     "DMEA_volcano_plot_Spearman.pdf" =
                       evalResults$DMEA.Spearman[[i]]$volcano.plot,
                     "DMEA_bar_plot.pdf" =
                       evalResults$DMEA[[i]]$bar.plot,
                     "DMEA_bar_plot_Spearman.pdf" =
                       evalResults$DMEA.Spearman[[i]]$bar.plot,
                     "DMEA_dot_plot.pdf" =
                       evalResults$DMEA[[i]]$dot.plot,
                     "DMEA_dot_plot_Spearman.pdf" =
                       evalResults$DMEA.Spearman[[i]]$dot.plot) 
  all.DMEA.files[[i]] <- DMEA.files
}
save_to_synapse_v2(all.DMEA.files)

# redo plots
setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/")
dir.create("Monocyte_vs_progenitor_signatures_beadOnly_LOO_2025-09-11")
setwd("Monocyte_vs_progenitor_signatures_beadOnly_LOO_2025-09-11")
drug.corr.df <- read.csv("drugAccuracy.csv")
pt.corr.df <- read.csv("patientAccuracy.csv")
p.df <- pt.corr.df
sigOrder <- p.df[order(p.df$Pearson.est),]$Signature
p.df$Signature <- factor(p.df$Signature, levels=unique(sigOrder))
ggplot2::ggplot(p.df, aes(x=Signature, y=Pearson.est)) + geom_violin(alpha=0) +
  geom_point(#aes(color=Drug)
    ) + 
  geom_boxplot(width=0.2, alpha = 0) + 
  labs(y="Pearson Correlation Estimate") + theme_classic(base_size = 12) +
  ggtitle(paste("Monocyte signatures predict drug sensitivity"))
ggsave(paste0("patientAccuracy","_bySignature.pdf"), width = 5, height = 5)

p.df <- drug.corr.df
sigOrder <- p.df[order(p.df$Pearson.est),]$Signature
p.df$Signature <- factor(p.df$Signature, levels=unique(sigOrder))
ggplot2::ggplot(p.df, aes(x=Signature, y=Pearson.est)) + geom_violin(alpha=0) +
  geom_point(#aes(color=Drug)
  ) + 
  geom_boxplot(width=0.2, alpha = 0) + 
  labs(y="Pearson Correlation Estimate") + theme_classic(base_size = 12) +
  ggtitle(paste("Monocyte signatures predict drug sensitivity"))
ggsave(paste0("drugAccuracy","_bySignature.pdf"), width = 5, height = 5)


#frac.corr.df <- read.csv("cellFractionCorrelations.csv")

drug.info <- read.csv("~/OneDrive - PNNL/Documents/PTRC2/BeatAML_single_drug_moa.csv",
                      stringsAsFactors = FALSE, fileEncoding = "latin1")
drug.info <- drug.info[,c("Drug","moa")]
drug.info[drug.info$Drug == "Ralimetinib (LY2228820)",]$moa <- "p38 MAPK inhibitor"
drug.info[drug.info$Drug == "Nilotinib",]$moa <- "Abl kinase inhibitor"
drug.info[drug.info$Drug == "AT-101",]$moa <- "BCL inhibitor"
drug.info[is.na(drug.info$moa),]$moa <- "Other"
library(patchwork); library(ggplot2)
pearson.plots <- NULL
spearman.plots <- NULL
#MOAsInTop50 <- names(gmt.drug$genesets)
#moaColors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(length(MOAsInTop50))
pearson.venn <- list()
spearman.venn <- list()
for (i in unique(drug.corr.df$Signature)) {
  p.df <- drug.corr.df[drug.corr.df$Signature == i,]
  p.df$Pearson.q <- qvalue::qvalue(p = p.df$Pearson.p, pi0 = 1)$qvalues
  p.df$Spearman.q <- qvalue::qvalue(p = p.df$Spearman.p, pi0 = 1)$qvalues
  plot.df <- merge(p.df, drug.info, by="Drug", all.x = TRUE)
  plot.df$Mechanism <- "Other"
  #plot.df[plot.df$moa %in% MOAsInTop50,]$Mechanism <- plot.df[plot.df$moa %in% MOAsInTop50,]$moa
  plot.df[grepl("Venetoclax",plot.df$Drug),]$Mechanism <- "BCL inhibitor"
  plot.df$Drug <- sub(" [(].*", "", plot.df$Drug) # shorten drug names for plot
  plot.df[plot.df$Drug == "NF-kB Activation Inhibitor",]$Drug <- "NFkB Inhibitor"
  
  rank.metrics <- c("Pearson.est", "Spearman.est")
  for (j in rank.metrics) {
    descr <- stringr::str_split_1(j, "[.]")[1]
    if ("Drug" %in% colnames(plot.df)) {
      if (j == "Pearson.est") {
        plot.df <- plot.df[plot.df$Pearson.est > 0 & plot.df$Pearson.q <= 0.05,]
        ylab <- paste0(descr," r")
        pearson.venn[[i]] <- unique(plot.df$Drug)
      } else {
        plot.df <- plot.df[plot.df$Spearman.est > 0 & plot.df$Spearman.q <= 0.05,]
        ylab <- paste0(descr," rho")
        spearman.venn[[i]] <- unique(plot.df$Drug)
      }
      plot.df$rank <- plot.df[,j]
      sigOrder <- na.omit(unique(plot.df[order(plot.df$rank, decreasing=TRUE),]$Drug))
      plot.annot <- paste0(i, "\n(", nrow(plot.df), " / ", nrow(p.df), " Drugs Positively Correlated)")
      corr.plot <- ggplot(plot.df, aes(x=Drug, y=rank, fill = Mechanism)) + 
        geom_col() + theme_minimal(base_size = 12) + ylab(ylab) + 
        ggplot2::scale_x_discrete(limits = sigOrder) +
        theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1),
              axis.title.x=element_blank()) +
        #scale_fill_manual(breaks=MOAsInTop50, values = moaColors) +
        ggtitle(plot.annot) + 
        theme(plot.title = element_text(hjust = 0.5, face="bold", size=16), legend.position="bottom")
      ggsave(paste0("Drug_DIA_WV_moaFill_",descr,"_", i, ".pdf"), corr.plot, width = 10, height = 5)
      if (is.null(pearson.plots) & j == "Pearson.est") {
        pearson.plots <- (corr.plot + theme(legend.position = "none"))
      } else if (j == "Pearson.est") {
        pearson.plots <- pearson.plots / (corr.plot + theme(legend.position = "none"))
      } else if (is.null(spearman.plots) & j == "Spearman.est") {
        spearman.plots <- (corr.plot + theme(legend.position = "none"))
      } else if (j == "Spearman.est") {
        spearman.plots <- spearman.plots / (corr.plot + theme(legend.position = "none"))
      }
    }
  }
}
#source("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/MPNST/Chr8/MPNST_Chr8_manuscript/Figure_3_Kinase/guides_build_mod.R")
#pearson.plots <- (pearson.plots / plot_spacer()) + plot_layout(guides='collect')
#pearson.plots <- pearson.plots + theme(legend.position = "none")
ggplot2::ggsave("Drug_DIA_WV_moaFill_Pearson_allSigs.pdf", pearson.plots, width=12, height=12)
ggplot2::ggsave("Drug_DIA_WV_moaFill_Spearman_allSigs.pdf", spearman.plots, width=12, height=12)
ggvenn::ggvenn(pearson.venn, show_percentage=FALSE, set_name_size=5, text_size=5)
ggsave("Drug_DIA_WV_Pearson_sigOverlap.pdf", width=5, height=5)
ggvenn::ggvenn(spearman.venn, show_percentage=FALSE, set_name_size=5, text_size=5)
ggsave("Drug_DIA_WV_Spearman_sigOverlap.pdf", width=5, height=5)

#### narrow down sorted signature ####
# "~/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Exp24_patient_cells/2025-01-20_v2_MonoSigComparison_cellFrac.R" suggests 
# NCF2, FCGRT, KCTD12, CD93 (unweighted) are best at predicting Ven AUC (all 4 quantified in all 122) (r=0.79, q=2E-6)
# NCF2 alone is also sufficient (r=0.73, q=6E-6)
setwd("~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/")
setwd("data")

# drug.info <- read.csv("~/OneDrive - PNNL/Documents/PTRC2/BeatAML_single_drug_moa_2025-01-20.csv",
#                       stringsAsFactors=FALSE, fileEncoding="latin1")
# gmt.drug <- DMEA::as_gmt(drug.info, sep=", ")
# saveRDS(gmt.drug, "gmt_BeatAML_drug_MOA_2025-01-20.rds")
gmt.drug <- readRDS("gmt_BeatAML_drug_MOA_2025-01-20.rds")

# import signatures and filter
#topPred <- read.csv("Monocyte_vs_progenitor_signatures_beadOnly_2025-06-30/topVenSensPredictions_2025-06-30.csv")
setwd("~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/")
topPred <- read.csv("Monocyte_vs_progenitor_signatures_beadOnly_2025-07-14/topVenSensPredictions_2025-07-15.csv")
sigs <- list()
setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/")
types <- c()
for (i in 1:nrow(topPred)) {
  temp.name <- ifelse(topPred$Signature[i]=="Sorted",paste0(topPred$Signature[i], "_", topPred$N_genes[i],"_Protein"),paste0(topPred$Signature[i], "_", topPred$N_genes[i],"_Gene"))
  if (topPred$N_genes[i]>1) {
    temp.name <- paste0(temp.name,'s') # make plural if appropriate
  }
  Gene <- strsplit(topPred$Genes[i], ", ")[[1]]
  sigs[[temp.name]] <- data.frame(Gene, Log2FC=1)
  if (topPred$Signature[i] == "Sorted") {
    types <- c(types, "global")
  } else {
    types <- c(types, "rna")
  }
}

# try just mixed sigs
topMixPred <- read.csv("Monocyte_vs_progenitor_signatures_beadOnly_2025-07-14/topMixVenSensPredictions_2025-07-15.csv")
sigs2 <- readRDS("Monocyte_vs_progenitor_signatures_beadOnly_2025-07-14/topMixVenSensSigs_2025-07-15.rds")
types <- c("global","rna")

setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/")
#dir.create("Monocyte_vs_progenitor_signatures_beadOnly_2025-06-30_LOO")
#setwd("Monocyte_vs_progenitor_signatures_beadOnly_2025-06-30_LOO")
#dir.create("Monocyte_vs_progenitor_signatures_beadOnly_2025-07-09_LOO")
#setwd("Monocyte_vs_progenitor_signatures_beadOnly_2025-07-09_LOO")
dir.create("Monocyte_vs_progenitor_signatures_beadOnly_2025-09-11_LOO")
setwd("Monocyte_vs_progenitor_signatures_beadOnly_2025-09-11_LOO")

dia.wo.out <- readRDS("~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/analysis/DIA_2batches_noOutliers.rds")
# sorted.patients <- c("18-00105", "21-00839", "22-00571", "22-00117", "16-01184",
#                      "19-00074", "18-00103", "21-00432", "17-01060", "22-00251")
sorted.patients <- unique(dia.wo.out$meta$patient)
BeatAML <- load_not_norm_BeatAML_for_DMEA3(exclude.samples=sorted.patients)

# evaluate signature
#evalResults <- compareSigs(sigs, BeatAML = BeatAML, types=c("global", "rna", "rna", "rna"), gmt = gmt.drug) 
#evalResults <- compareSigs(sigs, BeatAML = BeatAML, types=types, gmt = gmt.drug) 
evalResults <- compareSigs(sigs2, value.var="Pearson.est", BeatAML = BeatAML, types=types, gmt = gmt.drug) 
write.csv(evalResults$wv, "wv.csv", row.names = FALSE)
write.csv(evalResults$drug, "predictions.csv", row.names = FALSE)
write.csv(evalResults$drug.corr, "drugAccuracy.csv", row.names = FALSE)
write.csv(evalResults$pt.corr, "patientAccuracy.csv", row.names = FALSE)
saveRDS(evalResults$DMEA, "DMEA.rds")
saveRDS(evalResults$DMEA.Spearman, "DMEA_Spearman.rds")
all.DMEA.files <- list()
for (i in names(sigs)) {
  DMEA.files <- list("DMEA_results.csv" =
                       evalResults$DMEA[[i]]$result,
                     "DMEA_results_Spearman.csv" =
                       evalResults$DMEA.Spearman[[i]]$result,
                     "DMEA_volcano_plot.pdf" =
                       evalResults$DMEA[[i]]$volcano.plot,
                     "DMEA_volcano_plot_Spearman.pdf" =
                       evalResults$DMEA.Spearman[[i]]$volcano.plot,
                     "DMEA_bar_plot.pdf" =
                       evalResults$DMEA[[i]]$bar.plot,
                     "DMEA_bar_plot_Spearman.pdf" =
                       evalResults$DMEA.Spearman[[i]]$bar.plot,
                     "DMEA_dot_plot.pdf" =
                       evalResults$DMEA[[i]]$dot.plot,
                     "DMEA_dot_plot_Spearman.pdf" =
                       evalResults$DMEA.Spearman[[i]]$dot.plot) 
  all.DMEA.files[[i]] <- DMEA.files
}
save_to_synapse_v2(all.DMEA.files)

evalResults <- list("drug" = read.csv("predictions.csv"),
                    "drug.corr" = read.csv("drugAccuracy.csv"),
                    "pt.corr" = read.csv("patientAccuracy.csv"))

hist(evalResults$pt.corr$Pearson.est)
med.pt.r <- median(evalResults$pt.corr$Pearson.est)
med.pt.r # 0.786825
med.pt.r.s <- median(evalResults$pt.corr[evalResults$pt.corr$Signature=="Sorted: 25 proteins",]$Pearson.est)
med.pt.r.s # 0.7885036
med.pt.r.l <- median(evalResults$pt.corr[evalResults$pt.corr$Signature=="Lasry: 46 genes",]$Pearson.est)
med.pt.r.l # 0.7841898

sd.pt.r <- sd(evalResults$pt.corr$Pearson.est)
sd.pt.r # 0.1684279
frac.sd <- sd.pt.r/med.pt.r
frac.sd # 0.2140601
pt.r <- reshape2::dcast(evalResults$pt.corr, Barcode.ID ~ Signature, value.var="Pearson.est")

multi.cor <- cor.test(pt.r$`Lasry: 46 genes`, pt.r$`Sorted: 25 proteins`)
Pearson.est <- multi.cor$estimate
Pearson.p <- multi.cor$p.value
stats_pearson <- substitute(
  r == est * "," ~ ~"p" ~ "=" ~ p,
  list(
    est = format(as.numeric(Pearson.est), digits = 3),
    p = format(Pearson.p, digits = 3)
  )
)
ggplot(pt.r, aes(x=`Lasry: 46 genes`, y=`Sorted: 25 proteins`)) + geom_point() + theme_minimal() + 
  scale_x_continuous(limits=c(-1,1)) + scale_y_continuous(limits=c(-1,1)) + 
  labs(x="Lasry RNA-seq: 46 Genes", y = "Sorted Proteomics: 25 proteins", title="Drug Sensitivity Prediction Accuracy\nfor Each Patient (Pearson Correlation)") +
  geom_smooth(se=FALSE, linetype="dashed") + ggrepel::geom_label_repel(aes(label=Barcode.ID)) + theme(plot.title=element_text(hjust=0.5)) +
  ggplot2::geom_text(
    x = Inf, y = -Inf, vjust = "inward", hjust = "inward",
    colour = "blue", parse = TRUE,
    label = as.character(as.expression(stats_pearson)), size = 5
  ) 
ggsave("Lasry_vs_Sorted_multipleGenes_ptCorr.pdf", width=4, height=4)

multi.test <- t.test(evalResults$pt.corr[evalResults$pt.corr$Signature == "Lasry: 46 genes",]$Pearson.est, # mean 0.7444206
                     evalResults$pt.corr[evalResults$pt.corr$Signature == "Sorted: 25 proteins",]$Pearson.est, # mean 0.7473072 
                     #alternative = "less"
)
multi.test$estimate
multi.test$p.value
# greater p = 0.5685869; two-sided p = 0.8628261; less p = 0.4314131

ggplot(evalResults$pt.corr, aes(x=Signature, y=Pearson.est)) + geom_violin() + 
  geom_boxplot() + geom_point() + theme_classic() + 
  labs(y="Pearson r", title="Drug Sensitivity Prediction\nAccuracy for Each Patient") + 
  theme(plot.title=element_text(hjust=0.5), axis.text.x=element_text(angle=45, hjust=1, vjust=1), axis.title.x=element_blank()) #+ 
  #scale_y_continuous(limits=c(-1,1)) #+ geom_hline(yintercept=0,linetype="dashed", color="grey")
ggsave("Lasry_vs_Sorted_ptCorr.pdf", width=3, height=4)

# redo with metadata
meta.df <- BeatAML$meta
#meta.df[,c("SampleID.full","SampleID.abbrev", "Plex","Channel","Loading.Mass")] <- NULL
#meta.df[,c("Diagnosis","Recurrence","")] <- FALSE
meta.df <- meta.df[,c("Barcode.ID","FLT3.ITD","InitialAMLDiagnosis","PostChemotherapy")]
meta.long <- reshape2::melt(meta.df, id.var="Barcode.ID")
pt.r.meta <- merge(evalResults$pt.corr, meta.long, by="Barcode.ID")
ggplot(pt.r.meta, aes(x=value, y=Pearson.est)) + geom_violin() + 
  geom_boxplot() + geom_point() + theme_classic() + facet_grid(variable ~ Signature) +
  labs(y="Pearson r", title="Drug Sensitivity Prediction Accuracy for Each Patient") + 
  theme(plot.title=element_text(hjust=0.5), axis.text.x=element_text(angle=45, hjust=1, vjust=1), axis.title.x=element_blank()) #+ 
  #scale_y_continuous(limits=c(-1,1))  #+ geom_hline(yintercept=0,linetype="dashed", color="grey")
ggsave("Lasry_vs_Sorted_ptCorr_meta.pdf", width=5, height=6)

ggplot(pt.r.meta, aes(x=value, y=Pearson.est)) + geom_violin() + 
  geom_boxplot() + geom_point() + theme_classic() + facet_grid(. ~ variable) +
  labs(y="Pearson r", title="Drug Sensitivity Prediction Accuracy for Each Patient") + 
  theme(plot.title=element_text(hjust=0.5), axis.text.x=element_text(angle=45, hjust=1, vjust=1), axis.title.x=element_blank())# + scale_y_continuous(limits=c(-1,1))
ggsave("ptCorr_meta.pdf", width=5, height=4)
flt3.test <- t.test(pt.r.meta[pt.r.meta$variable == "FLT3.ITD" & pt.r.meta$value == "FALSE",]$Pearson.est, # mean 0.7583706
                    pt.r.meta[pt.r.meta$variable == "FLT3.ITD" & pt.r.meta$value == "TRUE",]$Pearson.est,  # mean 0.7243523
                    alternative = "greater"
)
flt3.test$p.value
flt3.test$estimate
# FALSE greater than TRUE p = 0.02509605 ***SIGNIFICANT***

initial.test <- t.test(pt.r.meta[pt.r.meta$variable == "InitialAMLDiagnosis" & pt.r.meta$value == "FALSE",]$Pearson.est, # mean 0.7260607
                       pt.r.meta[pt.r.meta$variable == "InitialAMLDiagnosis" & pt.r.meta$value == "TRUE",]$Pearson.est, # mean 0.7562088 
                       alternative = "less"
)
initial.test$p.value
initial.test$estimate
# less p =  0.05956627

chemo.test <- t.test(pt.r.meta[pt.r.meta$variable == "PostChemotherapy" & pt.r.meta$value == "FALSE",]$Pearson.est, # mean 0.7423232
                     pt.r.meta[pt.r.meta$variable == "PostChemotherapy" & pt.r.meta$value == "TRUE",]$Pearson.est, # mean 0.7532672 
                     alternative = "less"
)
chemo.test$p.value
chemo.test$estimate
# less p = 0.2439914

sigs.tested <- unique(pt.r.meta$Signature)
for (i in sigs.tested) {
  flt3.test <- t.test(pt.r.meta[pt.r.meta$variable == "FLT3.ITD" & pt.r.meta$value == "FALSE" & pt.r.meta$Signature==i,]$Pearson.est, # mean 0.755
                      pt.r.meta[pt.r.meta$variable == "FLT3.ITD" & pt.r.meta$value == "TRUE" & pt.r.meta$Signature==i,]$Pearson.est,  # mean 0.722
                      alternative = "greater"
  )
  if (flt3.test$p.value <= 0.05) {
    print(i,"FLT3.ITD greater pVal",flt3.test$p.value,"\n")
  }
  
  initial.test <- t.test(pt.r.meta[pt.r.meta$variable == "InitialAMLDiagnosis" & pt.r.meta$value == "FALSE" & pt.r.meta$Signature==i,]$Pearson.est, # mean 0.755
                         pt.r.meta[pt.r.meta$variable == "InitialAMLDiagnosis" & pt.r.meta$value == "TRUE" & pt.r.meta$Signature==i,]$Pearson.est,  # mean 0.722
                         alternative = "less"
  )
  if (initial.test$p.value <= 0.05) {
    print(i,"initialAMLDiagnosis less pVal",initial.test$p.value,"\n")
  }
  
  chemo.test <- t.test(pt.r.meta[pt.r.meta$variable == "PostChemotherapy" & pt.r.meta$value == "FALSE" & pt.r.meta$Signature==i,]$Pearson.est, # mean 0.755
                       pt.r.meta[pt.r.meta$variable == "PostChemotherapy" & pt.r.meta$value == "TRUE" & pt.r.meta$Signature==i,]$Pearson.est,  # mean 0.722
                       alternative = "less"
  )
  if (chemo.test$p.value <= 0.05) {
    print(i,"PostChemotherapy less pVal",chemo.test$p.value,"\n")
  }
}
# none pass p <=0.05

#### now just Ven
ven.corr <- evalResults$drug[evalResults$drug$Drug == "Venetoclax" & !is.na(evalResults$drug$delta_AUC_squared),]
pt.r <- reshape2::dcast(ven.corr, Barcode.ID ~ Signature, value.var="delta_AUC_squared")

multi.cor <- cor.test(pt.r$`Lasry: 46 genes`, pt.r$`Sorted: 25 proteins`)
Pearson.est <- multi.cor$estimate
Pearson.p <- multi.cor$p.value
stats_pearson <- substitute(
  r == est * "," ~ ~"p" ~ "=" ~ p,
  list(
    est = format(as.numeric(Pearson.est), digits = 3),
    p = format(Pearson.p, digits = 3)
  )
)

minVal <- 0
maxVal <- 0
for (i in 2:ncol(pt.r)) {
  temp.min <- min(pt.r[,i])
  temp.max <- max(pt.r[,i])
  if (temp.min < minVal) {
    minVal <- temp.min
  }
  if (temp.max > maxVal) {
    maxVal <- temp.max
  }
}
# minVal ends up as 0, maxVal ends up as 37965.320698934
ggplot(pt.r, aes(x=`Lasry: 46 genes`, y=`Sorted: 25 proteins`)) + geom_point() + theme_minimal() + 
  scale_x_continuous(limits=c(0,maxVal)) + scale_y_continuous(limits=c(0,maxVal)) + 
  labs(x="Lasry RNA-seq: 46 Genes", y = "Sorted Proteomics: 25 proteins", title="Ven Sensitivity Prediction Accuracy\nfor Each Patient (SSE)") +
  geom_smooth(method="lm",se=FALSE, linetype="dashed") + ggrepel::geom_label_repel(aes(label=Barcode.ID)) + theme(plot.title=element_text(hjust=0.5)) +
  ggplot2::geom_text(
    x = Inf, y = -Inf, vjust = "inward", hjust = "inward",
    colour = "blue", parse = TRUE,
    label = as.character(as.expression(stats_pearson)), size = 5
  ) 
ggsave("Lasry_vs_Sorted_multipleGenes_venSSE.pdf", width=4, height=4)

ggplot(ven.corr, aes(x=Signature, y=delta_AUC_squared)) + geom_violin() + 
  geom_boxplot() + geom_point() + theme_classic() + 
  labs(y="SSE", title="Ven Sensitivity Prediction\nAccuracy for Each Patient") + 
  theme(plot.title=element_text(hjust=0.5), axis.text.x=element_text(angle=45, hjust=1, vjust=1), axis.title.x=element_blank())# + scale_y_continuous(limits=c(-1,1))
ggsave("Lasry_vs_Sorted_venSSE.pdf", width=3, height=4)

multi.test <- t.test(ven.corr[ven.corr$Signature == "Lasry: 46 genes",]$delta_AUC_squared, # mean 2718.148
                     ven.corr[ven.corr$Signature == "Sorted: 25 proteins",]$delta_AUC_squared, # mean 2313.187
                     alternative = "greater"
)
multi.test$p.value
multi.test$estimate
# greater p = 0.2146595; two-sided p = 0.4293191; less p = 0.7853405

# redo with metadata
meta.df <- BeatAML$meta
#meta.df[,c("SampleID.full","SampleID.abbrev", "Plex","Channel","Loading.Mass")] <- NULL
#meta.df[,c("Diagnosis","Recurrence","")] <- FALSE
meta.df <- meta.df[,c("Barcode.ID","FLT3.ITD","InitialAMLDiagnosis","PostChemotherapy")]
meta.long <- reshape2::melt(meta.df, id.var="Barcode.ID")
pt.r.meta <- merge(ven.corr, meta.long, by="Barcode.ID")
ggplot(pt.r.meta, aes(x=value, y=delta_AUC_squared)) + geom_violin() + 
  geom_boxplot() + geom_point() + theme_classic() + facet_grid(variable ~ Signature) +
  labs(y="SSE", title="Ven Sensitivity Prediction Accuracy for Each Patient") + 
  theme(plot.title=element_text(hjust=0.5), axis.text.x=element_text(angle=45, hjust=1, vjust=1), axis.title.x=element_blank())# + scale_y_continuous(limits=c(-1,1))
ggsave("Lasry_vs_Sorted_venSSE_meta.pdf", width=5, height=6)

ggplot(pt.r.meta, aes(x=value, y=delta_AUC_squared)) + geom_violin() + 
  geom_boxplot() + geom_point() + theme_classic() + facet_grid(. ~ variable) +
  labs(y="SSE", title="Ven Sensitivity Prediction Accuracy for Each Patient") + 
  theme(plot.title=element_text(hjust=0.5), axis.text.x=element_text(angle=45, hjust=1, vjust=1), axis.title.x=element_blank())# + scale_y_continuous(limits=c(-1,1))
ggsave("venSSE_meta.pdf", width=5, height=4)
flt3.test <- t.test(pt.r.meta[pt.r.meta$variable == "FLT3.ITD" & pt.r.meta$value == "FALSE",]$delta_AUC_squared, # mean 2550.963
                    pt.r.meta[pt.r.meta$variable == "FLT3.ITD" & pt.r.meta$value == "TRUE",]$delta_AUC_squared,  # mean 2440.551 
                    alternative = "less"
)
flt3.test$p.value
flt3.test$estimate
# FALSE greater than TRUE p = 0.428956; less p = 0.571044

initial.test <- t.test(pt.r.meta[pt.r.meta$variable == "InitialAMLDiagnosis" & pt.r.meta$value == "FALSE",]$delta_AUC_squared, # mean 2578.667
                       pt.r.meta[pt.r.meta$variable == "InitialAMLDiagnosis" & pt.r.meta$value == "TRUE",]$delta_AUC_squared, # mean 2487.168
                       alternative = "greater"
)
initial.test$p.value
initial.test$estimate
# less p = 0.5579233; greater p = 0.4420767

chemo.test <- t.test(pt.r.meta[pt.r.meta$variable == "PostChemotherapy" & pt.r.meta$value == "FALSE",]$delta_AUC_squared, # mean 2162.294
                     pt.r.meta[pt.r.meta$variable == "PostChemotherapy" & pt.r.meta$value == "TRUE",]$delta_AUC_squared, # mean 3394.053
                     alternative = "less"
)
chemo.test$p.value
chemo.test$estimate
# greater p = 0.9477143; less p = 0.05228567

chemo.test.s <- t.test(pt.r.meta[pt.r.meta$variable == "PostChemotherapy" & pt.r.meta$value == "FALSE" & pt.r.meta$Signature=="Sorted: 25 proteins",]$delta_AUC_squared, # mean 2019.436 
                        pt.r.meta[pt.r.meta$variable == "PostChemotherapy" & pt.r.meta$value == "TRUE" & pt.r.meta$Signature=="Sorted: 25 proteins",]$delta_AUC_squared, # mean 3043.368 
                        alternative = "less"
)
chemo.test.s$p.value
chemo.test.s$estimate
# less p = 0.1475887

chemo.test.l <- t.test(pt.r.meta[pt.r.meta$variable == "PostChemotherapy" & pt.r.meta$value == "FALSE" & pt.r.meta$Signature=="Lasry: 46 genes",]$delta_AUC_squared, # mean 2305.152
                         pt.r.meta[pt.r.meta$variable == "PostChemotherapy" & pt.r.meta$value == "TRUE" & pt.r.meta$Signature=="Lasry: 46 genes",]$delta_AUC_squared, # mean 3744.737 
                         alternative = "less"
)
chemo.test.l$p.value
chemo.test.l$estimate
# less p = 0.110846

sigs.tested <- unique(pt.r.meta$Signature)
for (i in sigs.tested) {
  flt3.test <- t.test(pt.r.meta[pt.r.meta$variable == "FLT3.ITD" & pt.r.meta$value == "FALSE" & pt.r.meta$Signature==i,]$delta_AUC_squared, # mean 0.755
                      pt.r.meta[pt.r.meta$variable == "FLT3.ITD" & pt.r.meta$value == "TRUE" & pt.r.meta$Signature==i,]$delta_AUC_squared,  # mean 0.722
                      alternative = "greater"
  )
  if (flt3.test$p.value <= 0.05) {
    print(i,"FLT3.ITD greater pVal",flt3.test$p.value,"\n")
  }
  
  initial.test <- t.test(pt.r.meta[pt.r.meta$variable == "InitialAMLDiagnosis" & pt.r.meta$value == "FALSE" & pt.r.meta$Signature==i,]$delta_AUC_squared, # mean 0.755
                         pt.r.meta[pt.r.meta$variable == "InitialAMLDiagnosis" & pt.r.meta$value == "TRUE" & pt.r.meta$Signature==i,]$delta_AUC_squared,  # mean 0.722
                         alternative = "less"
  )
  if (initial.test$p.value <= 0.05) {
    print(i,"initialAMLDiagnosis less pVal",initial.test$p.value,"\n")
  }
  
  chemo.test <- t.test(pt.r.meta[pt.r.meta$variable == "PostChemotherapy" & pt.r.meta$value == "FALSE" & pt.r.meta$Signature==i,]$delta_AUC_squared, # mean 0.755
                       pt.r.meta[pt.r.meta$variable == "PostChemotherapy" & pt.r.meta$value == "TRUE" & pt.r.meta$Signature==i,]$delta_AUC_squared,  # mean 0.722
                       alternative = "less"
  )
  if (chemo.test$p.value <= 0.05) {
    print(i,"PostChemotherapy less pVal",chemo.test$p.value,"\n")
  }
}
# none passed p <= 0.05
# 
# base.path <- "~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/data/clinical_metadata/"
# setwd(base.path)
# patient.key <- readxl::read_xlsx("lab_dbgap_key.xlsx")
# patient.meta <- readxl::read_xlsx("Table_S1.xlsx",sheet=1)
# patient.meta <- as.data.frame(patient.meta)
# #rownames(patient.meta) <- patient.meta[,1]
# 
# # does metadata correlate Ven response?
# key.cols <- colnames(patient.meta)[colnames(patient.meta) %in% colnames(patient.key)]
# num.cols <- colnames(dplyr::select_if(patient.meta, is.numeric))
# num.meta <- dplyr::distinct(patient.meta[,c(key.cols, num.cols)])
# num.meta <- merge(patient.key, num.meta, by = key.cols)
# num.meta <- num.meta[,c("labId", num.cols[2:length(num.cols)])] # leave out dbgap_subject_id from numeric columns
# 
# # correlation between NCF2 protein or LRRC25 RNA and Ven sensitivity


# hist(evalResults$pt.corr$Pearson.est)
# med.pt.r <- median(evalResults$pt.corr$Pearson.est) # 0.785882478840521
# sd.pt.r <- sd(evalResults$pt.corr$Pearson.est) # 0.171287467988083
# frac.sd <- sd.pt.r/med.pt.r # 0.21795557554712
# pt.r <- reshape2::dcast(evalResults$pt.corr, Barcode.ID ~ Signature, value.var="Pearson.est")
# 
# multi.cor <- cor.test(pt.r$Lasry_16_Genes, pt.r$Sorted_4_Proteins)
# Pearson.est <- multi.cor$estimate
# Pearson.p <- multi.cor$p.value
# stats_pearson <- substitute(
#   r == est * "," ~ ~"p" ~ "=" ~ p,
#   list(
#     est = format(as.numeric(Pearson.est), digits = 3),
#     p = format(Pearson.p, digits = 3)
#   )
# )
# ggplot(pt.r, aes(x=Lasry_16_Genes, y=Sorted_4_Proteins)) + geom_point() + theme_minimal() + 
#   scale_x_continuous(limits=c(-1,1)) + scale_y_continuous(limits=c(-1,1)) + 
#   labs(x="Lasry RNA-seq: 16 Genes", y = "Sorted Proteomics: 4 Proteins", title="Drug Sensitivity Prediction Accuracy\nfor Each Patient (Pearson Correlation)") +
#   geom_smooth(se=FALSE, linetype="dashed") + ggrepel::geom_label_repel(aes(label=Barcode.ID)) + theme(plot.title=element_text(hjust=0.5)) +
#   ggplot2::geom_text(
#     x = Inf, y = -Inf, vjust = "inward", hjust = "inward",
#     colour = "blue", parse = TRUE,
#     label = as.character(as.expression(stats_pearson)), size = 5
#   ) 
# ggsave("Lasry_vs_Sorted_multipleGenes_ptCorr.pdf", width=4, height=4)
# 
# single.cor <- cor.test(pt.r$Lasry_1_Gene, pt.r$Sorted_1_Protein)
# Pearson.est <- single.cor$estimate
# Pearson.p <- single.cor$p.value
# stats_pearson <- substitute(
#   r == est * "," ~ ~"p" ~ "=" ~ p,
#   list(
#     est = format(as.numeric(Pearson.est), digits = 3),
#     p = format(Pearson.p, digits = 3)
#   )
# )
# ggplot(pt.r, aes(x=Lasry_1_Gene, y=Sorted_1_Protein)) + geom_point() + theme_minimal() + 
#   scale_x_continuous(limits=c(-1,1)) + scale_y_continuous(limits=c(-1,1)) + 
#   labs(x="Lasry RNA-seq: LRCC25", y = "Sorted Proteomics: NCF2", title="Drug Sensitivity Prediction Accuracy\nfor Each Patient (Pearson Correlation)") +
#   geom_smooth(se=FALSE, linetype="dashed") + ggrepel::geom_label_repel(aes(label=Barcode.ID)) + theme(plot.title=element_text(hjust=0.5)) +
#   ggplot2::geom_text(
#     x = Inf, y = -Inf, vjust = "inward", hjust = "inward",
#     colour = "blue", parse = TRUE,
#     label = as.character(as.expression(stats_pearson)), size = 5
#   ) 
# ggsave("Lasry_vs_Sorted_1Gene_ptCorr.pdf", width=4, height=4)
# 
# single.test <- t.test(evalResults$pt.corr[evalResults$pt.corr$Signature == "Lasry_1_Gene",]$Pearson.est, # mean 0.741
#                       evalResults$pt.corr[evalResults$pt.corr$Signature == "Sorted_1_Protein",]$Pearson.est, # mean 0.738
#                       #alternative = "greater"
# )
# single.test$estimate
# single.test$p.value
# # greater p = 0.445; two-sided p = 0.889; less p = 0.555
# multi.test <- t.test(evalResults$pt.corr[evalResults$pt.corr$Signature == "Lasry_16_Genes",]$Pearson.est, # mean 0.745
#                      evalResults$pt.corr[evalResults$pt.corr$Signature == "Sorted_4_Proteins",]$Pearson.est, # mean 0.747
#                      #alternative = "less"
# )
# multi.test$estimate
# multi.test$p.value
# # greater p = 0.544; two-sided p = 0.912; less p = 0.456
# 
# sorted.test <- t.test(evalResults$pt.corr[evalResults$pt.corr$Signature == "Sorted_1_Protein",]$Pearson.est, # mean 0.738
#                       evalResults$pt.corr[evalResults$pt.corr$Signature == "Sorted_4_Proteins",]$Pearson.est, # mean 0.747
#                       alternative = "less"
# )
# sorted.test$p.value
# sorted.test$estimate
# # greater p = 0.687; two-sided p = 0.626; less p = 0.313
# lasry.test <- t.test(evalResults$pt.corr[evalResults$pt.corr$Signature == "Lasry_1_Gene",]$Pearson.est, # mean 0.741
#                      evalResults$pt.corr[evalResults$pt.corr$Signature == "Lasry_16_Genes",]$Pearson.est, # mean 0.745
#                      alternative = "greater"
# )
# lasry.test$p.value
# lasry.test$estimate
# # greater p = 0.599; two-sided p = 0.803; less p = 0.4015
# 
# n.test <- t.test(evalResults$pt.corr[evalResults$pt.corr$Signature %in% c("Lasry_1_Gene","Sorted_1_Protein"),]$Pearson.est, # mean 0.740
#                  evalResults$pt.corr[evalResults$pt.corr$Signature %in% c("Lasry_16_Genes","Sorted_4_Proteins"),]$Pearson.est, # mean 0.746
#                  #alternative = "less"
# )
# n.test$p.value
# n.test$estimate
# # greater p = 0.700; two-sided p = 0.5996472; less p = 0.2998236
# 
# ggplot(evalResults$pt.corr, aes(x=Signature, y=Pearson.est)) + geom_violin() + 
#   geom_boxplot() + geom_point() + theme_classic() + 
#   labs(y="Pearson r", title="Drug Sensitivity Prediction\nAccuracy for Each Patient") + 
#   theme(plot.title=element_text(hjust=0.5), axis.text.x=element_text(angle=45, hjust=1, vjust=1))# + scale_y_continuous(limits=c(-1,1))
# ggsave("Lasry_vs_Sorted_ptCorr.pdf", width=4, height=4)
# 
# # redo with metadata
# meta.df <- BeatAML$meta
# #meta.df[,c("SampleID.full","SampleID.abbrev", "Plex","Channel","Loading.Mass")] <- NULL
# #meta.df[,c("Diagnosis","Recurrence","")] <- FALSE
# meta.df <- meta.df[,c("Barcode.ID","FLT3.ITD","InitialAMLDiagnosis","PostChemotherapy")]
# meta.long <- reshape2::melt(meta.df, id.var="Barcode.ID")
# pt.r.meta <- merge(evalResults$pt.corr, meta.long, by="Barcode.ID")
# ggplot(pt.r.meta, aes(x=value, y=Pearson.est)) + geom_violin() + 
#   geom_boxplot() + geom_point() + theme_classic() + facet_grid(variable ~ Signature) +
#   labs(y="Pearson r", title="Drug Sensitivity Prediction Accuracy for Each Patient") + 
#   theme(plot.title=element_text(hjust=0.5), axis.text.x=element_text(angle=45, hjust=1, vjust=1), axis.title.x=element_blank())# + scale_y_continuous(limits=c(-1,1))
# ggsave("Lasry_vs_Sorted_ptCorr_meta.pdf", width=6, height=6)
# 
# ggplot(pt.r.meta, aes(x=value, y=Pearson.est)) + geom_violin() + 
#   geom_boxplot() + geom_point() + theme_classic() + facet_grid(. ~ variable) +
#   labs(y="Pearson r", title="Drug Sensitivity Prediction Accuracy for Each Patient") + 
#   theme(plot.title=element_text(hjust=0.5), axis.text.x=element_text(angle=45, hjust=1, vjust=1), axis.title.x=element_blank())# + scale_y_continuous(limits=c(-1,1))
# ggsave("ptCorr_meta.pdf", width=5, height=4)
# flt3.test <- t.test(pt.r.meta[pt.r.meta$variable == "FLT3.ITD" & pt.r.meta$value == "FALSE",]$Pearson.est, # mean 0.755
#                     pt.r.meta[pt.r.meta$variable == "FLT3.ITD" & pt.r.meta$value == "TRUE",]$Pearson.est,  # mean 0.722
#                     #alternative = "less"
#                     )
# flt3.test$p.value
# flt3.test$estimate
# # two-sided p = 0.009206772; FALSE greater than TRUE p = 0.004603386; less p = 0.9953966 ***SIGNIFICANT***
# 
# initial.test <- t.test(pt.r.meta[pt.r.meta$variable == "InitialAMLDiagnosis" & pt.r.meta$value == "FALSE",]$Pearson.est, # mean 0.722
#                     pt.r.meta[pt.r.meta$variable == "InitialAMLDiagnosis" & pt.r.meta$value == "TRUE",]$Pearson.est, # mean 0.754
#                     #alternative = "greater"
#                     )
# initial.test$p.value
# initial.test$estimate
# # two-sided p = 0.02176726; less p = 0.01088363; greater p = 0.9891164 ***SIGNIFICANT***
# 
# chemo.test <- t.test(pt.r.meta[pt.r.meta$variable == "PostChemotherapy" & pt.r.meta$value == "FALSE",]$Pearson.est, # mean 0.740
#                     pt.r.meta[pt.r.meta$variable == "PostChemotherapy" & pt.r.meta$value == "TRUE",]$Pearson.est, # mean 0.750
#                     #alternative = "greater"
#                     )
# chemo.test$p.value
# chemo.test$estimate
# # two-sided p = 0.3898602; greater p = 0.8050699; less p = 0.1949301
# 
# sigs.tested <- unique(pt.r.meta$Signature)
# for (i in sigs.tested) {
#   flt3.test <- t.test(pt.r.meta[pt.r.meta$variable == "FLT3.ITD" & pt.r.meta$value == "FALSE" & pt.r.meta$Signature==i,]$Pearson.est, # mean 0.755
#                       pt.r.meta[pt.r.meta$variable == "FLT3.ITD" & pt.r.meta$value == "TRUE" & pt.r.meta$Signature==i,]$Pearson.est,  # mean 0.722
#                       alternative = "greater"
#   )
#   if (flt3.test$p.value <= 0.05) {
#     print(i,"FLT3.ITD greater pVal",flt3.test$p.value,"\n")
#   }
#   
#   initial.test <- t.test(pt.r.meta[pt.r.meta$variable == "InitialAMLDiagnosis" & pt.r.meta$value == "FALSE" & pt.r.meta$Signature==i,]$Pearson.est, # mean 0.755
#                          pt.r.meta[pt.r.meta$variable == "InitialAMLDiagnosis" & pt.r.meta$value == "TRUE" & pt.r.meta$Signature==i,]$Pearson.est,  # mean 0.722
#                          alternative = "less"
#   )
#   if (initial.test$p.value <= 0.05) {
#     print(i,"initialAMLDiagnosis less pVal",initial.test$p.value,"\n")
#   }
#   
#   chemo.test <- t.test(pt.r.meta[pt.r.meta$variable == "PostChemotherapy" & pt.r.meta$value == "FALSE" & pt.r.meta$Signature==i,]$Pearson.est, # mean 0.755
#                          pt.r.meta[pt.r.meta$variable == "PostChemotherapy" & pt.r.meta$value == "TRUE" & pt.r.meta$Signature==i,]$Pearson.est,  # mean 0.722
#                          alternative = "less"
#   )
#   if (chemo.test$p.value <= 0.05) {
#     print(i,"PostChemotherapy less pVal",chemo.test$p.value,"\n")
#   }
# }
# # none pass p <=0.05
# 
# #### now just Ven
# ven.corr <- evalResults$drug[evalResults$drug$Drug == "Venetoclax" & !is.na(evalResults$drug$delta_AUC_squared),]
# pt.r <- reshape2::dcast(ven.corr, Barcode.ID ~ Signature, value.var="delta_AUC_squared")
# 
# multi.cor <- cor.test(pt.r$Lasry_16_Genes, pt.r$Sorted_4_Proteins)
# Pearson.est <- multi.cor$estimate
# Pearson.p <- multi.cor$p.value
# stats_pearson <- substitute(
#   r == est * "," ~ ~"p" ~ "=" ~ p,
#   list(
#     est = format(as.numeric(Pearson.est), digits = 3),
#     p = format(Pearson.p, digits = 3)
#   )
# )
# 
# minVal <- 0
# maxVal <- 0
# for (i in 2:ncol(pt.r)) {
#   temp.min <- min(pt.r[,i])
#   temp.max <- max(pt.r[,i])
#   if (temp.min < minVal) {
#     minVal <- temp.min
#   }
#   if (temp.max > maxVal) {
#     maxVal <- temp.max
#   }
# }
# # minVal ends up as 0, maxVal ends up as 37965.320698934
# ggplot(pt.r, aes(x=Lasry_16_Genes, y=Sorted_4_Proteins)) + geom_point() + theme_minimal() + 
#   scale_x_continuous(limits=c(0,maxVal)) + scale_y_continuous(limits=c(0,maxVal)) + 
#   labs(x="Lasry RNA-seq: 16 Genes", y = "Sorted Proteomics: 4 Proteins", title="Ven Sensitivity Prediction Accuracy\nfor Each Patient (SSE)") +
#   geom_smooth(method="lm",se=FALSE, linetype="dashed") + ggrepel::geom_label_repel(aes(label=Barcode.ID)) + theme(plot.title=element_text(hjust=0.5)) +
#   ggplot2::geom_text(
#     x = Inf, y = -Inf, vjust = "inward", hjust = "inward",
#     colour = "blue", parse = TRUE,
#     label = as.character(as.expression(stats_pearson)), size = 5
#   ) 
# ggsave("Lasry_vs_Sorted_multipleGenes_venSSE.pdf", width=4, height=4)
# 
# single.cor <- cor.test(pt.r$Lasry_1_Gene, pt.r$Sorted_1_Protein)
# Pearson.est <- single.cor$estimate
# Pearson.p <- single.cor$p.value
# stats_pearson <- substitute(
#   r == est * "," ~ ~"p" ~ "=" ~ p,
#   list(
#     est = format(as.numeric(Pearson.est), digits = 3),
#     p = format(Pearson.p, digits = 3)
#   )
# )
# ggplot(pt.r, aes(x=Lasry_1_Gene, y=Sorted_1_Protein)) + geom_point() + theme_minimal() + 
#   scale_x_continuous(limits=c(0, maxVal)) + scale_y_continuous(limits=c(0,maxVal)) + 
#   labs(x="Lasry RNA-seq: LRCC25", y = "Sorted Proteomics: NCF2", title="Ven Sensitivity Prediction Accuracy\nfor Each Patient (SSE)") +
#   geom_smooth(method="lm", se=FALSE, linetype="dashed") + ggrepel::geom_label_repel(aes(label=Barcode.ID)) + theme(plot.title=element_text(hjust=0.5)) +
#   ggplot2::geom_text(
#     x = Inf, y = -Inf, vjust = "inward", hjust = "inward",
#     colour = "blue", parse = TRUE,
#     label = as.character(as.expression(stats_pearson)), size = 5
#   ) 
# ggsave("Lasry_vs_Sorted_1Gene_venSSE.pdf", width=4, height=4)
# 
# ggplot(ven.corr, aes(x=Signature, y=delta_AUC_squared)) + geom_violin() + 
#   geom_boxplot() + geom_point() + theme_classic() + 
#   labs(y="SSE", title="Ven Sensitivity Prediction\nAccuracy for Each Patient") + 
#   theme(plot.title=element_text(hjust=0.5), axis.text.x=element_text(angle=45, hjust=1, vjust=1))# + scale_y_continuous(limits=c(-1,1))
# ggsave("Lasry_vs_Sorted_venSSE.pdf", width=4, height=4)
# single.test <- t.test(ven.corr[ven.corr$Signature == "Lasry_1_Gene",]$delta_AUC_squared, # mean 3359.415
#                       ven.corr[ven.corr$Signature == "Sorted_1_Protein",]$delta_AUC_squared, # mean 2937.280
#                       alternative = "less"
#                       )
# single.test$p.value
# single.test$estimate
# # greater p = 0.2453893; two-sided p = 0.491; less p = 0.755
# multi.test <- t.test(ven.corr[ven.corr$Signature == "Lasry_16_Genes",]$delta_AUC_squared, # mean 2781.937
#                       ven.corr[ven.corr$Signature == "Sorted_4_Proteins",]$delta_AUC_squared, # mean 2302.425
#                       alternative = "greater"
#                      )
# multi.test$p.value
# multi.test$estimate
# # greater p = 0.158; two-sided p = 0.316; less p = 0.842
# 
# sorted.test <- t.test(ven.corr[ven.corr$Signature == "Sorted_1_Protein",]$delta_AUC_squared, # mean 2937.28
#                       ven.corr[ven.corr$Signature == "Sorted_4_Proteins",]$delta_AUC_squared, # mean 2302.425
#                       #alternative = "greater"
# )
# sorted.test$p.value
# sorted.test$estimate
# # greater p = 0.0972; two-sided p = 0.194; less p = 0.903
# lasry.test <- t.test(ven.corr[ven.corr$Signature == "Lasry_1_Gene",]$delta_AUC_squared, # mean 3359.415
#                      ven.corr[ven.corr$Signature == "Lasry_16_Genes",]$delta_AUC_squared, # mean 2781.937
#                      #alternative = "less"
# )
# lasry.test$p.value
# lasry.test$estimate
# # greater p = 0.1696971; two-sided p = 0.3393942; less p = 0.8303029
# 
# n.test <- t.test(ven.corr[ven.corr$Signature %in% c("Lasry_1_Gene","Sorted_1_Protein"),]$delta_AUC_squared, # mean 3148.348
#                      ven.corr[ven.corr$Signature %in% c("Lasry_16_Genes","Sorted_4_Proteins"),]$delta_AUC_squared, # mean 2542.181
#                      alternative = "greater"
# )
# n.test$p.value
# n.test$estimate
# # greater p = 0.05926081; two-sided p = 0.1185216; less p = 0.941
# 
# # redo with metadata
# meta.df <- BeatAML$meta
# #meta.df[,c("SampleID.full","SampleID.abbrev", "Plex","Channel","Loading.Mass")] <- NULL
# #meta.df[,c("Diagnosis","Recurrence","")] <- FALSE
# meta.df <- meta.df[,c("Barcode.ID","FLT3.ITD","InitialAMLDiagnosis","PostChemotherapy")]
# meta.long <- reshape2::melt(meta.df, id.var="Barcode.ID")
# pt.r.meta <- merge(ven.corr, meta.long, by="Barcode.ID")
# ggplot(pt.r.meta, aes(x=value, y=delta_AUC_squared)) + geom_violin() + 
#   geom_boxplot() + geom_point() + theme_classic() + facet_grid(variable ~ Signature) +
#   labs(y="SSE", title="Ven Sensitivity Prediction Accuracy for Each Patient") + 
#   theme(plot.title=element_text(hjust=0.5), axis.text.x=element_text(angle=45, hjust=1, vjust=1), axis.title.x=element_blank())# + scale_y_continuous(limits=c(-1,1))
# ggsave("Lasry_vs_Sorted_venSSE_meta.pdf", width=6, height=6)
# 
# ggplot(pt.r.meta, aes(x=value, y=delta_AUC_squared)) + geom_violin() + 
#   geom_boxplot() + geom_point() + theme_classic() + facet_grid(. ~ variable) +
#   labs(y="SSE", title="Ven Sensitivity Prediction Accuracy for Each Patient") + 
#   theme(plot.title=element_text(hjust=0.5), axis.text.x=element_text(angle=45, hjust=1, vjust=1), axis.title.x=element_blank())# + scale_y_continuous(limits=c(-1,1))
# ggsave("venSSE_meta.pdf", width=5, height=4)
# flt3.test <- t.test(pt.r.meta[pt.r.meta$variable == "FLT3.ITD" & pt.r.meta$value == "FALSE",]$delta_AUC_squared, # mean 2872.144
#                     pt.r.meta[pt.r.meta$variable == "FLT3.ITD" & pt.r.meta$value == "TRUE",]$delta_AUC_squared,  # mean 2841.264
#                     alternative = "greater"
# )
# flt3.test$p.value
# flt3.test$estimate
# # two-sided p = 0.9989766; FALSE greater than TRUE p = 0.494883; less p = 0.505117
# 
# initial.test <- t.test(pt.r.meta[pt.r.meta$variable == "InitialAMLDiagnosis" & pt.r.meta$value == "FALSE",]$delta_AUC_squared, # mean 2994.937
#                        pt.r.meta[pt.r.meta$variable == "InitialAMLDiagnosis" & pt.r.meta$value == "TRUE",]$delta_AUC_squared, # mean 2777.555
#                        alternative = "less"
# )
# initial.test$p.value
# initial.test$estimate
# # two-sided p = 0.6434979; less p = 0.6782511; greater p = 0.3217489
# 
# chemo.test <- t.test(pt.r.meta[pt.r.meta$variable == "PostChemotherapy" & pt.r.meta$value == "FALSE",]$delta_AUC_squared, # mean 2463.747; 348 data points
#                      pt.r.meta[pt.r.meta$variable == "PostChemotherapy" & pt.r.meta$value == "TRUE",]$delta_AUC_squared, # mean 3793.606; 140 data points
#                      alternative = "less"
# )
# chemo.test$p.value
# chemo.test$estimate
# # two-sided p = 0.01722762; greater p = 0.991; less p = 0.008613808 ***SIGNIFICANT***
# 
# chemo.test.s1 <- t.test(pt.r.meta[pt.r.meta$variable == "PostChemotherapy" & pt.r.meta$value == "FALSE" & pt.r.meta$Signature=="Sorted_1_Protein",]$delta_AUC_squared, # mean 2563
#                      pt.r.meta[pt.r.meta$variable == "PostChemotherapy" & pt.r.meta$value == "TRUE" & pt.r.meta$Signature=="Sorted_1_Protein",]$delta_AUC_squared, # mean 3867
#                      alternative = "less"
# )
# chemo.test.s1$p.value
# chemo.test.s1$estimate
# # less p = 0.11
# 
# chemo.test.s4 <- t.test(pt.r.meta[pt.r.meta$variable == "PostChemotherapy" & pt.r.meta$value == "FALSE" & pt.r.meta$Signature=="Sorted_4_Proteins",]$delta_AUC_squared, # mean 2046
#                         pt.r.meta[pt.r.meta$variable == "PostChemotherapy" & pt.r.meta$value == "TRUE" & pt.r.meta$Signature=="Sorted_4_Proteins",]$delta_AUC_squared, # mean 2940
#                         alternative = "less"
# )
# chemo.test.s4$p.value
# chemo.test.s4$estimate
# # less p = 0.116
# 
# chemo.test.l1 <- t.test(pt.r.meta[pt.r.meta$variable == "PostChemotherapy" & pt.r.meta$value == "FALSE" & pt.r.meta$Signature=="Lasry_1_Gene",]$delta_AUC_squared, # mean 2738
#                         pt.r.meta[pt.r.meta$variable == "PostChemotherapy" & pt.r.meta$value == "TRUE" & pt.r.meta$Signature=="Lasry_1_Gene",]$delta_AUC_squared, # mean 4904
#                         alternative = "less"
# )
# chemo.test.l1$p.value
# chemo.test.l1$estimate
# # less p = 0.0645
# 
# chemo.test.l16 <- t.test(pt.r.meta[pt.r.meta$variable == "PostChemotherapy" & pt.r.meta$value == "FALSE" & pt.r.meta$Signature=="Lasry_16_Genes",]$delta_AUC_squared, # mean 2507.863
#                         pt.r.meta[pt.r.meta$variable == "PostChemotherapy" & pt.r.meta$value == "TRUE" & pt.r.meta$Signature=="Lasry_16_Genes",]$delta_AUC_squared, # mean 3463.208
#                         alternative = "less"
# )
# chemo.test.l16$p.value
# chemo.test.l16$estimate
# # less p = 0.206
# 
# sigs.tested <- unique(pt.r.meta$Signature)
# for (i in sigs.tested) {
#   flt3.test <- t.test(pt.r.meta[pt.r.meta$variable == "FLT3.ITD" & pt.r.meta$value == "FALSE" & pt.r.meta$Signature==i,]$delta_AUC_squared, # mean 0.755
#                       pt.r.meta[pt.r.meta$variable == "FLT3.ITD" & pt.r.meta$value == "TRUE" & pt.r.meta$Signature==i,]$delta_AUC_squared,  # mean 0.722
#                       alternative = "greater"
#   )
#   if (flt3.test$p.value <= 0.05) {
#     print(i,"FLT3.ITD greater pVal",flt3.test$p.value,"\n")
#   }
#   
#   initial.test <- t.test(pt.r.meta[pt.r.meta$variable == "InitialAMLDiagnosis" & pt.r.meta$value == "FALSE" & pt.r.meta$Signature==i,]$delta_AUC_squared, # mean 0.755
#                          pt.r.meta[pt.r.meta$variable == "InitialAMLDiagnosis" & pt.r.meta$value == "TRUE" & pt.r.meta$Signature==i,]$delta_AUC_squared,  # mean 0.722
#                          alternative = "less"
#   )
#   if (initial.test$p.value <= 0.05) {
#     print(i,"initialAMLDiagnosis less pVal",initial.test$p.value,"\n")
#   }
#   
#   chemo.test <- t.test(pt.r.meta[pt.r.meta$variable == "PostChemotherapy" & pt.r.meta$value == "FALSE" & pt.r.meta$Signature==i,]$delta_AUC_squared, # mean 0.755
#                        pt.r.meta[pt.r.meta$variable == "PostChemotherapy" & pt.r.meta$value == "TRUE" & pt.r.meta$Signature==i,]$delta_AUC_squared,  # mean 0.722
#                        alternative = "less"
#   )
#   if (chemo.test$p.value <= 0.05) {
#     print(i,"PostChemotherapy less pVal",chemo.test$p.value,"\n")
#   }
# }
# # none passed p <= 0.05
# # 
dia.wo.out <- readRDS("~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/analysis/DIA_2batches_noOutliers.rds")
# sorted.patients <- c("18-00105", "21-00839", "22-00571", "22-00117", "16-01184",
#                      "19-00074", "18-00103", "21-00432", "17-01060", "22-00251")
sorted.patients <- unique(dia.wo.out$meta$patient)
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
# # 
# # # correlation between NCF2 protein or LRRC25 RNA and Ven sensitivity

fab.meta <- merge(patient.key, patient.meta, by=key.cols)
fab.meta <- fab.meta[,c("labId","fabBlastMorphology")]
fab.meta <- na.omit(fab.meta[fab.meta$labId %in% sorted.patients,]) # 9
# 3 are M1, 3 are M2, 1 is M4, 2 are M5

npm1.meta <- merge(patient.key, patient.meta, by=key.cols)
write.csv(npm1.meta, "Table_S1_merged.csv", row.names=FALSE)
beat.meta <- npm1.meta[!(npm1.meta$labId %in% sorted.patients),]
write.csv(beat.meta, "Table_S1_notSorted.csv", row.names=FALSE)
beat.meta.omics <- beat.meta[beat.meta$labId %in% unique(c(BeatAML$rna$Barcode.ID, BeatAML$global$Barcode.ID)),]
beat.meta.overlap <- beat.meta[beat.meta$labId %in% BeatAML$drug[!is.na(BeatAML$drug$Venetoclax),]$Barcode.ID,] # 106
write.csv(beat.meta.overlap, "Table_S1_notSorted_rnaProtVen.csv", row.names=FALSE)

# get numeric metadata
num.meta.overlap <- dplyr::distinct(beat.meta.overlap[,c("labId", num.cols)])
num.meta.overlap[,c("centerID","dbgap_subject_id")] <- NULL
write.csv(num.meta.overlap,"Table_S1_notSorted_rnaProtVen_numeric.csv", row.names=FALSE)

# correlations with pt corr for Sorted_25
corr.input <- pt.r[,c(1,3)]
colnames(corr.input)[1] <- "labId"
corr.input <- merge(corr.input, num.meta.overlap, by="labId")
write.csv(corr.input,"Table_S1_notSorted_rnaProtVen_numeric_corrInput.csv", row.names=FALSE)
sorted.meta.corr <- DMEA::rank_corr(corr.input, variable="Metadata", value="Metadata")
# No correlations met the FDR cut-off to produce scatter plots
write.csv(sorted.meta.corr$result,"Table_S1_notSorted_rnaProtVen_numeric_corr.csv", row.names=FALSE)

corr.input <- ven.corr[ven.corr$Signature=="Sorted: 25 proteins",c("Barcode.ID","delta_AUC_squared")]
colnames(corr.input)[1] <- "labId"
corr.input <- merge(corr.input, num.meta.overlap, by="labId")
write.csv(corr.input,"Table_S1_notSorted_rnaProtVen_numeric_corrInput_venSSE.csv", row.names=FALSE)
sorted.meta.corr <- DMEA::rank_corr(corr.input, variable="Metadata", value="Metadata", plots=FALSE)
# No correlations met the FDR cut-off to produce scatter plots
write.csv(sorted.meta.corr$result,"Table_S1_notSorted_rnaProtVen_numeric_corr_venSSE.csv", row.names=FALSE)
# better prediction with higher allelic ratio (Pearson r=0.3461, p=3.665E-4, q=1.026E-2)

factor.meta.overlap <- data.frame(labId=beat.meta.overlap$labId)
factor.meta.overlap[,c("Male","White","Hispanic","CEBPA_Biallelic","Relapse",
                       "Denovo","Transformed","cumulativeChemo",
                       "priorMalignancyNonMyeloid",
                       "priorMalignancyRadiationTx","priorMDS",
                       "priorMDSMoreThanTwoMths","priorMDSMPN",
                       "priorMDSMPNMoreThanTwoMths","priorMPN",
                       "priorMPNMoreThanTwoMths","BoneMarrow",
                       "CRWithInductionChemo")] <- NA

factor.meta.overlap[factor.meta.overlap$labId %in% 
                      beat.meta.overlap[beat.meta.overlap$reportedRace=="White",]$labId,]$White <- TRUE
factor.meta.overlap[factor.meta.overlap$labId %in% 
                      beat.meta.overlap[beat.meta.overlap$reportedRace!="White" &
                                          beat.meta.overlap$reportedRace!="Unknown",]$labId,]$White <- FALSE

factor.meta.overlap[factor.meta.overlap$labId %in% 
                      beat.meta.overlap[beat.meta.overlap$consensus_sex=="Male",]$labId,]$Male <- TRUE
factor.meta.overlap[factor.meta.overlap$labId %in% 
                      beat.meta.overlap[beat.meta.overlap$consensus_sex!="Male",]$labId,]$Male <- FALSE

factor.meta.overlap[factor.meta.overlap$labId %in% 
                      beat.meta.overlap[beat.meta.overlap$reportedEthnicity=="HISPANIC",]$labId,]$Hispanic <- TRUE
factor.meta.overlap[factor.meta.overlap$labId %in% 
                      beat.meta.overlap[beat.meta.overlap$reportedEthnicity=="NON-HISPANIC",]$labId,]$Hispanic <- FALSE

factor.meta.overlap[factor.meta.overlap$labId %in% 
                      beat.meta.overlap[beat.meta.overlap$CEBPA_Biallelic=="bi",]$labId,]$CEBPA_Biallelic <- TRUE
factor.meta.overlap[factor.meta.overlap$labId %in% 
                      beat.meta.overlap[beat.meta.overlap$CEBPA_Biallelic=="mono",]$labId,]$CEBPA_Biallelic <- FALSE

factor.meta.overlap[factor.meta.overlap$labId %in% 
                      beat.meta.overlap[beat.meta.overlap$isRelapse=="TRUE",]$labId,]$Relapse <- TRUE
factor.meta.overlap[factor.meta.overlap$labId %in% 
                      beat.meta.overlap[beat.meta.overlap$isRelapse=="FALSE",]$labId,]$Relapse <- FALSE

factor.meta.overlap[factor.meta.overlap$labId %in% 
                      beat.meta.overlap[beat.meta.overlap$isDenovo=="TRUE",]$labId,]$Denovo <- TRUE
factor.meta.overlap[factor.meta.overlap$labId %in% 
                      beat.meta.overlap[beat.meta.overlap$isDenovo=="FALSE",]$labId,]$Denovo <- FALSE

factor.meta.overlap[factor.meta.overlap$labId %in% 
                      beat.meta.overlap[beat.meta.overlap$isTransformed=="TRUE",]$labId,]$Transformed <- TRUE
factor.meta.overlap[factor.meta.overlap$labId %in% 
                      beat.meta.overlap[beat.meta.overlap$isTransformed=="FALSE",]$labId,]$Transformed <- FALSE

factor.meta.overlap[factor.meta.overlap$labId %in% 
                      beat.meta.overlap[beat.meta.overlap$specimenType=="Bone Marrow Aspirate",]$labId,]$BoneMarrow <- TRUE
factor.meta.overlap[factor.meta.overlap$labId %in% 
                      beat.meta.overlap[beat.meta.overlap$specimenType!="Bone Marrow Aspirate",]$labId,]$BoneMarrow <- FALSE

factor.meta.overlap[factor.meta.overlap$labId %in% 
                      beat.meta.overlap[grepl("Complete Response",beat.meta.overlap$responseToInductionTx) &
                                          beat.meta.overlap$typeInductionTx=="Standard Chemotherapy",]$labId,]$CRWithInductionChemo <- TRUE
factor.meta.overlap[factor.meta.overlap$labId %in% 
                      beat.meta.overlap[beat.meta.overlap$responseToInductionTx=="Refractory" &
                                          beat.meta.overlap$typeInductionTx=="Standard Chemotherapy",]$labId,]$CRWithInductionChemo <- FALSE


yn.factors <- c("cumulativeChemo","priorMalignancyRadiationTx","priorMDS",
                "priorMDSMoreThanTwoMths","priorMDSMPN","priorMalignancyNonMyeloid",
                "priorMDSMPNMoreThanTwoMths","priorMPN","priorMPNMoreThanTwoMths")
for (i in yn.factors) {
  factor.meta.overlap[factor.meta.overlap$labId %in% 
                        beat.meta.overlap[beat.meta.overlap[,i]=="y",]$labId,i] <- TRUE
  factor.meta.overlap[factor.meta.overlap$labId %in% 
                        beat.meta.overlap[beat.meta.overlap[,i]=="n",]$labId,i] <- FALSE
}

pn.factors <- c("FLT3-ITD","NPM1")
for (i in pn.factors) {
  factor.meta.overlap[factor.meta.overlap$labId %in% 
                        beat.meta.overlap[beat.meta.overlap[,i]=="positive",]$labId,i] <- TRUE
  factor.meta.overlap[factor.meta.overlap$labId %in% 
                        beat.meta.overlap[beat.meta.overlap[,i]=="negative",]$labId,i] <- FALSE
}

other.factors <- unique(beat.meta.overlap$specificDxAtAcquisition[duplicated(beat.meta.overlap$specificDxAtAcquisition)])
other.factors <- other.factors[other.factors!="Unknown"]
factor.meta.overlap[,paste0("specificDxAtAcquisition: ",other.factors)] <- NA
for (i in other.factors) {
  factor.meta.overlap[factor.meta.overlap$labId %in% 
                        beat.meta.overlap[beat.meta.overlap$specificDxAtAcquisition==i,]$labId,
                      paste0("specificDxAtAcquisition: ",i)] <- TRUE
  factor.meta.overlap[factor.meta.overlap$labId %in% 
                        beat.meta.overlap[beat.meta.overlap$specificDxAtAcquisition!=i &
                                            beat.meta.overlap$specificDxAtAcquisition!="Unknown",]$labId,
                      paste0("specificDxAtAcquisition: ",i)] <- FALSE
}

# later look at ELN2017 column
other.factors <- unique(beat.meta.overlap$ELN2017[duplicated(beat.meta.overlap$ELN2017)])
factor.meta.overlap[,paste0("ELN2017: ", other.factors)] <- NA
for (i in other.factors) {
  factor.meta.overlap[factor.meta.overlap$labId %in% 
                        beat.meta.overlap[beat.meta.overlap$ELN2017==i,]$labId,
                      paste0("ELN2017: ", i)] <- TRUE
  factor.meta.overlap[factor.meta.overlap$labId %in% 
                        beat.meta.overlap[beat.meta.overlap$ELN2017!=i,]$labId,
                      paste0("ELN2017: ", i)] <- FALSE
}

# later look at fabBlastMorphology column
other.factors <- na.omit(unique(beat.meta.overlap$fabBlastMorphology[duplicated(beat.meta.overlap$fabBlastMorphology)]))
factor.meta.overlap[,paste0("fabBlastMorphology: ", other.factors)] <- NA
for (i in other.factors) {
  factor.meta.overlap[factor.meta.overlap$labId %in% 
                        beat.meta.overlap[beat.meta.overlap$fabBlastMorphology==i,]$labId,
                      paste0("fabBlastMorphology: ",i)] <- TRUE
  factor.meta.overlap[factor.meta.overlap$labId %in% 
                        beat.meta.overlap[beat.meta.overlap$fabBlastMorphology!=i & 
                                            !is.na(beat.meta.overlap$fabBlastMorphology),]$labId,
                      paste0("fabBlastMorphology: ",i)] <- FALSE
}

# later look at diseaseStageAtSpecimenCollection column
other.factors <- unique(beat.meta.overlap$diseaseStageAtSpecimenCollection[duplicated(beat.meta.overlap$diseaseStageAtSpecimenCollection)])
other.factors <- other.factors[other.factors != "Unknown"]
factor.meta.overlap[,paste0("diseaseStageAtSpecimenCollection: ",other.factors)] <- NA
for (i in other.factors) {
  factor.meta.overlap[factor.meta.overlap$labId %in% 
                        beat.meta.overlap[beat.meta.overlap$diseaseStageAtSpecimenCollection==i,]$labId,
                      paste0("diseaseStageAtSpecimenCollection: ",i)] <- TRUE
  factor.meta.overlap[factor.meta.overlap$labId %in% 
                        beat.meta.overlap[beat.meta.overlap$diseaseStageAtSpecimenCollection!=i & 
                                            !is.na(beat.meta.overlap$diseaseStageAtSpecimenCollection) &
                                            beat.meta.overlap$diseaseStageAtSpecimenCollection!="Unknown",]$labId,
                      paste0("diseaseStageAtSpecimenCollection: ",i)] <- FALSE
}
reduced.meta <- meta.df[,c("Barcode.ID","InitialAMLDiagnosis","PostChemotherapy")] # FLT3-ITD is already included
colnames(reduced.meta)[1] <- "labId"
tf.meta.overlap <- merge(factor.meta.overlap, reduced.meta, by="labId")
write.csv(tf.meta.overlap,"Table_S1_notSorted_rnaProtVen_factors.csv", row.names=FALSE)

diffDx <- beat.meta.overlap[beat.meta.overlap$dxAtInclusion!=beat.meta.overlap$dxAtSpecimenAcquisition,]
# only sample where dxAtInclusion != dxAtSpecimenAcquisition: labId == "18-00149" dxAtInclusion == "MYELODYSPLASTIC SYNDROMES" and dxAtSpecimenAcquisition == "ACUTE MYELOID LEUKAEMIA (AML) AND RELATED PRECURSOR NEOPLASMS"

# t-tests
tf.test.input <- pt.r[,c(1,3)]
colnames(tf.test.input)[1] <- "labId"
tf.test.input <- merge(tf.test.input, tf.meta.overlap, by="labId")
write.csv(tf.test.input,"Table_S1_notSorted_rnaProtVen_factors_tTestInput.csv", row.names=FALSE)

metaFactor <- colnames(tf.test.input)[3:ncol(tf.test.input)]
tf.meta.test <- data.frame(metaFactor,meanTrue=NA,meanFalse=NA,medianTrue=NA,
                           medianFalse=NA,t=NA,p=NA,pTrueLess=NA,pTrueGreater=NA,
                           nTrue=NA,nFalse=NA)
for (i in metaFactor) {
  tf.meta.test[tf.meta.test$metaFactor==i,]$meanTrue <- mean(tf.test.input[tf.test.input[,i],2], na.rm = TRUE)
  tf.meta.test[tf.meta.test$metaFactor==i,]$meanFalse <- mean(tf.test.input[!tf.test.input[,i],2], na.rm = TRUE)
  
  tf.meta.test[tf.meta.test$metaFactor==i,]$medianTrue <- median(tf.test.input[tf.test.input[,i],2], na.rm = TRUE)
  tf.meta.test[tf.meta.test$metaFactor==i,]$medianFalse <- median(tf.test.input[!tf.test.input[,i],2], na.rm = TRUE)
  
  tf.meta.test[tf.meta.test$metaFactor==i,]$nTrue <- length(na.omit(tf.test.input[tf.test.input[,i],2]))
  tf.meta.test[tf.meta.test$metaFactor==i,]$nFalse <- length(na.omit(tf.test.input[!tf.test.input[,i],2]))
  
  tf.meta.test[tf.meta.test$metaFactor==i,]$p <- t.test(na.omit(tf.test.input[tf.test.input[,i],2]),
                      na.omit(tf.test.input[!tf.test.input[,i],2]))$p.value
  tf.meta.test[tf.meta.test$metaFactor==i,]$t <- t.test(na.omit(tf.test.input[tf.test.input[,i],2]),
                                                        na.omit(tf.test.input[!tf.test.input[,i],2]))$statistic
  tf.meta.test[tf.meta.test$metaFactor==i,]$pTrueLess <- t.test(na.omit(tf.test.input[tf.test.input[,i],2]),
                                                                na.omit(tf.test.input[!tf.test.input[,i],2]),
                      alternative = "less")$p.value
  tf.meta.test[tf.meta.test$metaFactor==i,]$pTrueGreater <- t.test(na.omit(tf.test.input[tf.test.input[,i],2]),
                                                                 na.omit(tf.test.input[!tf.test.input[,i],2]),
                      alternative = "greater")$p.value
} 
tf.meta.test$q <- qvalue::qvalue(tf.meta.test$p, pi0=1)$qvalues
write.csv(tf.meta.test,"Table_S1_notSorted_rnaProtVen_factors_tTest.csv", row.names=FALSE)
# M2, M4, AML with inv(16)(p13.1q22) or t(16;16)(p13.1;q22); CBFB-MYH11 are better predicted
# M1 is worse predicted

# bar plots: correlations, -logP; box plots for significantly differently predicted t-tests
dot.df <- sorted.meta.corr$result
dot.df$sig <- FALSE
if (any(dot.df$Pearson.q<0.05)) {
  dot.df[dot.df$Pearson.q<0.05,]$sig <- TRUE
}
cor.bar <- ggplot2::ggplot(dot.df, aes(x=Pearson.est, y=reorder(Metadata,Pearson.est), fill=sig)) +
  geom_bar(stat="identity") + scale_fill_manual(values=c("black","grey"), breaks=c(TRUE, FALSE), labels=c("q < 0.05","q > 0.05")) +
  theme_classic() + labs(x="Pearson r", fill="Significance") + theme(axis.title.y=element_blank())
cor.bar
ggsave("Table_S1_notSorted_rnaProtVen_numeric_corr_barPlot.pdf", cor.bar, width=5, height=5)
ggsave("Table_S1_notSorted_rnaProtVen_numeric_corr_barPlot_wh5.5.pdf", cor.bar, width=5.5, height=5.5)
ggsave("Table_S1_notSorted_rnaProtVen_numeric_corr_barPlot_w7.5h6.5.pdf", cor.bar, width=7.5, height=6.5)

dot.df <- tf.meta.test
dot.df$sig <- FALSE
if (any(dot.df$q<0.05)) {
  dot.df[dot.df$q<0.05,]$sig <- TRUE # only one is specificDxAtAcquisition: AML with inv(16)(p13.1q22) or t(16;16)(p13.1;q22); CBFB-MYH11 with 2.939324e-05
}
# test.bar <- ggplot2::ggplot(dot.df, aes(x=-log10(p), y=reorder(metaFactor,-p), fill=sig)) +
#   geom_bar(stat="identity") + scale_fill_manual(values=c("black","grey"), breaks=c(TRUE, FALSE)) +
#   theme_classic() + labs(x="-Log(P-value)") + theme(axis.title.y=element_blank())
test.bar <- ggplot2::ggplot(dot.df, aes(x=t, y=reorder(metaFactor,t), fill=sig)) +
  geom_bar(stat="identity") + scale_fill_manual(values=c("black","grey"), breaks=c(TRUE, FALSE), labels=c("q < 0.05","q > 0.05")) +
  theme_classic() + labs(x="t-score", fill="Significance") + theme(axis.title.y=element_blank())
test.bar
ggsave("Table_S1_notSorted_rnaProtVen_factor_tTest_barPlot.pdf", test.bar, width=7.5, height=6.5)

dot.df <- pt.r[,c(1,3)]
colnames(dot.df)[1] <- "labId"
dot.df <- merge(dot.df, beat.meta.overlap[,c("labId","specificDxAtAcquisition")], by="labId")
test.box <- ggplot2::ggplot(dot.df[dot.df$specificDxAtAcquisition!="Unknown",], 
                            aes(x=reorder(specificDxAtAcquisition,
                                          -`Sorted: 25 proteins`),
                                        y=`Sorted: 25 proteins`, 
                                        color=specificDxAtAcquisition)) +
  geom_violin(alpha=0) + geom_point() + geom_boxplot(width=0.2, alpha = 0) + 
  theme_classic() + labs(y="Pearson r", x="specificDxAtAcquisition") + 
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1), 
        legend.position="none", axis.title.x=element_blank())
test.box
ggsave("Table_S1_notSorted_rnaProtVen_factor_tTest_boxPlot_specificDxAtAcquisition_noUnknown.pdf", test.box, width=5, height=5)

#### use ven SSE ####
tf.test.input <- ven.corr[ven.corr$Signature=="Sorted: 25 proteins",c("Barcode.ID","delta_AUC_squared")]
colnames(tf.test.input)[1] <- "labId"
tf.test.input <- merge(tf.test.input, tf.meta.overlap, by="labId")
write.csv(tf.test.input,"Table_S1_notSorted_rnaProtVen_factors_tTestInput_venSSE.csv", row.names=FALSE)

metaFactor <- colnames(tf.test.input)[3:ncol(tf.test.input)]
tf.meta.test <- data.frame(metaFactor,meanTrue=NA,meanFalse=NA,medianTrue=NA,
                           medianFalse=NA,t=NA,p=NA,pTrueLess=NA,pTrueGreater=NA,
                           nTrue=NA,nFalse=NA)
for (i in metaFactor) {
  tf.meta.test[tf.meta.test$metaFactor==i,]$meanTrue <- mean(tf.test.input[tf.test.input[,i],2], na.rm = TRUE)
  tf.meta.test[tf.meta.test$metaFactor==i,]$meanFalse <- mean(tf.test.input[!tf.test.input[,i],2], na.rm = TRUE)
  
  tf.meta.test[tf.meta.test$metaFactor==i,]$medianTrue <- median(tf.test.input[tf.test.input[,i],2], na.rm = TRUE)
  tf.meta.test[tf.meta.test$metaFactor==i,]$medianFalse <- median(tf.test.input[!tf.test.input[,i],2], na.rm = TRUE)
  
  tf.meta.test[tf.meta.test$metaFactor==i,]$nTrue <- length(na.omit(tf.test.input[tf.test.input[,i],2]))
  tf.meta.test[tf.meta.test$metaFactor==i,]$nFalse <- length(na.omit(tf.test.input[!tf.test.input[,i],2]))
  
  tf.meta.test[tf.meta.test$metaFactor==i,]$p <- t.test(na.omit(tf.test.input[tf.test.input[,i],2]),
                                                        na.omit(tf.test.input[!tf.test.input[,i],2]))$p.value
  tf.meta.test[tf.meta.test$metaFactor==i,]$t <- t.test(na.omit(tf.test.input[tf.test.input[,i],2]),
                                                        na.omit(tf.test.input[!tf.test.input[,i],2]))$statistic
  tf.meta.test[tf.meta.test$metaFactor==i,]$pTrueLess <- t.test(na.omit(tf.test.input[tf.test.input[,i],2]),
                                                                na.omit(tf.test.input[!tf.test.input[,i],2]),
                                                                alternative = "less")$p.value
  tf.meta.test[tf.meta.test$metaFactor==i,]$pTrueGreater <- t.test(na.omit(tf.test.input[tf.test.input[,i],2]),
                                                                   na.omit(tf.test.input[!tf.test.input[,i],2]),
                                                                   alternative = "greater")$p.value
} 
tf.meta.test$q <- qvalue::qvalue(tf.meta.test$p, pi0=1)$qvalues
write.csv(tf.meta.test,"Table_S1_notSorted_rnaProtVen_factors_tTest_venSSE.csv", row.names=FALSE)

# bar plots: correlations, -logP; box plots for significantly differently predicted t-tests
dot.df <- sorted.meta.corr$result
dot.df$sig <- FALSE
if (any(dot.df$Pearson.q<0.05)) {
  dot.df[dot.df$Pearson.q<0.05,]$sig <- TRUE
}
cor.bar <- ggplot2::ggplot(dot.df, aes(x=Pearson.est, y=reorder(Metadata,Pearson.est), fill=sig)) +
  geom_bar(stat="identity") + scale_fill_manual(values=c("black","grey"), breaks=c(TRUE, FALSE), labels=c("q < 0.05","q > 0.05")) +
  theme_classic() + labs(x="Pearson r", fill="Significance") + theme(axis.title.y=element_blank())
cor.bar
ggsave("Table_S1_notSorted_rnaProtVen_numeric_corr_barPlot_venSSE.pdf", cor.bar, width=5, height=5)
ggsave("Table_S1_notSorted_rnaProtVen_numeric_corr_barPlot_wh5.5_venSSE.pdf", cor.bar, width=5.5, height=5.5)
ggsave("Table_S1_notSorted_rnaProtVen_numeric_corr_barPlot_w7.5h6.5_venSSE.pdf", cor.bar, width=7.5, height=6.5)

dot.df <- sorted.meta.corr$result
dot.df$sig <- FALSE
if (any(dot.df$Spearman.q<0.05)) {
  dot.df[dot.df$Spearman.q<0.05,]$sig <- TRUE
}
cor.bar <- ggplot2::ggplot(dot.df, aes(x=Spearman.est, y=reorder(Metadata,Spearman.est), fill=sig)) +
  geom_bar(stat="identity") + scale_fill_manual(values=c("black","grey"), breaks=c(TRUE, FALSE), labels=c("q < 0.05","q > 0.05")) +
  theme_classic() + labs(x="Spearman rho", fill="Significance") + theme(axis.title.y=element_blank())
cor.bar
ggsave("Table_S1_notSorted_rnaProtVen_numeric_corr_barPlot_venSSE_spearman.pdf", cor.bar, width=5, height=5)
ggsave("Table_S1_notSorted_rnaProtVen_numeric_corr_barPlot_wh5.5_venSSE_spearman.pdf", cor.bar, width=5.5, height=5.5)
ggsave("Table_S1_notSorted_rnaProtVen_numeric_corr_barPlot_w7.5h6.5_venSSE_spearman.pdf", cor.bar, width=7.5, height=6.5)


dot.df <- tf.meta.test
dot.df$sig <- FALSE
if (any(dot.df$q<0.05)) {
  dot.df[dot.df$q<0.05,]$sig <- TRUE # only one is specificDxAtAcquisition: AML with inv(16)(p13.1q22) or t(16;16)(p13.1;q22); CBFB-MYH11 with 2.939324e-05
}
# test.bar <- ggplot2::ggplot(dot.df, aes(x=-log10(p), y=reorder(metaFactor,-p), fill=sig)) +
#   geom_bar(stat="identity") + scale_fill_manual(values=c("black","grey"), breaks=c(TRUE, FALSE)) +
#   theme_classic() + labs(x="-Log(P-value)") + theme(axis.title.y=element_blank())
test.bar <- ggplot2::ggplot(dot.df, aes(x=t, y=reorder(metaFactor,t), fill=sig)) +
  geom_bar(stat="identity") + scale_fill_manual(values=c("black","grey"), breaks=c(TRUE, FALSE), labels=c("q < 0.05","q > 0.05")) +
  theme_classic() + labs(x="t-score", fill="Significance") + theme(axis.title.y=element_blank())
test.bar
ggsave("Table_S1_notSorted_rnaProtVen_factor_tTest_barPlot_venSSE.pdf", test.bar, width=7.5, height=6.5)

dot.df <- corr.input[,c("labId","delta_AUC_squared","allelic_ratio")]
Pearson.est <- cor.test(dot.df$delta_AUC_squared, dot.df$allelic_ratio)$estimate
Pearson.p <- cor.test(dot.df$delta_AUC_squared, dot.df$allelic_ratio)$p.value
stats_pearson <- substitute(
  r == est * "," ~ ~"p" ~ "=" ~ p,
  list(
    est = as.numeric(format(Pearson.est, digits = 3)),
    p = format(Pearson.p, digits = 3)
  )
)
Spearman.est <- cor.test(dot.df$delta_AUC_squared, dot.df$allelic_ratio, method="spearman")$estimate
Spearman.p <- cor.test(dot.df$delta_AUC_squared, dot.df$allelic_ratio, method="spearman")$p.value
stats_spearman <- substitute(
  rho == est * "," ~ ~"p" ~ "=" ~ p,
  list(
    est = as.numeric(format(Spearman.est, digits = 3)),
    p = format(Spearman.p, digits = 3)
  )
)
scatter.plot <- ggplot2::ggplot(data = dot.df,
                aes(x = allelic_ratio, y = delta_AUC_squared)) +
  ggplot2::geom_point() +
  ggplot2::labs(x = "Allelic Ratio", y = "Prediction Error") +
  ggplot2::geom_smooth(method = "lm", size = 1.5,
                       linetype = "dashed", color = "blue",
                       se = FALSE, na.rm = TRUE) +
  ggplot2::geom_text(
    x = 15, y = 20000, vjust = "inward", hjust = "inward",
    colour = "blue", parse = TRUE,
    label = as.character(as.expression(stats_pearson))) + 
  ggplot2::geom_text(
    x = 15, y = 18000, vjust = "inward", hjust = "inward",
    colour = "blue", parse = TRUE,
    label = as.character(as.expression(stats_spearman))) + theme_classic()
scatter.plot
ggsave("Table_S1_notSorted_rnaProtVen_factor_corr_scatterPlot_venSSE.pdf", scatter.plot, width=3, height=3)


npm1.meta <- npm1.meta[npm1.meta$labId %in% sorted.patients,] # 9 patients
write.csv(npm1.meta, "Table_S1_sortedDIApatients.csv", row.names=FALSE)

missing.pts <- sort(sorted.patients[!(sorted.patients %in% npm1.meta$labId)])
missing.pts # "16-01184" "19-00074" "19-00406" "20-00083" "20-00450" "21-00034" "21-00176" "21-00432" "21-00839" "22-00117" "22-00251" "22-00571" "22-00697" "23-00083"
found.pts <- sort(sorted.patients[(sorted.patients %in% npm1.meta$labId)])
found.pts # "17-01060" "18-00103" "18-00105" "18-00190" "18-00260" "18-00290" "18-00390" "19-00019" "19-00431"

#### check for genes of interest ####
# from papers Anupriya sent
# BCL2 and MDM inhibitor resistance in monocytic leukemia cells: https://ashpublications.org/blood/article/142/Supplement%201/2937/502904/Cebp-IL1-TNF-Positive-Feedback-Loop-Drives-Drug
# don't have access to full article
bcl2.mdm.res.up <- c("CEBPB", "MCL1","BCL2A1") # also need to add NFKB/IL1/TNF pathway
bcl2.mdm.res.dn <- c("CASP3","CASP6","BCL2","CDKN1A","PMAIP1","BBC3","BMF","TP53")

# https://pmc.ncbi.nlm.nih.gov/articles/PMC9131911/#sec2
ven.res <- c("CD11B","CD16","CD56","CD64","HLADR") # these are protein names - check for gene symbols
ven.res <- c("ITGAM","FCGR3A", "FCGR3B","NCAM1","FCGR1A","HLA-DRB1") # gene symbols
ven.res <- c("ITGAM","FCGR3A", "FCGR1A","HLA-DRB1") # FCGR3B and NCAM1 are not detected in our proteomics
ven.sens <- c("CD117") # protein name
ven.sens <- c("KIT") # gene symbol
dora.sens <- c("HLADR") # protein name
dora.sens <- c("HLA-DRB1") # gene symbol

# https://pmc.ncbi.nlm.nih.gov/articles/PMC10618724/#_ad93_
av.res <- c("CD14", "RAS") # which RAS?
av.sens <- c("CD117","IDH1","NPM1") # protein name
av.sens <- c("KIT", "IDH1","NPM1") # gene symbol
#ven.rux.sig <- readxl::read_excel("/Users/gara093/Downloads/bcd-23-0014_table_s7_suppst7.xlsx")

# https://aacrjournals.org/cancerdiscovery/article/13/6/1408/726964/Combinatorial-BCL2-Family-Expression-in-Acute
m5.dn <- c("CD117") # protein name
m5.dn <- c("KIT") # gene symbol
m5.up <- c("CD11b", "CD68","CD64") # protein name
m5.up <- c("ITGAM", "CD68", "FCGR1A") # gene symbol

# LSC signature different between sens/res AML
# https://aacrjournals.org/cancerdiscovery/article/13/6/1408/726964/Combinatorial-BCL2-Family-Expression-in-Acute
# https://pmc.ncbi.nlm.nih.gov/articles/PMC7124979/#_ad93_
lsc <- readxl::read_excel("/Users/gara093/Downloads/NIHMS1551453-supplement-2.xlsx", sheet="Table S4")
colnames(lsc) <- lsc[2,]
lsc <- lsc[-c(1,2),]
lsc <- as.list(lsc)

# get our sorted signature
sorted <- read.csv("~/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/analysis/combined24-27/using_cellType-sortType-patient_factors/DIA_2batches_noOutliers_noMSC_Bead/no_filter/CD14_Pos_vs_Neg/global/Differential_expression/Differential_expression_results.csv")
sorted <- sorted[sorted$adj.P.Val <= 0.05,]

anupriya.sigs <- append(lsc, list("BCL2 Res Down" = bcl2.mdm.res.dn, 
                                  "BCL2 Res Up" = bcl2.mdm.res.up, 
                                  "Ven Res" = ven.res, "Ven Sens" = ven.sens, 
                                  "Dora Sens" = dora.sens, "AzaVen Res" = av.res, 
                                  "AzaVen Sens" = av.sens, "M5 Down" = m5.dn,
                                  "M5 Up" = m5.up))
dir.create("Anupriya_sigs")
setwd("Anupriya_sigs")
for (i in names(anupriya.sigs)) {
  if (any(sorted$Gene %in% anupriya.sigs[[i]])) {
    ggvenn::ggvenn(list("Sorted" = sorted$Gene, i = anupriya.sigs[[i]]), 
                   show_elements=TRUE, show_percentage=FALSE)
    ggsave(paste0("venn_sorted_vs_",i,"_wElements.pdf"), width=5, height=5) 
    ggvenn::ggvenn(list("Sorted" = sorted$Gene, i = anupriya.sigs[[i]]), 
                   show_percentage=FALSE)
    ggsave(paste0("venn_sorted_vs_",i,".pdf"), width=5, height=5)
  }
  
  if (any(sorted[sorted$Log2FC>0,]$Gene %in% anupriya.sigs[[i]])) {
  ggvenn::ggvenn(list("Sorted" = sorted[sorted$Log2FC>0,]$Gene, i = anupriya.sigs[[i]]), 
                 show_elements=TRUE, show_percentage=FALSE)
  ggsave(paste0("venn_sorted_up_vs_",i,"_wElements.pdf"), width=5, height=5)
  ggvenn::ggvenn(list("Sorted" = sorted[sorted$Log2FC>0,]$Gene, i = anupriya.sigs[[i]]), 
                 show_percentage=FALSE)
  ggsave(paste0("venn_sorted_up_vs_",i,".pdf"), width=5, height=5)
  }
  
  if (any(sorted[sorted$Log2FC<0,] %in% anupriya.sigs[[i]])) {
    ggvenn::ggvenn(list("Sorted" = sorted[sorted$Log2FC<0,]$Gene, i = anupriya.sigs[[i]]), 
                   show_elements=TRUE, show_percentage=FALSE)
    ggsave(paste0("venn_sorted_dn_vs_",i,"_wElements.pdf"), width=5, height=5) 
    ggvenn::ggvenn(list("Sorted" = sorted[sorted$Log2FC<0,]$Gene, i = anupriya.sigs[[i]]), 
                  show_percentage=FALSE)
    ggsave(paste0("venn_sorted_dn_vs_",i,".pdf"), width=5, height=5) 
  }
}

#anu.sigs <- data.table::rbindlist(anupriya.sigs, use.names=TRUE, idcol="Signature")
#anu.sigs <- as.data.frame(anupriya.sigs)
anu.sigs <- rbind(anupriya.sigs)

anu.genes <- unique(c(unlist(anupriya.sigs), sorted$Gene))
anu.sigs <- data.frame(anu.genes, N_sigs = 0, Sigs = "")
anu.sigs[anu.sigs$anu.genes %in% sorted[sorted$Log2FC > 0,]$Gene,]$N_sigs <- 1
anu.sigs[anu.sigs$anu.genes %in% sorted[sorted$Log2FC > 0,]$Gene,]$Sigs <- "Sorted CD14+"

anu.sigs[anu.sigs$anu.genes %in% sorted[sorted$Log2FC < 0,]$Gene,]$N_sigs <- 1
anu.sigs[anu.sigs$anu.genes %in% sorted[sorted$Log2FC < 0,]$Gene,]$Sigs <- "Sorted CD34+"

for (i in names(anupriya.sigs)) {
  anu.sigs[anu.sigs$anu.genes %in% anupriya.sigs[[i]],]$N_sigs <- anu.sigs[anu.sigs$anu.genes %in% anupriya.sigs[[i]],]$N_sigs + 1
  anu.sigs[anu.sigs$anu.genes %in% anupriya.sigs[[i]],]$Sigs <- paste0(anu.sigs[anu.sigs$anu.genes %in% anupriya.sigs[[i]],]$Sigs, ", ", i)
}

anu.sigs[startsWith(anu.sigs$Sigs, ", "),]$Sigs <- sub(", ", "", anu.sigs[startsWith(anu.sigs$Sigs, ", "),]$Sigs)
write.csv(anu.sigs, "Overlap_across_signatures_from_Anupriya.csv", row.names=FALSE)

cd14.overlap.sigs <- unique(unlist(strsplit(unlist(anu.sigs[grepl("Sorted CD14+",anu.sigs$Sigs),]$Sigs),", "))) # 10 including CD14+
cd34.overlap.sigs <- unique(unlist(strsplit(unlist(anu.sigs[grepl("Sorted CD34+",anu.sigs$Sigs),]$Sigs),", "))) # 9 including CD34+
anupriya.sigs2 <- append(anupriya.sigs, list("Sorted CD14+" = sorted[sorted$Log2FC > 0,]$Gene,
                                             "Sorted CD34+" = sorted[sorted$Log2FC < 0,]$Gene))
saveRDS(anupriya.sigs2, "Signatures_from_Anupriya.rds")
# convert to matrix
myPlot = UpSetR::fromList(anupriya.sigs2)
pdf("sigsFromAnupriya_upsetPlot.pdf", onefile=FALSE)
UpSetR::upset(myPlot, order.by="freq")
dev.off()

write.csv(anu.sigs[anu.sigs$N_sigs > 1 & grepl("Sorted",anu.sigs$Sigs) &
                     (grepl("M5",anu.sigs$Sigs) | grepl("Res",anu.sigs$Sigs) | grepl("Sens",anu.sigs$Sigs)),], 
          "Overlap_across_signatures_from_Anupriya_sortedM5ResSens.csv", row.names=FALSE)

