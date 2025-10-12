# Differential expression & enrichment analyses: global & phospho
# PTRC2: exp 21
# Author: Belinda B. Garana
# Created: 2023-12-06
# Last edit: 2025-01-31

library(readxl); library(panSEA); library(synapser); library(plotly);
library(stringr); library(tidyr); library(ggplot2); library(MSnSet.utils)
#source("https://raw.githubusercontent.com/PNNL-CompBio/MPNST_Chr8/refs/heads/main/panSEA_helper_20240913.R")
#source("https://raw.githubusercontent.com/PNNL-CompBio/Exp21_NRAS_ASO_treated_patients/refs/heads/main/panSEA_helper_20240508_updated20250131.R")
source("~/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Chr8/panSEA_helper_20240913.R")
source("~/OneDrive - PNNL/Documents/GitHub/Exp21_NRAS_ASO_treated_patients/proteomics/panSEA_helper_20240508_updated20250131.R")
source("~/OneDrive - PNNL/Documents/GitHub/Exp21_NRAS_ASO_treated_patients/proteomics/customPCA.R")
source("~/OneDrive - PNNL/Documents/GitHub/Exp21_NRAS_ASO_treated_patients/proteomics/PCA3D.R")
synapser::synLogin()

rasPlots <- function(temp.expr, cc.df, meta.df.nras) {
  # violin plots of ras genes
  long.temp.expr <- reshape2::melt(temp.expr)
  long.temp.expr$variable <- as.character(long.temp.expr$variable)
  long.temp.expr$Tube <- substr(long.temp.expr$variable, 2, nchar(long.temp.expr$variable))
  long.temp.expr <- merge(long.temp.expr, meta.df.nras, by="Tube")
  
  Gene <- unique(temp.expr$Gene)
  p.vals <- data.frame(Gene, Lower = NA, Higher = NA, Two.sided = NA)
  for (temp.ras in Gene) {
    p.df <- long.temp.expr[long.temp.expr$Gene == temp.ras,]
    p.df$Treatment <- factor(p.df$Treatment, levels=unique(p.df$Treatment))
    temp.plot <- ggplot2::ggplot(p.df, aes(x=Treatment, y=value)) + geom_violin(alpha=0) +
      geom_point(aes(color=PatientID)) + geom_boxplot(width=0.2, alpha = 0) +
      labs(y=paste("Relative Abundance of",temp.ras)) + theme_classic(base_size = 12) +
      theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),
            axis.title.x = element_blank())
    
    # t-tests
    temp.nras <- as.numeric(p.df[grepl("NRAS", p.df$Treatment, ignore.case = TRUE),]$value, na.rm = TRUE)
    temp.other <- as.numeric(p.df[!grepl("NRAS", p.df$Treatment, ignore.case = TRUE),]$value, na.rm = TRUE)
    temp.higher <- t.test(temp.nras, temp.other, alternative="greater")$p.value
    p.vals[p.vals$Gene == temp.ras,]$Higher <- temp.higher
    temp.lower <- t.test(temp.nras, temp.other, alternative="less")$p.value
    p.vals[p.vals$Gene == temp.ras,]$Lower <- temp.lower
    temp2 <- t.test(temp.nras, temp.other, alternative="two.sided")$p.value
    p.vals[p.vals$Gene == temp.ras,]$Two.sided <- temp2
    if (temp.higher <= 0.05) {
      temp.title <- paste0(temp.ras," expression is higher with NRAS ASO (p = ", signif(temp.higher,2),")")
      temp.plot <- temp.plot + ggtitle(temp.title)
    } else if (temp.lower <= 0.05) {
      temp.title <- paste0(temp.ras," expression is lower with NRAS ASO (p = ", signif(temp.lower,2),")")
      temp.plot <- temp.plot + ggtitle(temp.title)
    }
    ggsave(paste0(temp.ras,"_expression.pdf"), temp.plot, width = 8, height = 5)
  }
  write.csv(p.vals, "RAS_expression_pVals.csv", row.names = FALSE)
  
  if (any(p.vals$Lower <= 0.05) | any(p.vals$Higher <= 0.05) | any(p.vals$Two.sided <= 0.05)) {
    long.p.vals <- reshape2::melt(p.vals)
    colnames(long.p.vals)[2:3] <- c("Direction", "p")
    long.p.vals$minusLogP <- -log10(long.p.vals$p)
    long.p.vals <- long.p.vals[long.p.vals$p <= 0.05,]
    geneOrder <- long.p.vals[order(long.p.vals$p),]$Gene
    ggplot2::ggplot(long.p.vals, aes(x=Gene, y=minusLogP, fill=Direction)) +
      geom_bar(stat="identity", position="dodge") + theme_classic(base_size=12) + labs(y="-Log(P-value)") +
      scale_x_discrete(limits=geneOrder)+
      theme(axis.title.x = element_blank(), axis.text.x=element_text(angle=45,vjust=1,hjust=1))
    ggsave("Altered_RAS_expression_minusLogP.pdf", width=ifelse(0.35*nrow(long.p.vals)<3,3,0.35*nrow(long.p.vals)), height=3)
    
    dn.genes <- long.p.vals[long.p.vals$Direction=="Lower",]$Gene
    up.genes <- long.p.vals[long.p.vals$Direction=="Higher",]$Gene
    long.p.vals2 <- long.p.vals[long.p.vals$Direction=="Two.sided",]
    if (any(long.p.vals2$Gene %in% dn.genes)) {
      long.p.vals2[long.p.vals2$Gene %in% dn.genes,]$Direction <- "Lower"
    }
    if (any(long.p.vals2$Gene %in% up.genes)) {
      long.p.vals2[long.p.vals2$Gene %in% up.genes,]$Direction <- "Higher"
    }
    geneOrder <- long.p.vals2[order(long.p.vals2$p),]$Gene
    ggplot2::ggplot(long.p.vals2, aes(x=Gene, y=minusLogP, fill=Direction)) +
      geom_bar(stat="identity", position=position_dodge2(preserve="single")
      ) + theme_classic(base_size=12) + labs(y="-Log(P-value)") +
      scale_x_discrete(limits=geneOrder)+
      theme(axis.title.x = element_blank(), axis.text.x=element_text(angle=45,vjust=1,hjust=1))
    ggsave("Altered_RAS_expression_minusLogP_2sided.pdf", width=ifelse(0.35*nrow(long.p.vals)<3,3,0.35*nrow(long.p.vals)), height=3)
    
    long.p.vals2 <- long.p.vals[long.p.vals$Direction!="Two.sided",]
    geneOrder <- long.p.vals2[order(long.p.vals2$p),]$Gene
    ggplot2::ggplot(long.p.vals2, aes(x=Gene, y=minusLogP, fill=Direction)) +
      geom_bar(stat="identity", position=position_dodge2(preserve="single")
      ) + theme_classic(base_size=12) + labs(y="-Log(P-value)") +
      scale_x_discrete(limits=geneOrder)+
      theme(axis.title.x = element_blank(), axis.text.x=element_text(angle=45,vjust=1,hjust=1))
    ggsave("Altered_RAS_expression_minusLogP_1sided.pdf", width=ifelse(0.35*nrow(long.p.vals)<3,3,0.35*nrow(long.p.vals)), height=3)
  }
  
  # heatmap of ras genes
  expr.mat <- dplyr::select_if(temp.expr, is.numeric)
  expr.mat <- filter_for_hclust(expr.mat)
  
  if (nrow(expr.mat) > 1) {
    # create heatmaps
    expr.mat <- as.matrix(expr.mat)
    if (TRUE) {
      temp.heatmap <- try(R.utils::withTimeout(pheatmap::pheatmap(expr.mat, color =
                                                                    colorRampPalette(
                                                                      c("navy", "white", "firebrick3"))(50),
                                                                  cluster_rows = TRUE, cluster_cols = TRUE,
                                                                  scale = "row", annotation_col = cc.df,
                                                                  angle_col = "45",
                                                                  show_colnames = FALSE,
                                                                  fontsize = 10), timeout = 300, onTimeout="error"), silent = TRUE)
      if (!inherits(temp.heatmap, "try-error")) {
        my.clust.heatmaps <- temp.heatmap
      }
    }
    temp.heatmap <- try(R.utils::withTimeout(pheatmap::pheatmap(expr.mat, color =
                                                                  colorRampPalette(
                                                                    c("navy", "white", "firebrick3"))(50),
                                                                cluster_rows = TRUE, cluster_cols = TRUE,
                                                                annotation_col = cc.df,
                                                                angle_col = "45",
                                                                show_colnames = FALSE,
                                                                fontsize = 10), timeout = 300, onTimeout="error"), silent = TRUE)
    if (!inherits(temp.heatmap, "try-error")) {
      my.abs.heatmaps <- temp.heatmap
    }
  }
  all.heatmaps <- list("Not_scaled.bp" = my.abs.heatmaps,
                       "Scaled.bp" = my.clust.heatmaps)
  poi.files <- list("Expression_of_RAS_genes.csv" =
                      temp.expr,
                    "Heatmaps" = all.heatmaps)
  save_to_synapse(poi.files)
}
#### 1. Import BeatAML data for DMEA ####
setwd("~/OneDrive - PNNL/Documents/GitHub/Exp21_NRAS_ASO_treated_patients/proteomics/data")
# load data from Synapse and store locally
BeatAML <- load_BeatAML_for_DMEA()

# check patients included in exp 21
meta.df <- readxl::read_excel(synapser::synGet("syn51421353")$path, 
                              sheet = 2) # 80 samples

# remove samples that Sunil said to exclude
excluded.patients <- c("20-00576", "20-00432", "18-00101", "17-01083", "16-00815", "16-01253")
meta.df <- meta.df[!(meta.df$PatientID %in% excluded.patients),] # 50 samples
included.patients <- unique(meta.df$PatientID)

# remove patients included in this study
global.BeatAML <- BeatAML$global[!(BeatAML$global$Barcode.ID %in% included.patients),] # 210 patients, 7084 proteins
phospho.BeatAML <- BeatAML$phospho[!(BeatAML$phospho$Barcode.ID %in% included.patients),] # 210 patients, 4077 phosphoproteins
drug.BeatAML <- BeatAML$drug[!(BeatAML$drug$Barcode.ID %in% included.patients),] # 210 patients, 321 drugs
temp.expr <- list("Global" = global.BeatAML,
                  "Phospho" = phospho.BeatAML)
# temp.expr <- list("Phospho" = phospho.BeatAML)

#### 2. Run panSEA across experiments, omics for each contrast (NRAS, Gilt) ####
setwd("~/OneDrive - PNNL/Documents/GitHub/Exp21_NRAS_ASO_treated_patients/proteomics/")
dir.create("analysis")
setwd("analysis")
base.path <- "~/OneDrive - PNNL/Documents/GitHub/Exp21_NRAS_ASO_treated_patients/proteomics/analysis/"

pt <- list("Meta" = "syn64732333", # exp 21
           "Global" = "syn52751324",
           "Phospho" = "syn52751606")
cell <- list("Meta" = "syn52653527", # exp 20
           "Global" = "syn52623904",
           "Phospho" = "syn52623908")
cellKO <- list("Meta" = "syn64731726", # exp 22
             "Global" = "syn53212316",
             "Phospho" = "syn53212321")
inputs <- list("Patient ASO" = pt,
               "Cell Line ASO" = cell,
               "Cell Line KO" = cellKO)
meta.sheet <- c(2,1,2)
names(meta.sheet) <- names(inputs)
synIDs <- c("syn53180709","syn52256487","syn53212344")
names(synIDs) <- names(inputs)
contrast2 <- list("Patient ASO" = "PatientID",
                  "Cell Line ASO" = c("Time", "AML_L"),
                  "Cell Line KO" = c("Volume"))

for (i in names(inputs)) {
  setwd(base.path)
  dir.create(i)
  temp.path <- file.path(base.path,i)
  setwd(temp.path)
  ## read in data
  meta.df <- readxl::read_excel(synapser::synGet(inputs[[i]]$Meta)$path, 
                                sheet = meta.sheet[[i]])
  global.df <- read.table(synapser::synGet(inputs[[i]]$Global)$path, sep = "\t")
  phospho.df <- read.table(synapser::synGet(inputs[[i]]$Phospho)$path, sep = "\t")
  
  # # check ksdb coverage
  # gmt2 <- get_gmt2(gmt.list2 = "ksdb_human", force = TRUE)
  # coverage.df <- get_coverage(phospho.df, gmt2[[1]], "SUB_SITE")
  # write.csv(coverage.df, paste0(i, "_ksdb_coverage.csv"), row.names = FALSE)
  # hist(coverage.df$N_features_covered[coverage.df$N_features_covered<=10 & coverage.df$N_features_covered>0])
  
  # remove samples that Sunil said to exclude
  if (grepl("patient", i, ignore.case = TRUE)) {
    excluded.patients <- c("20-00576", "20-00432", "18-00101", "17-01083", "16-00815", "16-01253") # these are Ven resistant
    meta.df <- meta.df[!(meta.df$PatientID %in% excluded.patients),] 
  }
  
  # make sure contrasts are valid
  for (j in contrast2[[i]]) {
    for (k in 1:nrow(meta.df)) {
      meta.df[k,j] <- make.names(meta.df[k,j])
    }
  }
  
  # add column for feature names and later make it the first column
  global.df$Gene <- rownames(global.df)
  phospho.df$SUB_SITE <- rownames(phospho.df)
  
  ## prep contrasts
  temp.meta <- list()
  if (grepl("patient", i, ignore.case = TRUE)) { # exp 21
    contrasts <- c("NRAS", "NRAS_w_Gilt_ctl", "Gilt", "Gilt_w_NRAS_ctl", "Gilt_w_NRAS")
    # NRAS ASO vs. NEG ASO
    aso.treatments <- c("120HR 5uM NRAS ASO", "120HR 5uM NEG ASO")
    meta.df.nras <- as.data.frame(meta.df[meta.df$Treatment %in% aso.treatments,])
    meta.df.nras$NRAS <- "Untreated"
    meta.df.nras[meta.df.nras$Treatment == "120HR 5uM NRAS ASO",]$NRAS <- "Treated"
    rownames(meta.df.nras) <- paste0("X",meta.df.nras$Tube)
    temp.meta[[contrasts[1]]] <- meta.df.nras

    # Gilt + NRAS vs. Gilt
    nras.gilt.treatments <- c("120HR 5uM NRAS ASO + 6HR 100nM GILT", "120HR NEG ASO + 6HR 100nM GILT")
    meta.df.nras.gilt <- as.data.frame(meta.df[meta.df$Treatment %in% nras.gilt.treatments,])
    meta.df.nras.gilt$NRAS_w_Gilt_ctl <- "Untreated"
    meta.df.nras.gilt[meta.df.nras.gilt$Treatment == "120HR 5uM NRAS ASO + 6HR 100nM GILT",]$NRAS_w_Gilt_ctl <- "Treated"
    rownames(meta.df.nras.gilt) <- paste0("X",meta.df.nras.gilt$Tube)
    temp.meta[[contrasts[2]]] <- meta.df.nras.gilt

    # Gilt + NEG ASO vs. NEG ASO
    gilt.treatments <- c("120HR 5uM NEG ASO", "120HR NEG ASO + 6HR 100nM GILT")
    meta.df.gilt <- as.data.frame(meta.df[meta.df$Treatment %in% gilt.treatments,])
    meta.df.gilt$Gilt <- "Untreated"
    meta.df.gilt[meta.df.gilt$Treatment == "120HR NEG ASO + 6HR 100nM GILT",]$Gilt <- "Treated"
    rownames(meta.df.gilt) <- paste0("X",meta.df.gilt$Tube)
    temp.meta[[contrasts[3]]] <- meta.df.gilt
    
    # Gilt + NRAS vs. NRAS
    nras.gilt.treatments <- c("120HR 5uM NRAS ASO + 6HR 100nM GILT", "120HR 5uM NRAS ASO")
    meta.df.nras.gilt <- as.data.frame(meta.df[meta.df$Treatment %in% nras.gilt.treatments,])
    meta.df.nras.gilt$Gilt_w_NRAS_ctl <- "Untreated"
    meta.df.nras.gilt[meta.df.nras.gilt$Treatment == "120HR 5uM NRAS ASO + 6HR 100nM GILT",]$Gilt_w_NRAS_ctl <- "Treated"
    rownames(meta.df.nras.gilt) <- paste0("X",meta.df.nras.gilt$Tube)
    temp.meta[[contrasts[4]]] <- meta.df.nras.gilt
    
    # Gilt + NRAS
    nras.gilt.treatments <- c("120HR 5uM NRAS ASO + 6HR 100nM GILT", "120HR 5uM NEG ASO")
    meta.df.nras.gilt <- as.data.frame(meta.df[meta.df$Treatment %in% nras.gilt.treatments,])
    meta.df.nras.gilt$Gilt_w_NRAS <- "Untreated"
    meta.df.nras.gilt[meta.df.nras.gilt$Treatment == "120HR 5uM NRAS ASO + 6HR 100nM GILT",]$Gilt_w_NRAS <- "Treated"
    rownames(meta.df.nras.gilt) <- paste0("X",meta.df.nras.gilt$Tube)
    temp.meta[[contrasts[5]]] <- meta.df.nras.gilt
    
    # PCA
    meta.df$Tube <- as.numeric(meta.df$Tube)
    meta.df <- meta.df[!is.na(meta.df$Tube) & meta.df$Include,]
    meta.df[grepl("BASELINE",meta.df$Treatment),]$Treatment <- "BASELINE (No ASO)"
    meta.df <- as.data.frame(meta.df)
    rownames(meta.df) <- paste0("X", meta.df$Tube)
    temp.omics <- list("Global" = global.df,
                       "Phospho" = phospho.df)
    for (k in names(temp.omics)) {
      m <- MSnSet(exprs = temp.omics[[k]][,rownames(meta.df)] %>% as.matrix(),
                            pData = meta.df)
      MSnSet.utils::plot_pca(m, phenotype = "Treatment", label = "PatientID") #+ ggtitle(paste0(k," PCA"))
      ggsave(paste0(i, "_", k,"_corrected_PCA_by_", "Treatment", "_wPatientIDLabels.pdf"), width=7, height=7)
      MSnSet.utils::plot_pca(m, phenotype = "Treatment") #+ ggtitle(paste0(k," PCA"))
      ggsave(paste0(i, "_", k,"_corrected_PCA_by_", "Treatment", ".pdf"), width=7, height=7)
      MSnSet.utils::plot_pca(m, phenotype = "Treatment", shape=as.factor(meta.df$PatientID)) #+ ggtitle(paste0(k," PCA"))
      ggsave(paste0(i, "_", k,"_corrected_PCA_by_", "Treatment", "_wPatientIDShape.pdf"), width=7, height=7)
      
      m_corrected <- correct_batch_effect_NA(m, "PatientID", par.prior = T)
      MSnSet.utils::plot_pca(m_corrected, phenotype = "Treatment", label = "PatientID") #+ ggtitle(paste0(k," PCA"))
      ggsave(paste0(i, "_", k,"_patientCorrected_PCA_by_", "Treatment", "_wPatientIDLabels.pdf"), width=7, height=7)
      MSnSet.utils::plot_pca(m_corrected, phenotype = "Treatment") #+ ggtitle(paste0(k," PCA"))
      ggsave(paste0(i, "_", k,"_patientCorrected_PCA_by_", "Treatment", ".pdf"), width=7, height=7)
      MSnSet.utils::plot_pca(m_corrected, phenotype = "Treatment", shape=as.factor(meta.df$PatientID)) #+ ggtitle(paste0(k," PCA"))
      ggsave(paste0(i, "_", k,"_patientCorrected_PCA_by_", "Treatment", "_wPatientIDShape.pdf"), width=7, height=7)
    }

    expr.df <- global.df[,c("Gene", rownames(meta.df.nras))]
    cc.df <- meta.df.nras[,c("Treatment", contrast2[[i]])]
    temp.expr <- expr.df[grepl("RAS", expr.df$Gene, ignore.case=TRUE), ]
    rownames(temp.expr) <- temp.expr$Gene

    rasPlots(temp.expr, cc.df, meta.df.nras)
  } else if (grepl("ASO", i, ignore.case = TRUE)) { # exp 20
    contrasts <- c("NRAS", "NRAS_w_Gilt_ctl", "Gilt", "Gilt_w_NRAS_ctl", "Gilt_w_NRAS")
    # NRAS ASO vs. NEG ASO
    aso.treatments <- c("NRAS ASO", "CTRL ASO")
    meta.df.nras <- as.data.frame(meta.df[meta.df$Group %in% aso.treatments,])
    meta.df.nras$NRAS <- "Untreated"
    meta.df.nras[grepl("NRAS ASO", meta.df.nras$Sample),]$NRAS <- "Treated"
    rownames(meta.df.nras) <- meta.df.nras$MeasurementName
    temp.meta[[contrasts[1]]] <- meta.df.nras

    # NRAS w Gilt ctl
    nras.gilt.treatments <- c("NRAS ASO + GILT", "CTRL ASO + GILT")
    meta.df.nras.gilt <- as.data.frame(meta.df[meta.df$Group %in% nras.gilt.treatments,])
    meta.df.nras.gilt$NRAS_w_Gilt_ctl <- "Untreated"
    meta.df.nras.gilt[meta.df.nras.gilt$Group == "NRAS ASO + GILT",]$NRAS_w_Gilt_ctl <- "Treated"
    rownames(meta.df.nras.gilt) <- meta.df.nras.gilt$MeasurementName
    temp.meta[[contrasts[2]]] <- meta.df.nras.gilt

    # Gilt
    gilt.treatments <- c("CTRL ASO + GILT", "CTRL ASO")
    meta.df.gilt <- as.data.frame(meta.df[meta.df$Group %in% gilt.treatments,])
    meta.df.gilt$Gilt <- "Untreated"
    meta.df.gilt[meta.df.gilt$Group == "CTRL ASO + GILT",]$Gilt <- "Treated"
    rownames(meta.df.gilt) <- meta.df.gilt$MeasurementName
    temp.meta[[contrasts[3]]] <- meta.df.gilt
    
    # Gilt w NRAS ctl
    nras.gilt.treatments <- c("NRAS ASO + GILT", "NRAS ASO")
    meta.df.nras.gilt <- as.data.frame(meta.df[meta.df$Group %in% nras.gilt.treatments,])
    meta.df.nras.gilt$Gilt_w_NRAS_ctl <- "Untreated"
    meta.df.nras.gilt[meta.df.nras.gilt$Group == "NRAS ASO + GILT",]$Gilt_w_NRAS_ctl <- "Treated"
    rownames(meta.df.nras.gilt) <- meta.df.nras.gilt$MeasurementName
    temp.meta[[contrasts[4]]] <- meta.df.nras.gilt
    
    # Gilt + NRAS ASO
    nras.gilt.treatments <- c("NRAS ASO + GILT", "CTRL ASO")
    meta.df.nras.gilt <- as.data.frame(meta.df[meta.df$Group %in% nras.gilt.treatments,])
    meta.df.nras.gilt$Gilt_w_NRAS <- "Untreated"
    meta.df.nras.gilt[meta.df.nras.gilt$Group == "NRAS ASO + GILT",]$Gilt_w_NRAS <- "Treated"
    rownames(meta.df.nras.gilt) <- meta.df.nras.gilt$MeasurementName
    temp.meta[[contrasts[5]]] <- meta.df.nras.gilt
    
    meta.df[grepl("baseline",meta.df$Group),]$Group <- "BASELINE (No ASO)"
    meta.df[grepl("baseline",meta.df$Time),]$Time <- "0h"
    meta.df[grepl("2.HR",meta.df$Time),]$Time <- "2h"
    meta.df[grepl("6.HR",meta.df$Time),]$Time <- "6h"
    meta.df[grepl("24.HR",meta.df$Time),]$Time <- "24h"
    meta.df$Time <- factor(meta.df$Time, levels=c("0h", "2h", "6h", "24h"))
    meta.df$AML_L <- substr(meta.df$AML_L, 2, nchar(meta.df$AML_L))
    meta.df$TimeAndL <- paste0(meta.df$Time, ", ", meta.df$AML_L)
    meta.df$GroupAndTime <- paste0(meta.df$Group, ", ", meta.df$Time)
    meta.df <- as.data.frame(meta.df)
    rownames(meta.df) <- meta.df$MeasurementName
    temp.omics <- list("Global" = global.df,
                       "Phospho" = phospho.df)
    temp.omics.corr <- list()
    for (k in names(temp.omics)) { # 6854 complete rows for global, 25025 complete phospho rows
      m <- MSnSet(exprs = temp.omics[[k]][,rownames(meta.df)] %>% as.matrix(),
                            pData = meta.df)
      MSnSet.utils::plot_pca(m, phenotype = "Group", label = "AML_L") #+ ggtitle(paste0(k," PCA"))
      ggsave(paste0(i, "_", k,"_corrected_PCA_by_", "Group", "_wLLabels.pdf"), width=10, height=5)
      MSnSet.utils::plot_pca(m, phenotype = "Group", label = "Time") #+ ggtitle(paste0(k," PCA"))
      ggsave(paste0(i, "_", k,"_corrected_PCA_by_", "Group", "_wTimeLabels.pdf"), width=5, height=5)
      MSnSet.utils::plot_pca(m, phenotype = "Group", shape = as.factor(meta.df$Time)) #+ ggtitle(paste0(k," PCA"))
      ggsave(paste0(i, "_", k,"_corrected_PCA_by_", "Group", "_wTimeShape.pdf"), width=5, height=5)
      MSnSet.utils::plot_pca(m, phenotype = "Group", label = "TimeAndL") #+ ggtitle(paste0(k," PCA"))
      ggsave(paste0(i, "_", k,"_corrected_PCA_by_", "Group", "_wTimeAndLLabels.pdf"), width=5, height=5) # was w10 h10
      MSnSet.utils::plot_pca(m, phenotype = "GroupAndTime", shape = as.factor(meta.df$Time)) #+ ggtitle(paste0(k," PCA"))
      ggsave(paste0(i, "_", k,"_corrected_PCA_by_", "GroupAndTime", "_wTimeShape.pdf"), width=5, height=5) # was w10 h10
      customPCA(m, phenotype = "Group", phenotype2="Time", shape = as.factor(meta.df$Time)) #+ ggtitle(paste0(k," PCA"))
      ggsave(paste0(i, "_", k,"_corrected_PCA_by_GroupAndTime.pdf"), width=5, height=4) # was w10 h10
      PCA3D(m, phenotype = "Group", phenotype2="Time") #+ ggtitle(paste0(k," PCA"))
      ggsave(paste0(i, "_", k,"_corrected_3DPCA_by_GroupAndTime.pdf"), width=5, height=4) # was w10 h10
      
      # m_corrected <- correct_batch_effect_NA(m, "TimeAndL", par.prior = T)
      # MSnSet.utils::plot_pca(m_corrected, phenotype = "Group", label = "AML_L") #+ ggtitle(paste0(k," PCA"))
      # ggsave(paste0(i, "_", k,"_volAndTimeCorrected_PCA_by_", "Group", "_wLLabels.pdf"), width=10, height=5)
      # MSnSet.utils::plot_pca(m_corrected, phenotype = "Group", label = "Time") #+ ggtitle(paste0(k," PCA"))
      # ggsave(paste0(i, "_", k,"_volAndTimeCorrected_PCA_by_", "Group", "_wTimeLabels.pdf"), width=5, height=5)
      # MSnSet.utils::plot_pca(m_corrected, phenotype = "Group", shape = as.factor(meta.df$Time)) #+ ggtitle(paste0(k," PCA"))
      # ggsave(paste0(i, "_", k,"_volAndTimeCorrected_PCA_by_", "Group", "_wTimeShape.pdf"), width=5, height=5)
      # MSnSet.utils::plot_pca(m_corrected, phenotype = "Group", label = "TimeAndL") #+ ggtitle(paste0(k," PCA"))
      # ggsave(paste0(i, "_", k,"_volAndTimeCorrected_PCA_by_", "Group", "_wTimeAndLLabels.pdf"), width=5, height=5) # was w10 h10
      # MSnSet.utils::plot_pca(m_corrected, phenotype = "GroupAndTime", shape = as.factor(meta.df$Time)) #+ ggtitle(paste0(k," PCA"))
      # ggsave(paste0(i, "_", k,"_volAndTimeCorrected_PCA_by_", "GroupAndTime", "_wTimeShape.pdf"), width=5, height=5) # was w10 h10
      # customPCA(m_corrected, phenotype = "Group", phenotype2="Time", shape = as.factor(meta.df$Time)) #+ ggtitle(paste0(k," PCA"))
      # ggsave(paste0(i, "_", k,"_volAndTimeCorrected_PCA_by_GroupAndTime.pdf"), width=5, height=4) # was w10 h10
      # PCA3D(m_corrected, phenotype = "Group", phenotype2="Time") #+ ggtitle(paste0(k," PCA"))
      # ggsave(paste0(i, "_", k,"_volAndTimeCorrected_3DPCA_by_GroupAndTime.pdf"), width=5, height=4) # was w10 h10
      
      m_corrected <- correct_batch_effect_NA(m, "AML_L", par.prior = T)
      temp.omics.corr[[k]] <- exprs(m_corrected)
      write.table(exprs(m_corrected),
                  file = paste0(i, "_", k,"_ligandCorrected.txt"),
                  quote=F, sep="\t")
      
      MSnSet.utils::plot_pca(m_corrected, phenotype = "Group", label = "AML_L") #+ ggtitle(paste0(k," PCA"))
      ggsave(paste0(i, "_", k,"_ligandCorrected_PCA_by_", "Group", "_wLLabels.pdf"), width=10, height=5)
      MSnSet.utils::plot_pca(m_corrected, phenotype = "Group", label = "Time") #+ ggtitle(paste0(k," PCA"))
      ggsave(paste0(i, "_", k,"_ligandCorrected_PCA_by_", "Group", "_wTimeLabels.pdf"), width=5, height=5)
      MSnSet.utils::plot_pca(m_corrected, phenotype = "Group", shape = meta.df$Time)
      ggsave(paste0(i, "_", k,"_ligandCorrected_PCA_by_", "Group", "_wTimeShape.pdf"), width=5, height=5)
      MSnSet.utils::plot_pca(m_corrected, phenotype = "Group",alpha=0) +
        geom_point(aes(shape=meta.df$Time), alpha=1) + labs(shape="Time (hr)") + 
        scale_shape_manual(values=c(16, 3, 4, 17), breaks=c("0h","2h","6h","24h"), labels=c(0, 2, 6, 24)) + 
        scale_color_manual(values=c(rgb(246,203, 82, maxColorValue=255), rgb(7, 6, 5, maxColorValue=255), 
                                    rgb(97, 181,125, maxColorValue=255), rgb(57, 112, 176, maxColorValue=255), 
                                    rgb(168, 80, 152, maxColorValue=255)), breaks=levels(m_corrected$Group)) + 
        scale_fill_manual(values=c(rgb(246,203, 82, maxColorValue=255), rgb(7, 6, 5, maxColorValue=255), 
                                    rgb(97, 181,125, maxColorValue=255), rgb(57, 112, 176, maxColorValue=255), 
                                    rgb(168, 80, 152, maxColorValue=255)), breaks=levels(m_corrected$Group)) +
        ggtitle(ifelse(k=="Global","Proteomics","Phosphoproteomics")) + theme(plot.title=element_text(face="bold"), legend.key=element_rect(fill="transparent")) # legend not fixed?
      ggsave(paste0(i, "_", k,"_ligandCorrected_PCA_by_", "Group", "_wTimeShape_v2.pdf"), width=5, height=5)
      MSnSet.utils::plot_pca(m_corrected, phenotype = "Group", label = "TimeAndL") #+ ggtitle(paste0(k," PCA"))
      ggsave(paste0(i, "_", k,"_ligandCorrected_PCA_by_", "Group", "_wTimeAndLLabels.pdf"), width=5, height=5) # was w10 h10
      MSnSet.utils::plot_pca(m_corrected, phenotype = "GroupAndTime", shape = as.factor(meta.df$Time))
      ggsave(paste0(i, "_", k,"_ligandCorrected_PCA_by_", "GroupAndTime", "_wTimeShape.pdf"), width=5, height=5)# was w10 h10
      customPCA(m_corrected, phenotype = "Group", phenotype2="Time", shape = as.factor(meta.df$Time)) #+ ggtitle(paste0(k," PCA"))
      ggsave(paste0(i, "_", k,"_ligandCorrected_PCA_by_GroupAndTime.pdf"), width=5, height=4) # was w10 h10
      PCA3D(m_corrected, phenotype = "Group", phenotype2="Time") #+ ggtitle(paste0(k," PCA"))
      ggsave(paste0(i, "_", k,"_ligandCorrected_3DPCA_by_GroupAndTime.pdf"), width=5, height=4) # was w10 h10
    }

    meta.df.nras2 <- meta.df[rownames(meta.df.nras),]
    #expr.df <- global.df[,c("Gene", rownames(meta.df.nras2))]
    expr.df <- temp.omics.corr[["Global"]][,rownames(meta.df.nras2)]
    expr.df$Gene <- rownames(expr.df)
    cc.df <- meta.df.nras2[,c("Group", contrast2[[i]])]
    temp.expr <- expr.df[grepl("RAS", expr.df$Gene, ignore.case=TRUE), ]
    rownames(temp.expr) <- temp.expr$Gene

    # violin plots of ras genes
    long.temp.expr <- reshape2::melt(temp.expr)
    colnames(long.temp.expr)[2] <- "MeasurementName"
    long.temp.expr <- merge(long.temp.expr, meta.df.nras2, by="MeasurementName")

    Gene <- unique(temp.expr$Gene)
    p.vals <- data.frame(Gene, Lower = NA, Higher = NA, Two.sided = NA)
    for (temp.ras in Gene) {
      p.df <- long.temp.expr[long.temp.expr$Gene == temp.ras,]
      p.df$Group <- factor(p.df$Group, levels=unique(p.df$Group))
      p.df$Time <- factor(p.df$Time, levels=c("0h", "2h", "6h", "24h"))
      temp.plot <- ggplot2::ggplot(p.df, aes(x=Group, y=value)) + geom_violin(alpha=0) +
        geom_point(aes(color=Time, shape=AML_L)) + geom_boxplot(width=0.2, alpha = 0) +
        labs(y=paste("Relative Abundance of",temp.ras)) + theme_classic(base_size = 12) +
        theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),
              axis.title.x = element_blank())

      # t-tests
      temp.nras <- as.numeric(p.df[grepl("NRAS", p.df$Group, ignore.case = TRUE),]$value, na.rm = TRUE)
      temp.other <- as.numeric(p.df[!grepl("NRAS", p.df$Group, ignore.case = TRUE),]$value, na.rm = TRUE)
      temp.higher <- t.test(temp.nras, temp.other, alternative="greater")$p.value
      p.vals[p.vals$Gene == temp.ras,]$Higher <- temp.higher
      temp.lower <- t.test(temp.nras, temp.other, alternative="less")$p.value
      p.vals[p.vals$Gene == temp.ras,]$Lower <- temp.lower
      temp2 <- t.test(temp.nras, temp.other, alternative="two.sided")$p.value
      p.vals[p.vals$Gene == temp.ras,]$Two.sided <- temp2
      if (temp.higher <= 0.05) {
        temp.title <- paste0(temp.ras," expression is\nhigher with NRAS ASO (p = ", signif(temp.higher,2),")")
        temp.plot <- temp.plot + ggtitle(temp.title)
      } else if (temp.lower <= 0.05) {
        temp.title <- paste0(temp.ras," expression is\nlower with NRAS ASO (p = ", signif(temp.lower,2),")")
        temp.plot <- temp.plot + ggtitle(temp.title)
      }
      ggsave(paste0(temp.ras,"_expression.pdf"), temp.plot, width = 5, height = 5)
    }
    write.csv(p.vals, "RAS_expression_pVals.csv", row.names = FALSE)

    if (any(p.vals$Lower <= 0.05) | any(p.vals$Higher <= 0.05) | any(p.vals$Two.sided <= 0.05)) {
      long.p.vals <- reshape2::melt(p.vals)
      colnames(long.p.vals)[2:3] <- c("Direction", "p")
      long.p.vals$minusLogP <- -log10(long.p.vals$p)
      long.p.vals <- long.p.vals[long.p.vals$p <= 0.05,]
      geneOrder <- long.p.vals[order(long.p.vals$p),]$Gene
      ggplot2::ggplot(long.p.vals, aes(x=Gene, y=minusLogP, fill=Direction)) +
        geom_bar(stat="identity", position="dodge") + theme_classic(base_size=12) + labs(y="-Log(P-value)") +
        scale_x_discrete(limits=geneOrder)+
        theme(axis.title.x = element_blank(), axis.text.x=element_text(angle=45,vjust=1,hjust=1))
      ggsave("Altered_RAS_expression_minusLogP.pdf", width=ifelse(0.35*nrow(long.p.vals)<3,3,0.35*nrow(long.p.vals)), height=3)

      dn.genes <- long.p.vals[long.p.vals$Direction=="Lower",]$Gene
      up.genes <- long.p.vals[long.p.vals$Direction=="Higher",]$Gene
      long.p.vals2 <- long.p.vals[long.p.vals$Direction=="Two.sided",]
      if (any(long.p.vals2$Gene %in% dn.genes)) {
        long.p.vals2[long.p.vals2$Gene %in% dn.genes,]$Direction <- "Lower"
      }
      if (any(long.p.vals2$Gene %in% up.genes)) {
        long.p.vals2[long.p.vals2$Gene %in% up.genes,]$Direction <- "Higher"
      }
      geneOrder <- long.p.vals2[order(long.p.vals2$p),]$Gene
      ggplot2::ggplot(long.p.vals2, aes(x=Gene, y=minusLogP, fill=Direction)) +
        geom_bar(stat="identity", position=position_dodge2(preserve="single")
        ) + theme_classic(base_size=12) + labs(y="-Log(P-value)") +
        scale_x_discrete(limits=geneOrder)+
        theme(axis.title.x = element_blank(), axis.text.x=element_text(angle=45,vjust=1,hjust=1))
      ggsave("Altered_RAS_expression_minusLogP_2sided.pdf", width=ifelse(0.35*nrow(long.p.vals)<3,3,0.35*nrow(long.p.vals)), height=3)

      long.p.vals2 <- long.p.vals[long.p.vals$Direction!="Two.sided",]
      geneOrder <- long.p.vals2[order(long.p.vals2$p),]$Gene
      ggplot2::ggplot(long.p.vals2, aes(x=Gene, y=minusLogP, fill=Direction)) +
        geom_bar(stat="identity", position=position_dodge2(preserve="single")
        ) + theme_classic(base_size=12) + labs(y="-Log(P-value)") +
        scale_x_discrete(limits=geneOrder)+
        theme(axis.title.x = element_blank(), axis.text.x=element_text(angle=45,vjust=1,hjust=1))
      ggsave("Altered_RAS_expression_minusLogP_1sided.pdf", width=ifelse(0.35*nrow(long.p.vals)<3,3,0.35*nrow(long.p.vals)), height=3)
    }

    # heatmap of ras genes
    expr.mat <- dplyr::select_if(temp.expr, is.numeric)
    expr.mat <- filter_for_hclust(expr.mat)

    if (nrow(expr.mat) > 1) {
      # create heatmaps
      expr.mat <- as.matrix(expr.mat)
      if (TRUE) {
        temp.heatmap <- try(R.utils::withTimeout(pheatmap::pheatmap(expr.mat, color =
                                                                      colorRampPalette(
                                                                        c("navy", "white", "firebrick3"))(50),
                                                                    cluster_rows = TRUE, cluster_cols = TRUE,
                                                                    scale = "row", annotation_col = cc.df,
                                                                    angle_col = "45",
                                                                    show_colnames = FALSE,
                                                                    fontsize = 10), timeout = 300, onTimeout="error"), silent = TRUE)
        if (!inherits(temp.heatmap, "try-error")) {
          my.clust.heatmaps <- temp.heatmap
        }
      }
      temp.heatmap <- try(R.utils::withTimeout(pheatmap::pheatmap(expr.mat, color =
                                                                    colorRampPalette(
                                                                      c("navy", "white", "firebrick3"))(50),
                                                                  cluster_rows = TRUE, cluster_cols = TRUE,
                                                                  annotation_col = cc.df,
                                                                  angle_col = "45",
                                                                  show_colnames = FALSE,
                                                                  fontsize = 10), timeout = 300, onTimeout="error"), silent = TRUE)
      if (!inherits(temp.heatmap, "try-error")) {
        my.abs.heatmaps <- temp.heatmap
      }
    }
    all.heatmaps <- list("Not_scaled.bp" = my.abs.heatmaps,
                         "Scaled.bp" = my.clust.heatmaps)
    poi.files <- list("Expression_of_RAS_genes.csv" =
                        temp.expr,
                      "Heatmaps" = all.heatmaps)
    save_to_synapse(poi.files)
  } else { # exp 22
    contrasts <- c("NRAS")
    # NRAS vs. NT guide
    aso.treatments <- c("NRAS", "NT guide")
    meta.df.nras <- as.data.frame(meta.df[grepl(aso.treatments[1], meta.df$`Sample Name`) |
                                            grepl(aso.treatments[2], meta.df$`Sample Name`),])
    meta.df.nras$NRAS <- "Untreated"
    meta.df.nras[grepl("NRAS",meta.df.nras$`Sample Name`),]$NRAS <- "Treated"
    rownames(meta.df.nras) <- paste0("X",meta.df.nras$Tube)
    temp.meta[[contrasts[1]]] <- meta.df.nras
    
    meta.df$Tube <- as.numeric(meta.df$Tube)
    meta.df <- meta.df[!is.na(meta.df$Tube),]
    meta.df <- as.data.frame(meta.df)
    rownames(meta.df) <- paste0("X",meta.df$Tube)
    temp.omics <- list("Global" = global.df,
                       "Phospho" = phospho.df)
    meta.df$Treatment <- "BASELINE (No guide)"
    meta.df[grepl("guide", meta.df$`Sample Name`),]$Treatment <- "NT guide"
    meta.df[grepl("NRAS", meta.df$`Sample Name`),]$Treatment <- "NRAS KO"
    meta.df$Volume <- substr(meta.df$Volume, 2, nchar(meta.df$Volume))
    for (k in names(temp.omics)) {
      m <- MSnSet(exprs = temp.omics[[k]][,rownames(meta.df)] %>% as.matrix(),
                            pData = meta.df)
      MSnSet.utils::plot_pca(m, phenotype = "Treatment", label = "Volume") #+ ggtitle(paste0(k," PCA"))
      ggsave(paste0(i, "_", k,"_corrected_PCA_by_", "Treatment", "_wVolumeLabels.pdf"), width=5, height=5)
      
      m_corrected <- correct_batch_effect_NA(m, "Volume", par.prior = T)
      MSnSet.utils::plot_pca(m_corrected, phenotype = "Treatment", label = "Volume") #+ ggtitle(paste0(k," PCA"))
      ggsave(paste0(i, "_", k,"_volCorrected_PCA_by_", "Treatment", "_wVolumeLabels.pdf"), width=5, height=5)
    }

    meta.df.nras2 <- meta.df[rownames(meta.df.nras),]
    expr.df <- global.df[,c("Gene", rownames(meta.df.nras2))]
    cc.df <- meta.df.nras2[,c("Treatment", contrast2[[i]])]
    temp.expr <- expr.df[grepl("RAS", expr.df$Gene, ignore.case=TRUE), ]
    rownames(temp.expr) <- temp.expr$Gene

    # violin plots of ras genes
    long.temp.expr <- reshape2::melt(temp.expr)
    long.temp.expr$variable <- as.character(long.temp.expr$variable)
    long.temp.expr$Tube <- substr(long.temp.expr$variable, 2, nchar(long.temp.expr$variable))
    long.temp.expr <- merge(long.temp.expr, meta.df.nras2, by="Tube")

    Gene <- unique(temp.expr$Gene)
    p.vals <- data.frame(Gene, Lower = NA, Higher = NA, Two.sided = NA)
    for (temp.ras in Gene) {
      p.df <- long.temp.expr[long.temp.expr$Gene == temp.ras,]
      p.df$Treatment <- factor(p.df$Treatment, levels=c("NT guide", "NRAS KO"))
      temp.plot <- ggplot2::ggplot(p.df, aes(x=Treatment, y=value)) + geom_violin(alpha=0) +
        geom_point(aes(color=Volume)) + geom_boxplot(width=0.2, alpha = 0) +
        labs(y=paste("Relative Abundance of",temp.ras),color="L") + theme_classic(base_size = 12) +
        theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),
              axis.title.x = element_blank())

      # t-tests
      temp.nras <- as.numeric(p.df[grepl("NRAS", p.df$Treatment, ignore.case = TRUE),]$value, na.rm = TRUE)
      temp.other <- as.numeric(p.df[!grepl("NRAS", p.df$Treatment, ignore.case = TRUE),]$value, na.rm = TRUE)
      temp.higher <- t.test(temp.nras, temp.other, alternative="greater")$p.value
      p.vals[p.vals$Gene == temp.ras,]$Higher <- temp.higher
      temp.lower <- t.test(temp.nras, temp.other, alternative="less")$p.value
      p.vals[p.vals$Gene == temp.ras,]$Lower <- temp.lower
      temp2 <- t.test(temp.nras, temp.other, alternative="two.sided")$p.value
      p.vals[p.vals$Gene == temp.ras,]$Two.sided <- temp2
      if (temp.higher <= 0.05) {
        temp.title <- paste0(temp.ras," expression is\nhigher with NRAS KO (p = ", signif(temp.higher,2),")")
        temp.plot <- temp.plot + ggtitle(temp.title)
      } else if (temp.lower <= 0.05) {
        temp.title <- paste0(temp.ras," expression is\nlower with NRAS KO (p = ", signif(temp.lower,2),")")
        temp.plot <- temp.plot + ggtitle(temp.title)
      }
      ggsave(paste0(temp.ras,"_expression.pdf"), temp.plot, width = 5, height = 5)
    }
    write.csv(p.vals, "RAS_expression_pVals.csv", row.names = FALSE)

    if (any(p.vals$Lower <= 0.05) | any(p.vals$Higher <= 0.05) | any(p.vals$Two.sided <= 0.05)) {
      long.p.vals <- reshape2::melt(p.vals)
      colnames(long.p.vals)[2:3] <- c("Direction", "p")
      long.p.vals$minusLogP <- -log10(long.p.vals$p)
      long.p.vals <- long.p.vals[long.p.vals$p <= 0.05,]
      geneOrder <- long.p.vals[order(long.p.vals$p),]$Gene
      ggplot2::ggplot(long.p.vals, aes(x=Gene, y=minusLogP, fill=Direction)) +
        geom_bar(stat="identity", position="dodge") + theme_classic(base_size=12) + labs(y="-Log(P-value)") +
        scale_x_discrete(limits=geneOrder)+
        theme(axis.title.x = element_blank(), axis.text.x=element_text(angle=45,vjust=1,hjust=1))
      ggsave("Altered_RAS_expression_minusLogP.pdf", width=ifelse(0.35*nrow(long.p.vals)<3,3,0.35*nrow(long.p.vals)), height=3)

      dn.genes <- long.p.vals[long.p.vals$Direction=="Lower",]$Gene
      up.genes <- long.p.vals[long.p.vals$Direction=="Higher",]$Gene
      long.p.vals2 <- long.p.vals[long.p.vals$Direction=="Two.sided",]
      if (any(long.p.vals2$Gene %in% dn.genes)) {
        long.p.vals2[long.p.vals2$Gene %in% dn.genes,]$Direction <- "Lower"
      }
      if (any(long.p.vals2$Gene %in% up.genes)) {
        long.p.vals2[long.p.vals2$Gene %in% up.genes,]$Direction <- "Higher"
      }
      geneOrder <- long.p.vals2[order(long.p.vals2$p),]$Gene
      ggplot2::ggplot(long.p.vals2, aes(x=Gene, y=minusLogP, fill=Direction)) +
        geom_bar(stat="identity", position=position_dodge2(preserve="single")
        ) + theme_classic(base_size=12) + labs(y="-Log(P-value)") +
        scale_x_discrete(limits=geneOrder)+
        theme(axis.title.x = element_blank(), axis.text.x=element_text(angle=45,vjust=1,hjust=1))
      ggsave("Altered_RAS_expression_minusLogP_2sided.pdf", width=ifelse(0.35*nrow(long.p.vals)<3,3,0.35*nrow(long.p.vals)), height=3)

      long.p.vals2 <- long.p.vals[long.p.vals$Direction!="Two.sided",]
      geneOrder <- long.p.vals2[order(long.p.vals2$p),]$Gene
      ggplot2::ggplot(long.p.vals2, aes(x=Gene, y=minusLogP, fill=Direction)) +
        geom_bar(stat="identity", position=position_dodge2(preserve="single")
        ) + theme_classic(base_size=12) + labs(y="-Log(P-value)") +
        scale_x_discrete(limits=geneOrder)+
        theme(axis.title.x = element_blank(), axis.text.x=element_text(angle=45,vjust=1,hjust=1))
      ggsave("Altered_RAS_expression_minusLogP_1sided.pdf", width=ifelse(0.35*nrow(long.p.vals)<3,3,0.35*nrow(long.p.vals)), height=3)
    }

    # heatmap of ras genes
    expr.mat <- dplyr::select_if(temp.expr, is.numeric)
    expr.mat <- filter_for_hclust(expr.mat)

    if (nrow(expr.mat) > 1) {
      # create heatmaps
      expr.mat <- as.matrix(expr.mat)
      if (TRUE) {
        temp.heatmap <- try(R.utils::withTimeout(pheatmap::pheatmap(expr.mat, color =
                                                                      colorRampPalette(
                                                                        c("navy", "white", "firebrick3"))(50),
                                                                    cluster_rows = TRUE, cluster_cols = TRUE,
                                                                    scale = "row", annotation_col = cc.df,
                                                                    angle_col = "45",
                                                                    show_colnames = FALSE,
                                                                    fontsize = 10), timeout = 300, onTimeout="error"), silent = TRUE)
        if (!inherits(temp.heatmap, "try-error")) {
          my.clust.heatmaps <- temp.heatmap
        }
      }
      temp.heatmap <- try(R.utils::withTimeout(pheatmap::pheatmap(expr.mat, color =
                                                                    colorRampPalette(
                                                                      c("navy", "white", "firebrick3"))(50),
                                                                  cluster_rows = TRUE, cluster_cols = TRUE,
                                                                  annotation_col = cc.df,
                                                                  angle_col = "45",
                                                                  show_colnames = FALSE,
                                                                  fontsize = 10), timeout = 300, onTimeout="error"), silent = TRUE)
      if (!inherits(temp.heatmap, "try-error")) {
        my.abs.heatmaps <- temp.heatmap
      }
    }
    all.heatmaps <- list("Not_scaled.bp" = my.abs.heatmaps,
                         "Scaled.bp" = my.clust.heatmaps)
    poi.files <- list("Expression_of_RAS_genes.csv" =
                        temp.expr,
                      "Heatmaps" = all.heatmaps)
    save_to_synapse(poi.files)
  }

  for (j in contrasts) { #was j in contrasts
    setwd(base.path)
    temp.omics <- list("Global" = global.df[,c("Gene", rownames(temp.meta[[j]]))],
                       "Phospho" = phospho.df[,c("SUB_SITE", rownames(temp.meta[[j]]))])
    # temp.omics <- list("Phospho" = phospho.df[,c("SUB_SITE", rownames(temp.meta[[j]]))])

    # create PCAs
    for (k in names(temp.omics)) {
      m_corrected <- MSnSet(exprs = temp.omics[[k]][,rownames(temp.meta[[j]])] %>% as.matrix(),
                            pData = temp.meta[[j]])
      if (i == "Cell Line ASO") {
        MSnSet.utils::plot_pca(m_corrected, phenotype = j, label = "Time") #+ ggtitle(paste0(k," PCA"))
        ggsave(paste0(i, "_", j, "_", k,"_corrected_PCA_by_", j, "_wTimeLabels.pdf"), width=5, height=5)
        MSnSet.utils::plot_pca(m_corrected, phenotype = j, label = "AML_L") #+ ggtitle(paste0(k," PCA"))
        ggsave(paste0(i, "_", j, "_", k,"_corrected_PCA_by_", j, "_wLLabels.pdf"), width=5, height=5)
      } else if (i == "Cell Line KO") {
        MSnSet.utils::plot_pca(m_corrected, phenotype = j, label = "Volume") #+ ggtitle(paste0(k," PCA"))
        ggsave(paste0(i, "_", j, "_", k,"_corrected_PCA_by_", j, "_wVolumeLabels.pdf"), width=5, height=5)
      } else if (i == "Patient ASO") {
        MSnSet.utils::plot_pca(m_corrected, phenotype = j, label = "PatientID") #+ ggtitle(paste0(k," PCA"))
        ggsave(paste0(i, "_", j, "_", k,"_corrected_PCA_by_", j, "_wPatientIDLabels.pdf"), width=5, height=5)
      }

      for (m in contrast2[[i]]) {
        MSnSet.utils::plot_pca(m_corrected, phenotype = m, label = j) #+ ggtitle(paste0(k," PCA"))
        ggsave(paste0(i, "_", j, "_", k,"_corrected_PCA_by_", m, "_wLabels.pdf"), width=5, height=5)
      }
    }

    # run panSEA
    panSEA2(j, contrast2[[i]], meta.df = temp.meta[[j]], omics = temp.omics,
            annotations = c(), expr = temp.expr, gmt.drug = BeatAML$gmt,
            drug.sens = drug.BeatAML, base.path = base.path, GSEA = FALSE, DMEA = FALSE,
            temp.path = temp.path, synapse_id = synIDs[[i]], plots=FALSE)
  }
}

#### 3. Compile results across exp 20-22 for each omics, contrast ####
base.path <- "~/OneDrive - PNNL/Documents/GitHub/Exp21_NRAS_ASO_treated_patients/proteomics/analysis/"
setwd(base.path)
synapser::synLogin()
source("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/helperScripts/compile_mCorr.R")
source("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/helperScripts/compile_mGSEA.R")
source("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/helperScripts/compile_mDMEA.R")
library(plyr);library(dplyr)
# compile results for:
all.files <- list()
contrasts <- c("NRAS", "NRAS_w_Gilt_ctl", "Gilt", "Gilt_w_NRAS_ctl", "Gilt_w_NRAS")
omics.types <- c("Global", "Phospho")
for (k in contrasts) { # each contrast
  contrast.files <- list()
  for (i in omics.types) { # each omics type
    temp.degs <- list()
    temp.gsea <- list()
    temp.gsea.kegg <- list()
    temp.drug <- list()
    temp.dmea <- list()
    omics.files <- list()
    for (j in names(synIDs)) { # each exp
      synFolders <- unlist(as.list(synapser::synGetChildren(synIDs[[j]], list("folder"), sortBy = 'NAME')))
      if (paste0(k,"_Treated_vs_Untreated") %in% synFolders) {
        # locate contrast folder
        contrastSyn <- synapser::synStore(synapser::Folder(paste0(k,"_Treated_vs_Untreated"), parent = synIDs[[j]])) 
        
        # locate omics folder
        omicsSyn <- synapser::synStore(synapser::Folder(i, parent = contrastSyn)) 
        
        # locate differential expression
        degSyn <- synapser::synStore(synapser::Folder("Differential_expression", parent = omicsSyn))
        degSynFiles <- unlist(as.list(synapser::synGetChildren(degSyn, list("file"), sortBy = 'NAME')))
        degSynResultFile <- degSynFiles[grep("Differential_expression_results.csv", degSynFiles)+1] # synapse ID is immediately after file name
        temp.degs[[j]] <- read.csv(synapser::synGet(degSynResultFile)$path)

        # locate GSEA hallmark or KSEA
        if (grepl("phospho", i, ignore.case = TRUE)) {
          kseaSyn <- synapser::synStore(synapser::Folder("KSEA", parent = omicsSyn))
          kseaSynFiles <- unlist(as.list(synapser::synGetChildren(kseaSyn, list("file"), sortBy = 'NAME')))
          kseaSynResultFile <- kseaSynFiles[grep("KSEA_results.csv", kseaSynFiles)+1] # synapse ID is immediately after file name
          temp.gsea[[j]] <- read.csv(synapser::synGet(kseaSynResultFile)$path)
        } else {
          gseaSyn <- synapser::synStore(synapser::Folder("GSEA", parent = omicsSyn)) 
          hallSyn <- synapser::synStore(synapser::Folder("GSEA_Hallmark", parent = gseaSyn))
          hallSynFiles <- unlist(as.list(synapser::synGetChildren(hallSyn, list("file"), sortBy = 'NAME')))
          if (!is.null(hallSynFiles)) {
            hallSynResultFile <- hallSynFiles[grep("GSEA_results.csv", hallSynFiles)+1] # synapse ID is immediately after file name
            temp.gsea[[j]] <- read.csv(synapser::synGet(hallSynResultFile)$path)
          }
          
          keggSyn <- synapser::synStore(synapser::Folder("GSEA_KEGG", parent = gseaSyn)) 
          keggSynFiles <- unlist(as.list(synapser::synGetChildren(keggSyn, list("file"), sortBy = 'NAME')))
          if (!is.null(keggSynFiles)) {
            keggSynResultFile <- keggSynFiles[grep("GSEA_results.csv", keggSynFiles)+1] # synapse ID is immediately after file name
            temp.gsea.kegg[[j]] <- read.csv(synapser::synGet(keggSynResultFile)$path) 
          }
        }
        
        # locate DMEA
        dmeaSyn <- synapser::synStore(synapser::Folder("DMEA", parent = omicsSyn))
        dmeaSynFiles <- unlist(as.list(synapser::synGetChildren(dmeaSyn, list("file"), sortBy = 'NAME')))
        dmeaSynResultFile <- dmeaSynFiles[grep("DMEA_results.csv", dmeaSynFiles)+1] # synapse ID is immediately after file name
        temp.dmea[[j]] <- read.csv(synapser::synGet(dmeaSynResultFile)$path)
        drugSynResultFile <- dmeaSynFiles[grep("DMEA_correlation_results.csv", dmeaSynFiles)+1] # synapse ID is immediately after file name
        temp.drug[[j]] <- read.csv(synapser::synGet(drugSynResultFile)$path)
      }
    }
    
    # compile DEGs
    if (length(temp.degs) > 1) {
      degResults <- panSEA::compile_mDEG(temp.degs)
      omics.files[["Differential_expression"]] <- list("Differential_expression_results.csv" =
                                                         degResults$results,
                                                       "Differential_expression_mean_results.csv" =
                                                         degResults$mean.results,
                                                       "Differential_expression_venn_diagram.pdf" =
                                                         degResults$venn.diagram,
                                                       "Differential_expression_correlation_matrix.pdf" =
                                                         degResults$corr.matrix,
                                                       "Differential_expression_dot_plot.pdf" =
                                                         degResults$dot.plot)
    }

    # compile GSEA results
    if (length(temp.gsea) > 1) {
      gseaResults <- compile_mGSEA(temp.gsea)
      if (grepl("phospho", i, ignore.case=TRUE)) {
        omics.files[["KSEA"]] <- list("KSEA_results.csv" =
                                        gseaResults$results,
                                      "KSEA_mean_results.csv" =
                                        gseaResults$mean.results,
                                      "KSEA_venn_diagram.pdf" =
                                        gseaResults$venn.diagram,
                                      "KSEA_correlation_matrix.pdf" =
                                        gseaResults$corr.matrix,
                                      "KSEA_dot_plot.pdf" =
                                        gseaResults$dot.plot)
      } else {
        omics.files[["GSEA"]] <- list("GSEA_results.csv" =
                                        gseaResults$results,
                                      "GSEA_mean_results.csv" =
                                        gseaResults$mean.results,
                                      "GSEA_venn_diagram.pdf" =
                                        gseaResults$venn.diagram,
                                      "GSEA_correlation_matrix.pdf" =
                                        gseaResults$corr.matrix,
                                      "GSEA_dot_plot.pdf" =
                                        gseaResults$dot.plot)
      }
    }
    
    if (length(temp.gsea.kegg) > 1) {
      gseaResults <- compile_mGSEA(temp.gsea.kegg)
      omics.files[["GSEA_KEGG"]] <- list("GSEA_results.csv" =
                                      gseaResults$results,
                                    "GSEA_mean_results.csv" =
                                      gseaResults$mean.results,
                                    "GSEA_venn_diagram.pdf" =
                                      gseaResults$venn.diagram,
                                    "GSEA_correlation_matrix.pdf" =
                                      gseaResults$corr.matrix,
                                    "GSEA_dot_plot.pdf" =
                                      gseaResults$dot.plot) 
    }
    
    # compile DMEA results
    if (length(temp.drug) > 1 & length(temp.dmea) > 1) {
      drugResults <- compile_mCorr(temp.drug)
      drug.files <- list("DMEA_correlation_results.csv" =
                           drugResults$results,
                         "DMEA_mean_correlation_results.csv" =
                           drugResults$mean.results,
                         "DMEA_correlation_venn_diagram.pdf" =
                           drugResults$venn.diagram,
                         "DMEA_correlation_correlation_matrix.pdf" =
                           drugResults$corr.matrix,
                         "DMEA_correlation_dot_plot.pdf" =
                           drugResults$dot.plot)
      dmeaResults <- compile_mDMEA(temp.dmea)
      omics.files[["DMEA"]] <- list("DMEA_results.csv" =
                                      dmeaResults$results,
                                    "DMEA_mean_results.csv" =
                                      dmeaResults$mean.results,
                                    "DMEA_venn_diagram.pdf" =
                                      dmeaResults$venn.diagram,
                                    "DMEA_correlation_matrix.pdf" =
                                      dmeaResults$corr.matrix,
                                    "DMEA_dot_plot.pdf" =
                                      dmeaResults$dot.plot,
                                    "Drug_correlations" = drug.files)
    }
    
    # save results
    contrast.files[[i]] <- omics.files
  } 
  all.files[[k]] <- contrast.files
}
dir.create("Compiled_results_exp20-22")
setwd("Compiled_results_exp20-22")
compiledSyn <- synapser::synStore(synapser::Folder("Compiled_results_exp20-22",
                                                   parent = synIDs[["Patient ASO"]])) 
save_to_synapse(all.files, compiledSyn)

#### 4. redo panSEA for exp 20 for each timepoint ####
setwd("~/OneDrive - PNNL/Documents/GitHub/Exp21_NRAS_ASO_treated_patients/proteomics/")
dir.create("analysis")
setwd("analysis")
base.path <- "~/OneDrive - PNNL/Documents/GitHub/Exp21_NRAS_ASO_treated_patients/proteomics/analysis/"

cell <- list("Meta" = "syn52653527", # exp 20
             "Global" = "syn52623904",
             "Phospho" = "syn52623908")
inputs <- list("Cell Line ASO by Time" = cell)
meta.sheet <- c(1)
names(meta.sheet) <- names(inputs)
synIDs <- c("syn52256487")
names(synIDs) <- names(inputs)
contrast2 <- list("Cell Line ASO by Time" = c("AML_L"))

for (i in names(inputs)) {
  setwd(base.path)
  dir.create(i)
  
  ## read in data
  meta.df <- readxl::read_excel(synapser::synGet(inputs[[i]]$Meta)$path, 
                                sheet = meta.sheet[[i]])
  global.df <- read.table(synapser::synGet(inputs[[i]]$Global)$path, sep = "\t")
  phospho.df <- read.table(synapser::synGet(inputs[[i]]$Phospho)$path, sep = "\t")
  
  # make sure contrasts are valid
  for (j in contrast2[[i]]) {
    for (k in 1:nrow(meta.df)) {
      meta.df[k,j] <- make.names(meta.df[k,j])
    }
  }
  
  # add column for feature names and later make it the first column
  #global.df$Gene <- rownames(global.df)
  phospho.df$SUB_SITE <- rownames(phospho.df)
  
  ## prep contrasts
  timepoints <- na.omit(unique(meta.df$Time[meta.df$Time != "baseline"]))
  for (k in timepoints) {
    setwd(file.path(base.path, i))
    dir.create(k)
    setwd(k)
    temp.path <- file.path(base.path,i, k)
    timeSyn <- synapser::synStore(synapser::Folder(k, parent = synIDs[[i]])) 
    
    temp.meta.df <- meta.df[meta.df$Time == k,]
    temp.meta <- list()
    contrasts <- c("NRAS", "NRAS_w_Gilt_ctl", "Gilt", "Gilt_w_NRAS_ctl", "Gilt_w_NRAS")
    # NRAS ASO vs. NEG ASO
    aso.treatments <- c("NRAS ASO", "CTRL ASO")
    meta.df.nras <- as.data.frame(temp.meta.df[temp.meta.df$Group %in% aso.treatments,])
    meta.df.nras$NRAS <- "Untreated"
    meta.df.nras[grepl("NRAS ASO", meta.df.nras$Sample),]$NRAS <- "Treated"
    rownames(meta.df.nras) <- meta.df.nras$MeasurementName
    temp.meta[[contrasts[1]]] <- meta.df.nras

    # Gilt + NRAS vs. Gilt
    nras.gilt.treatments <- c("NRAS ASO + GILT", "CTRL ASO + GILT")
    meta.df.nras.gilt <- as.data.frame(temp.meta.df[temp.meta.df$Group %in% nras.gilt.treatments,])
    meta.df.nras.gilt$NRAS_w_Gilt_ctl <- "Untreated"
    meta.df.nras.gilt[meta.df.nras.gilt$Group == "NRAS ASO + GILT",]$NRAS_w_Gilt_ctl <- "Treated"
    rownames(meta.df.nras.gilt) <- meta.df.nras.gilt$MeasurementName
    temp.meta[[contrasts[2]]] <- meta.df.nras.gilt

    # Gilt + Neg ASO vs. Neg ASO
    gilt.treatments <- c("CTRL ASO + GILT", "CTRL ASO")
    meta.df.gilt <- as.data.frame(temp.meta.df[temp.meta.df$Group %in% gilt.treatments,])
    meta.df.gilt$Gilt <- "Untreated"
    meta.df.gilt[meta.df.gilt$Group == "CTRL ASO + GILT",]$Gilt <- "Treated"
    rownames(meta.df.gilt) <- meta.df.gilt$MeasurementName
    temp.meta[[contrasts[3]]] <- meta.df.gilt
    
    # Gilt + NRAS vs. NRAS
    nras.gilt.treatments <- c("NRAS ASO + GILT", "NRAS ASO")
    meta.df.nras.gilt <- as.data.frame(temp.meta.df[temp.meta.df$Group %in% nras.gilt.treatments,])
    meta.df.nras.gilt$Gilt_w_NRAS_ctl <- "Untreated"
    meta.df.nras.gilt[meta.df.nras.gilt$Group == "NRAS ASO + GILT",]$Gilt_w_NRAS_ctl <- "Treated"
    rownames(meta.df.nras.gilt) <- meta.df.nras.gilt$MeasurementName
    temp.meta[[contrasts[4]]] <- meta.df.nras.gilt
    
    # Gilt + NRAS
    nras.gilt.treatments <- c("NRAS ASO + GILT", "CTRL ASO")
    meta.df.nras.gilt <- as.data.frame(temp.meta.df[temp.meta.df$Group %in% nras.gilt.treatments,])
    meta.df.nras.gilt$Gilt_w_NRAS <- "Untreated"
    meta.df.nras.gilt[meta.df.nras.gilt$Group == "NRAS ASO + GILT",]$Gilt_w_NRAS <- "Treated"
    rownames(meta.df.nras.gilt) <- meta.df.nras.gilt$MeasurementName
    temp.meta[[contrasts[5]]] <- meta.df.nras.gilt
    
    
    # run panSEA for Gilt effect
    for (j in contrasts) {
      setwd(base.path)
      temp.omics <- list("Global" = global.df[,c("Gene", rownames(temp.meta[[j]]))],
                         "Phospho" = phospho.df[,c("SUB_SITE", rownames(temp.meta[[j]]))])
      # temp.omics <- list("Phospho" = phospho.df[,c("SUB_SITE", rownames(temp.meta[[j]]))])
      panSEA2(j, contrast2[[i]], meta.df = temp.meta[[j]], omics = temp.omics, 
              annotations = c(), expr = temp.expr, gmt.drug = BeatAML$gmt, 
              drug.sens = drug.BeatAML, base.path = base.path, #GSEA=FALSE, DMEA=FALSE, plots=FALSE,
              temp.path = temp.path, synapse_id = timeSyn) 
    } 
  }
}

#### 5. compile exp 20 results over time for each contrast ####
base.path <- "~/OneDrive - PNNL/Documents/GitHub/Exp21_NRAS_ASO_treated_patients/proteomics/analysis/"
setwd(base.path)
synapser::synLogin()
source("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/helperScripts/compile_mCorr.R")
source("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/helperScripts/compile_mGSEA.R")
source("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/helperScripts/compile_mDMEA.R")
source("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/helperScripts/compile_mDEG.R")
library(plyr);library(dplyr)
# compile results for:
all.files <- list()
contrasts <- c("NRAS", "NRAS_w_Gilt_ctl", "Gilt", "Gilt_w_NRAS_ctl", "Gilt_w_NRAS")
omics.types <- c("Global", "Phospho")
timepoints <- c("2 HR", "6 HR", "24 HR") # impose sequential order
for (k in contrasts) { # each contrast
  contrast.files <- list()
  for (i in omics.types) { # each omics type
    temp.degs <- list()
    temp.gsea <- list()
    temp.gsea.kegg <- list()
    temp.drug <- list()
    temp.dmea <- list()
    omics.files <- list()
    for (j in timepoints) { # each exp
      timeSyn <- synapser::synStore(synapser::Folder(j, parent = synIDs[["Cell Line ASO by Time"]]))
      timeSynFolders <- unlist(as.list(synapser::synGetChildren(timeSyn, list("folder"), sortBy = 'NAME')))
      if (paste0(k,"_Treated_vs_Untreated") %in% synFolders) {
        # locate contrast folder
        contrastSyn <- synapser::synStore(synapser::Folder(paste0(k,"_Treated_vs_Untreated"), parent = timeSyn)) 
        
        # locate omics folder
        omicsSyn <- synapser::synStore(synapser::Folder(i, parent = contrastSyn)) 
        
        # locate differential expression
        degSyn <- synapser::synStore(synapser::Folder("Differential_expression", parent = omicsSyn))
        degSynFiles <- unlist(as.list(synapser::synGetChildren(degSyn, list("file"), sortBy = 'NAME')))
        degSynResultFile <- degSynFiles[grep("Differential_expression_results.csv", degSynFiles)+1] # synapse ID is immediately after file name
        temp.degs[[j]] <- read.csv(synapser::synGet(degSynResultFile)$path)

        # locate GSEA hallmark or KSEA
        if (grepl("phospho", i, ignore.case = TRUE)) {
          kseaSyn <- synapser::synStore(synapser::Folder("KSEA", parent = omicsSyn))
          kseaSynFiles <- unlist(as.list(synapser::synGetChildren(kseaSyn, list("file"), sortBy = 'NAME')))
          kseaSynResultFile <- kseaSynFiles[grep("KSEA_results.csv", kseaSynFiles)+1] # synapse ID is immediately after file name
          temp.gsea[[j]] <- read.csv(synapser::synGet(kseaSynResultFile)$path)
        } else {
          gseaSyn <- synapser::synStore(synapser::Folder("GSEA", parent = omicsSyn)) 
          hallSyn <- synapser::synStore(synapser::Folder("GSEA_Hallmark", parent = gseaSyn))
          hallSynFiles <- unlist(as.list(synapser::synGetChildren(hallSyn, list("file"), sortBy = 'NAME')))
          if (!is.null(hallSynFiles)) {
            hallSynResultFile <- hallSynFiles[grep("GSEA_results.csv", hallSynFiles)+1] # synapse ID is immediately after file name
            temp.gsea[[j]] <- read.csv(synapser::synGet(hallSynResultFile)$path)
          }
          
          keggSyn <- synapser::synStore(synapser::Folder("GSEA_KEGG", parent = gseaSyn)) 
          keggSynFiles <- unlist(as.list(synapser::synGetChildren(keggSyn, list("file"), sortBy = 'NAME')))
          if (!is.null(keggSynFiles)) {
            keggSynResultFile <- keggSynFiles[grep("GSEA_results.csv", keggSynFiles)+1] # synapse ID is immediately after file name
            temp.gsea.kegg[[j]] <- read.csv(synapser::synGet(keggSynResultFile)$path) 
          }
        }
        
        # locate DMEA
        dmeaSyn <- synapser::synStore(synapser::Folder("DMEA", parent = omicsSyn))
        dmeaSynFiles <- unlist(as.list(synapser::synGetChildren(dmeaSyn, list("file"), sortBy = 'NAME')))
        dmeaSynResultFile <- dmeaSynFiles[grep("DMEA_results.csv", dmeaSynFiles)+1] # synapse ID is immediately after file name
        temp.dmea[[j]] <- read.csv(synapser::synGet(dmeaSynResultFile)$path)
        drugSynResultFile <- dmeaSynFiles[grep("DMEA_correlation_results.csv", dmeaSynFiles)+1] # synapse ID is immediately after file name
        temp.drug[[j]] <- read.csv(synapser::synGet(drugSynResultFile)$path)
      }
    }
    
    # compile DEGs
    if (length(temp.degs) > 1) {
      degResults <- panSEA::compile_mDEG(temp.degs)
      omics.files[["Differential_expression"]] <- list("Differential_expression_results.csv" =
                                                         degResults$results,
                                                       "Differential_expression_mean_results.csv" =
                                                         degResults$mean.results,
                                                       "Differential_expression_venn_diagram.pdf" =
                                                         degResults$venn.diagram,
                                                       "Differential_expression_correlation_matrix.pdf" =
                                                         degResults$corr.matrix,
                                                       "Differential_expression_dot_plot.pdf" =
                                                         degResults$dot.plot)
    }

    # compile GSEA results
    if (length(temp.gsea) > 1) {
      gseaResults <- compile_mGSEA(temp.gsea)
      if (grepl("phospho", i, ignore.case=TRUE)) {
        omics.files[["KSEA"]] <- list("KSEA_results.csv" =
                                        gseaResults$results,
                                      "KSEA_mean_results.csv" =
                                        gseaResults$mean.results,
                                      "KSEA_venn_diagram.pdf" =
                                        gseaResults$venn.diagram,
                                      "KSEA_correlation_matrix.pdf" =
                                        gseaResults$corr.matrix,
                                      "KSEA_dot_plot.pdf" =
                                        gseaResults$dot.plot)
      } else {
        omics.files[["GSEA"]] <- list("GSEA_results.csv" =
                                        gseaResults$results,
                                      "GSEA_mean_results.csv" =
                                        gseaResults$mean.results,
                                      "GSEA_venn_diagram.pdf" =
                                        gseaResults$venn.diagram,
                                      "GSEA_correlation_matrix.pdf" =
                                        gseaResults$corr.matrix,
                                      "GSEA_dot_plot.pdf" =
                                        gseaResults$dot.plot)
      }
    }
    
    if (length(temp.gsea.kegg) > 1) {
      gseaResults <- compile_mGSEA(temp.gsea.kegg)
      omics.files[["GSEA_KEGG"]] <- list("GSEA_results.csv" =
                                           gseaResults$results,
                                         "GSEA_mean_results.csv" =
                                           gseaResults$mean.results,
                                         "GSEA_venn_diagram.pdf" =
                                           gseaResults$venn.diagram,
                                         "GSEA_correlation_matrix.pdf" =
                                           gseaResults$corr.matrix,
                                         "GSEA_dot_plot.pdf" =
                                           gseaResults$dot.plot) 
    }
    
    # compile DMEA results
    if (length(temp.drug) > 1 & length(temp.dmea) > 1) {
      drugResults <- compile_mCorr(temp.drug)
      drug.files <- list("DMEA_correlation_results.csv" =
                           drugResults$results,
                         "DMEA_mean_correlation_results.csv" =
                           drugResults$mean.results,
                         "DMEA_correlation_venn_diagram.pdf" =
                           drugResults$venn.diagram,
                         "DMEA_correlation_correlation_matrix.pdf" =
                           drugResults$corr.matrix,
                         "DMEA_correlation_dot_plot.pdf" =
                           drugResults$dot.plot)
      dmeaResults <- compile_mDMEA(temp.dmea)
      omics.files[["DMEA"]] <- list("DMEA_results.csv" =
                                      dmeaResults$results,
                                    "DMEA_mean_results.csv" =
                                      dmeaResults$mean.results,
                                    "DMEA_venn_diagram.pdf" =
                                      dmeaResults$venn.diagram,
                                    "DMEA_correlation_matrix.pdf" =
                                      dmeaResults$corr.matrix,
                                    "DMEA_dot_plot.pdf" =
                                      dmeaResults$dot.plot,
                                    "Drug_correlations" = drug.files)
    }
    
    # save results
    contrast.files[[i]] <- omics.files
  } 
  all.files[[k]] <- contrast.files
}
setwd(base.path)
dir.create("Compiled_exp20_results_over_time")
setwd("Compiled_exp20_results_over_time")
compiledSyn <- synapser::synStore(synapser::Folder("Compiled_results_over_time",
                                                   parent = synIDs[["Cell Line ASO by Time"]])) 
save_to_synapse(all.files, compiledSyn)

#### plot drugs for NRAS KO ####
# load data
base.path <- "~/OneDrive - PNNL/Documents/GitHub/Exp21_NRAS_ASO_treated_patients/proteomics/analysis/"
setwd(base.path)
dir.create("Compiled_exp20_results_over_time")
setwd("Compiled_exp20_results_over_time")
setwd("NRAS/Global/DMEA/Drug_correlations")
drug.corr.df <- read.csv("DMEA_correlation_results.csv")
colnames(drug.corr.df)[2] <- "Drug"
drug.info <- read.csv("~/OneDrive - PNNL/Documents/PTRC2/BeatAML_single_drug_moa.csv",
                      stringsAsFactors = FALSE, fileEncoding = "latin1")
drug.info <- drug.info[,c("Drug","moa")]
drug.info[drug.info$Drug == "Ralimetinib (LY2228820)",]$moa <- "p38 MAPK inhibitor"
drug.info[drug.info$Drug == "Nilotinib",]$moa <- "Abl kinase inhibitor"
drug.info[drug.info$Drug == "AT-101",]$moa <- "BCL inhibitor"
drug.info[is.na(drug.info$Drug),]$moa <- "Other"
library(patchwork); library(ggplot2)
dir.create("negativelyCorr")
setwd("negativelyCorr")
pearson.plots <- NULL
spearman.plots <- NULL
#MOAsInTop50 <- names(gmt.drug$genesets)
#moaColors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(length(MOAsInTop50))
pearson.venn <- list()
spearman.venn <- list()
for (i in unique(drug.corr.df$type)) {
  p.df <- drug.corr.df[drug.corr.df$type == i,]
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
        plot.df <- plot.df[plot.df$Pearson.est < 0 & plot.df$Pearson.q <= 0.05,]
        ylab <- paste0(descr," r")
        pearson.venn[[i]] <- unique(plot.df$Drug)
      } else {
        plot.df <- plot.df[plot.df$Spearman.est < 0 & plot.df$Spearman.q <= 0.05,]
        ylab <- paste0(descr," rho")
        spearman.venn[[i]] <- unique(plot.df$Drug)
      }
      plot.df$rank <- plot.df[,j]
      sigOrder <- na.omit(unique(plot.df[order(plot.df$rank, decreasing=TRUE),]$Drug))
      plot.annot <- paste0(i, "\n(", nrow(plot.df), " / ", nrow(p.df), " Drugs Negatively Correlated)")
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
ggplot2::ggsave("Drug_DIA_WV_moaFill_Pearson_allSigs.pdf", pearson.plots, width=10, height=8)
ggplot2::ggsave("Drug_DIA_WV_moaFill_Spearman_allSigs.pdf", spearman.plots, width=10, height=8)
ggvenn::ggvenn(pearson.venn, show_percentage=FALSE, set_name_size=5, text_size=5)
ggsave("Drug_DIA_WV_Pearson_sigOverlap.pdf", width=5, height=5)
ggvenn::ggvenn(spearman.venn, show_percentage=FALSE, set_name_size=5, text_size=5)
ggsave("Drug_DIA_WV_Spearman_sigOverlap.pdf", width=5, height=5)

# load data
base.path <- "~/OneDrive - PNNL/Documents/GitHub/Exp21_NRAS_ASO_treated_patients/proteomics/analysis/"
setwd(base.path)
dir.create("Compiled_results_exp20-22")
setwd("Compiled_results_exp20-22")
setwd("NRAS/Global/DMEA/Drug_correlations")
drug.corr.df <- read.csv("DMEA_correlation_results.csv")
colnames(drug.corr.df)[2] <- "Drug"
drug.info <- read.csv("~/OneDrive - PNNL/Documents/PTRC2/BeatAML_single_drug_moa.csv",
                      stringsAsFactors = FALSE, fileEncoding = "latin1")
drug.info <- drug.info[,c("Drug","moa")]
drug.info[drug.info$Drug == "Ralimetinib (LY2228820)",]$moa <- "p38 MAPK inhibitor"
drug.info[drug.info$Drug == "Nilotinib",]$moa <- "Abl kinase inhibitor"
drug.info[drug.info$Drug == "AT-101",]$moa <- "BCL inhibitor"
drug.info[is.na(drug.info$Drug),]$moa <- "Other"
library(patchwork); library(ggplot2)
dir.create("negativelyCorr")
setwd("negativelyCorr")
pearson.plots <- NULL
spearman.plots <- NULL
#MOAsInTop50 <- names(gmt.drug$genesets)
#moaColors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(length(MOAsInTop50))
pearson.venn <- list()
spearman.venn <- list()
for (i in unique(drug.corr.df$type)) {
  p.df <- drug.corr.df[drug.corr.df$type == i,]
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
        plot.df <- plot.df[plot.df$Pearson.est < 0 & plot.df$Pearson.q <= 0.05,]
        ylab <- paste0(descr," r")
        pearson.venn[[i]] <- unique(plot.df$Drug)
      } else {
        plot.df <- plot.df[plot.df$Spearman.est < 0 & plot.df$Spearman.q <= 0.05,]
        ylab <- paste0(descr," rho")
        spearman.venn[[i]] <- unique(plot.df$Drug)
      }
      plot.df$rank <- plot.df[,j]
      sigOrder <- na.omit(unique(plot.df[order(plot.df$rank, decreasing=TRUE),]$Drug))
      plot.annot <- paste0(i, "\n(", nrow(plot.df), " / ", nrow(p.df), " Drugs Negatively Correlated)")
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
ggplot2::ggsave("Drug_DIA_WV_moaFill_Pearson_allSigs.pdf", pearson.plots, width=12, height=8)
ggplot2::ggsave("Drug_DIA_WV_moaFill_Spearman_allSigs.pdf", spearman.plots, width=12, height=8)
ggvenn::ggvenn(pearson.venn, show_percentage=FALSE, set_name_size=5, text_size=5)
ggsave("Drug_DIA_WV_Pearson_sigOverlap.pdf", width=5, height=5)
ggvenn::ggvenn(spearman.venn, show_percentage=FALSE, set_name_size=5, text_size=5)
ggsave("Drug_DIA_WV_Spearman_sigOverlap.pdf", width=5, height=5)

##### plot drugs for NRAS KO - highlight BCL2, FLT3i ####
# load data
base.path <- "~/OneDrive - PNNL/Documents/GitHub/Exp21_NRAS_ASO_treated_patients/proteomics/analysis/"
setwd(base.path)
dir.create("Compiled_results_exp20-22")
setwd("Compiled_results_exp20-22")
setwd("NRAS/Global/DMEA/Drug_correlations")
drug.corr.df <- read.csv("DMEA_correlation_results.csv")
colnames(drug.corr.df)[2] <- "Drug"
drug.corr.df[drug.corr.df$Drug == "NF-kB Activation Inhibitor",]$Drug <- "NFkB Inhibitor"
drug.corr.df$Drug <- sub(" [(].*", "", drug.corr.df$Drug) # shorten drug names for plot
drug.info <- read.csv("~/OneDrive - PNNL/Documents/PTRC2/BeatAML_single_drug_moa.csv",
                      stringsAsFactors = FALSE, fileEncoding = "latin1")
drug.info[drug.info$Drug == "NF-kB Activation Inhibitor",]$Drug <- "NFkB Inhibitor"
drug.info$Drug <- sub(" [(].*", "", drug.info$Drug) # shorten drug names for plot
BCL2i <- unique(drug.info[grepl("BCL", drug.info$moa),]$Drug)
FLT3i <- unique(drug.info[grepl("FLT", drug.info$moa),]$Drug)
drug.info <- drug.info[,c("Drug","moa")]
library(patchwork); library(ggplot2); library(plyr); library(dplyr)
dir.create("negativelyCorr")
setwd("negativelyCorr")
pearson.plots <- NULL
spearman.plots <- NULL
pearson.venn <- list()
spearman.venn <- list()
for (i in unique(drug.corr.df$type)) {
  p.df <- drug.corr.df[drug.corr.df$type == i,]
  p.df$Pearson.q <- qvalue::qvalue(p = p.df$Pearson.p, pi0 = 1)$qvalues
  p.df$Spearman.q <- qvalue::qvalue(p = p.df$Spearman.p, pi0 = 1)$qvalues
  p.df <- merge(p.df, drug.info, by="Drug", all.x = TRUE)
  p.df$Mechanism <- "Other"
  for (j in 1:nrow(p.df)) {
    temp.drugs <- strsplit(p.df$Drug[j], " - ")[[1]]
    if (any(temp.drugs %in% BCL2i) & any(temp.drugs %in% FLT3i)) {
      p.df$Mechanism[j] <- "BCL2 & FLT3 inhibitor"
    } else if (any(temp.drugs %in% BCL2i)) {
      p.df$Mechanism[j] <- "BCL2 inhibitor"
    } else if (any(temp.drugs %in% FLT3i)) {
      p.df$Mechanism[j] <- "FLT3 inhibitor"
    }
  }
  
  rank.metrics <- c("Pearson.est", "Spearman.est")
  for (j in rank.metrics) {
    descr <- stringr::str_split_1(j, "[.]")[1]
    if ("Drug" %in% colnames(p.df)) {
      if (j == "Pearson.est") {
        plot.df <- p.df[p.df$Pearson.est < 0 & p.df$Pearson.q <= 0.05,]
        ylab <- paste0(descr," r")
        pearson.venn[[i]] <- unique(plot.df$Drug)
      } else {
        plot.df <- p.df[p.df$Spearman.est < 0 & p.df$Spearman.q <= 0.05,]
        ylab <- paste0(descr," rho")
        spearman.venn[[i]] <- unique(plot.df$Drug)
      }
      plot.df$rank <- plot.df[,j]
      sigOrder <- unique(plot.df[order(plot.df$rank),]$Drug)
      plot.annot <- paste0(i, "\n(", nrow(plot.df), " / ", nrow(p.df), " Drugs Negatively Correlated)")
      corr.plot <- ggplot(plot.df, aes(x=Drug, y=rank, fill = Mechanism)) + 
        geom_col() + theme_minimal(base_size = 12) + ylab(ylab) + 
        ggplot2::scale_x_discrete(limits = sigOrder) +
        theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1),
              axis.title.x=element_blank()) +
        ggtitle(plot.annot) + 
        theme(plot.title = element_text(hjust = 0.5, face="bold", size=16), legend.position="bottom")
      ggsave(paste0("Drug_DIA_WV_moaFill_",descr,"_", i, "_BCL2-FTL3i.pdf"), corr.plot, width = 12, height = 3.25)
      
      top20.df <- plot.df %>% slice_min(Pearson.est, n=20)
      corr.plot <- ggplot(top20.df, aes(x=Drug, y=rank, fill = Mechanism)) + 
        geom_col() + theme_minimal(base_size = 12) + ylab(ylab) + 
        ggplot2::scale_x_discrete(limits = sigOrder[sigOrder %in% top20.df$Drug]) +
        theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1),
              axis.title.x=element_blank()) +
        ggtitle(plot.annot) + 
        theme(plot.title = element_text(hjust = 0.5, face="bold", size=16), legend.position="bottom")
      ggsave(paste0("Drug_DIA_WV_moaFill_",descr,"_", i, "_BCL2-FTL3i_top20.pdf"), corr.plot, width = 5, height = 3.25)
      
      top20.df <- plot.df %>% slice_min(Pearson.est, n=10)
      corr.plot <- ggplot(top20.df, aes(x=Drug, y=rank, fill = Mechanism)) + 
        geom_col() + theme_minimal(base_size = 12) + ylab(ylab) + 
        ggplot2::scale_x_discrete(limits = sigOrder[sigOrder %in% top20.df$Drug]) +
        theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1),
              axis.title.x=element_blank()) +
        ggtitle(plot.annot) + 
        theme(plot.title = element_text(hjust = 0.5, face="bold", size=16), legend.position="bottom")
      ggsave(paste0("Drug_DIA_WV_moaFill_",descr,"_", i, "_BCL2-FTL3i_top10.pdf"), corr.plot, width = 3, height = 3.25)
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
ggplot2::ggsave("Drug_DIA_WV_moaFill_Pearson_allSigs_BCL2-FTL3i.pdf", pearson.plots, width=12, height=8)
ggplot2::ggsave("Drug_DIA_WV_moaFill_Spearman_allSigs_BCL2-FTL3i.pdf", spearman.plots, width=12, height=8)
ggvenn::ggvenn(pearson.venn, show_percentage=FALSE, set_name_size=5, text_size=5)
ggsave("Drug_DIA_WV_Pearson_sigOverlap_BCL2-FTL3i.pdf", width=5, height=5)
ggvenn::ggvenn(spearman.venn, show_percentage=FALSE, set_name_size=5, text_size=5)
ggsave("Drug_DIA_WV_Spearman_sigOverlap_BCL2-FTL3i.pdf", width=5, height=5)
