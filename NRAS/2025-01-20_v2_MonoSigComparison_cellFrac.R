# score samples using various mono vs. prog signatures

library(synapser)
library(DMEA)
library(tidyr)
library(tibble)
library(reshape2)
library(plyr)
library(dplyr)
synapser::synLogin()

evalOneMonoSig <- function(global.df100, temp.sig) {
  temp.sig <- na.omit(temp.sig)
  if (nrow(temp.sig) > 0) {
    # perform weighted voting
    temp.sig2 <- temp.sig[temp.sig$Gene %in% colnames(global.df100)[2:ncol(global.df100)],]
    if (nrow(temp.sig2) > 0) {
      global.df100 <- global.df100[,c("Sample",temp.sig2$Gene)]
      if (nrow(global.df100) > 0) {
        sorted.wv <- panSEA::WV(global.df100, temp.sig2)
        temp.wv <- sorted.wv$scores 
      }
    } else {
      temp.wv <- data.frame()
    }
  } else {
    temp.wv <- data.frame()
  }
  
  return(temp.wv)
}

evalMonoSig <- function(global.df100, sig.matrix) {
  # evaluate each signature
  wv.df <- data.frame()
  for (i in 2:ncol(sig.matrix)) {
    cat("evaluating",names(sig.matrix)[i],"\n")
    temp.wv.df <- evalOneMonoSig(global.df100, sig.matrix[,c("Gene",names(sig.matrix)[i])])
    temp.wv.df$Signature <- names(sig.matrix)[i]
    wv.df <- rbind(wv.df, temp.wv.df)
  }
  
  return(wv.df)
}

compareSigs <- function(global.df100, sigs, value.var = "Log2FC",  
                        fillVals = RColorBrewer::brewer.pal(length(sigs), "Set2")) {
  # combine signatures into matrix
  filtered.sigs.df <- data.table::rbindlist(sigs, use.names = TRUE, idcol = "Signature")
  sig.matrix <- reshape2::dcast(filtered.sigs.df, Gene ~ Signature, mean,
                                value.var = value.var)
  rownames(sig.matrix) <- sig.matrix$Gene
  sig.matrix <- sig.matrix[,c("Gene",names(sigs))]
  
  # test signature matrix
  wv.df <- evalMonoSig(global.df100, sig.matrix)
  
  return(wv.df)
}

# load sorted proteomics signature
sig.paths <- list("Sorted" = "analysis/combined24-27/DIA_2batches_noOutliers_noMSC/Sort Type_Bead/CD14_Pos_vs_Neg/global/Differential_expression/Differential_expression_results.csv",
                  "van Galen" = "data/externalSignatures/formatted/notFilteredForMalignant/Differential_expression_van_Galen_AML_D0_Mono-like_vs_Prog-like_noNA.csv",
                  "Triana" = "data/externalSignatures/formatted/Triana_RNA_AML_100PercentCells_Classical-Monocytes_vs_HSCs-and-MPPs_differentialExpression.csv",
                  "Lasry" = "data/externalSignatures/formatted/notFilteredForMalignant/Differential_expression_Lasry_AML_CD14PosMonocyte_vs_HSC_protein-coding.csv")

# import signatures and filter
sigs <- list()
setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/")
for (i in names(sig.paths)) {
  sigs[[i]] <- read.csv(sig.paths[[i]])
  sigs[[i]] <- na.omit(sigs[[i]][sigs[[i]]$adj.P.Val <= 0.05,c("Gene","Log2FC")])
}

#### look at van Galen monocyte and progenitor markers ####
### get van Galen markers
synapser::synLogin()
vg <- readxl::read_excel(synapser::synGet("syn61831560")$path)
vg.tumor <- vg[,8:(ncol(vg)-1)]
colnames(vg.tumor) <- vg.tumor[1,]
vg.tumor <- vg.tumor[2:(nrow(vg.tumor)-2),]

# prepare progenitor markers
vg.prog <- as.data.frame(vg.tumor[,"Progenitor-like"])
vg.prog$Marker <- "Progenitor"
colnames(vg.prog)[1] <- "Gene"

# prepare monocytic markers
vg.mono <- as.data.frame(vg.tumor[,"Monocyte-like"])
vg.mono$Marker <- "Monocyte"
colnames(vg.mono)[1] <- "Gene"

# combine prog and mono markers
vg.prog.mono <- rbind(vg.prog, vg.mono)

### get global diffexp for each exp: NRAS ASO or KO vs. Neg ASO or KO
# global.inputs <- list("Patient ASO" = "syn64732753",
#                "Cell Line ASO" = "syn64734615",
#                "Cell Line KO" = "syn64735095",
#                "Cell Line ASO 2HR" = "syn64736130",
#                "Cell Line ASO 6HR" = "syn64738709",
#                "Cell Line ASO 24HR" = "syn64739111")
# phospho.inputs <- list("Patient ASO" = "syn64732862",
#                        "Cell Line ASO" = "syn64734712",
#                        "Cell Line KO" = "syn64735201",
#                        "Cell Line ASO 2HR" = "syn64736236",
#                        "Cell Line ASO 6HR" = "syn64738821",
#                        "Cell Line ASO 24HR" = "syn64739245")
global.inputs <- list("Patient ASO" = "syn66726954",
                      "Cell Line ASO" = "syn66726844",
                      "Cell Line KO" = "syn66727085",
                      "Cell Line ASO 2HR" = "syn66727770",
                      "Cell Line ASO 6HR" = "syn66729260",
                      "Cell Line ASO 24HR" = "syn66730823")
phospho.inputs <- list("Patient ASO" = "syn66726980",
                       "Cell Line ASO" = "syn66726866",
                       "Cell Line KO" = "syn66727197",
                       "Cell Line ASO 2HR" = "syn66727869",
                       "Cell Line ASO 6HR" = "syn66729468",
                       "Cell Line ASO 24HR" = "syn66730918")
all.inputs <- list("Global" = global.inputs, "Phospho" = phospho.inputs)
library(ggplot2)
for (omics in names(all.inputs)) {
  inputs <- all.inputs[[omics]]
  setwd("~/OneDrive - PNNL/Documents/GitHub/Exp21_NRAS_ASO_treated_patients/proteomics/")
  dir.create("analysis")
  setwd("analysis")
  dir.create("Monocyte_vs_progenitor_signatures_20250530")
  setwd("Monocyte_vs_progenitor_signatures_20250530")
  degs <- list()
  for (i in names(inputs)) {
    # read in global diff exp
    plot.df <- read.csv(synapser::synGet(inputs[[i]])$path)
    sig.plot.df <- na.omit(plot.df[plot.df$adj.P.Val <= 0.05,])
    if (nrow(sig.plot.df) > 0) {
      degs[[i]] <- sig.plot.df
    }
    
    if (omics != "Phospho") {
      # label for van Galen mono or prog markers
      plot.df$Marker <- NA
      if (any(plot.df$Gene %in% vg.prog$Gene)) {
        plot.df[plot.df$Gene %in% vg.prog$Gene,]$Marker <- "Progenitor-like"
        nProgMarkers <- length(which(plot.df$Gene %in% vg.prog$Gene))
        nSigProgMarkers <- nrow(plot.df[plot.df$Gene %in% vg.prog$Gene &
                                          plot.df$adj.P.Val <= 0.05,])
      }

      if (any(plot.df$Gene %in% vg.mono$Gene)) {
        plot.df[plot.df$Gene %in% vg.mono$Gene,]$Marker <- "Monocyte-like"
        nMonoMarkers <- length(which(plot.df$Gene %in% vg.mono$Gene))
        nSigMonoMarkers <- nrow(plot.df[plot.df$Gene %in% vg.mono$Gene &
                                          plot.df$adj.P.Val <= 0.05,])
      }
      nMarkers <- nrow(plot.df[!is.na(plot.df$Marker),])
      nSigMarkers <- nrow(plot.df[!is.na(plot.df$Marker) & plot.df$adj.P.Val <= 0.05,])
      # temp.title <- paste0(nSigMonoMarkers, " / ", nMonoMarkers, " monocyte-like genes and ",
      #                      nSigProgMarkers, " / ", nProgMarkers, " progenitor-like genes ",
      #                      "differentially expressed (adjusted p <= 0.05)")
      temp.title <- "NRAS knockdown effect on AML cell state, from proteomics"
      
      # # plot all features
      # limit.x <- ceiling(max(abs(plot.df$Log2FC)))
      # limit.y <- ceiling(max(abs(-log10(plot.df$P.Value))))
      # title.x <- "Log2 Fold Change (NRAS ASO / CTRL ASO)"
      # if (grepl(" KO",i)) {
      #   title.x <- "Log2 Fold Change (NRAS single guide / non-target guide)"
      # }
      # ggplot2::ggplot(plot.df, aes(x=Log2FC, y = -log10(P.Value), color = Marker)) +
      #   ggplot2::geom_point(size = 4) + theme_classic(base_size = 12) +
      #   ggrepel::geom_text_repel(
      #     data = subset(plot.df, !is.na(Marker)),
      #     mapping = aes(label = Gene, size = I(4)), nudge_y = 0.25
      #   ) + ggtitle(temp.title) +
      #   ggplot2::scale_color_manual(
      #     values = c(scales::hue_pal()(2), "grey"), name = "van Galen gene set",
      #     breaks = c("Progenitor-like", "Monocyte-like", "Other")
      #   ) +
      #   ggplot2::xlim(-limit.x, limit.x) +
      #   ggplot2::ylim(0, limit.y) +
      #   ggplot2::xlab(title.x) +
      #   ggplot2::ylab("-Log10 P-value") +
      #   ggplot2::geom_vline(xintercept = 0, linetype = "solid", color = "black",
      #                       linewidth = 0.5) +
      #   ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black",
      #                       linewidth = 0.5)
      # ggplot2::ggsave(paste0(i, "_vanGalenMarkers_diffexp_volcanoPlot.pdf"), width=5, height=3)
      
      # plot only markers
      plot.df <- na.omit(plot.df)
      limit.x <- ceiling(max(abs(plot.df$Log2FC)))
      limit.y <- ceiling(max(abs(-log10(plot.df$P.Value))))
      title.x <- expression(paste("Log"[2]," Fold Change (NRAS ASO / CTRL ASO)"))
      temp.sub <- "n=9 Late Gilterinib Resistant Molm14 samples per condition"
      if (grepl(" KO",i)) {
        title.x <- expression(paste("Log"[2]," Fold Change (NRAS CRISPR KO / CTRL)"))
        temp.sub <- "n=9 NRAS CRISPR KO and n=3 CTRL Late Gilterinib Resistant Molm14 samples"
      } else if (grepl("Patient",i)) {
        temp.sub <- "n=11 Gilterinib Resistant patient samples per condition"
      }
      ggplot2::ggplot(plot.df, aes(x=Log2FC, y = -log10(P.Value) #, label=Gene)
                                   ), guide="legend") +
        ggplot2::geom_point(aes(fill = Marker), color="black", pch=21, size = 3
          ) + theme_classic(base_size = 16 # was 20
            ) +
        # ggrepel::geom_label_repel(
        #   data = subset(plot.df, adj.P.Val <= 0.05),
        #   mapping = aes(label = Gene#, size = I(4)
        #                 ), nudge_y = 0.25
        # ) + 
        ggplot2::geom_vline(xintercept = 0, linetype = "solid", color = "black",
                            linewidth = 1) +
        ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black",
                            linewidth = 1) + theme(axis.line=element_line(size=1), plot.title=element_text(size=22)
                                                   ) + # was title size 27
        ggrepel::geom_label_repel(data = subset(plot.df, adj.P.Val <= 0.05),
                                  aes(fill=factor(Marker), label=Gene), color="black", show.legend=FALSE) + ggtitle(label=temp.title, subtitle=temp.sub) +
        ggplot2::scale_fill_manual(
          #values = c("orangered", "turquoise3", "grey"), 
          values = c("#8DD0BB", "#A158A3", "grey"), 
          name = "van Galen gene set",
          breaks = c("Progenitor-like", "Monocyte-like", "Other")
        ) +
        ggplot2::xlim(-limit.x, limit.x) +
        ggplot2::ylim(0, limit.y) + theme(legend.position="inside", legend.position.inside=c(0.8,0.8),
                                          #axis.text = element_text(size=22), 
                                          #axis.title = element_text(size=20),
                                          #legend.text = element_text(size=15),
                                          #legend.title = element_text(size=15), 
                                          plot.subtitle=element_text(face="italic")
                                          ) +
        labs(x=title.x, y = expression(paste("-Log"[10]," P-value")), label=element_blank()) #+
         #+
        #geom_point(data = plot.df, col = "black", stroke = 1, shape = 21)
      ggplot2::ggsave(paste0(omics, "_", i, "_vanGalenMarkersOnly_diffexp_volcanoPlot_labels_20250625.pdf"), width=8.5, height=7.5) 
    }
  }
  
  # compare exps
  venn.list <- list()
  for (i in c("Patient ASO", "Cell Line ASO", "Cell Line KO")) {
    if(omics == "Phospho") {
      venn.list[[i]] <- degs[[i]]$SUB_SITE
    } else {
      venn.list[[i]] <- degs[[i]]$Gene
    }
  }
  ggvenn::ggvenn(venn.list, show_percentage = FALSE, text_size = 5, set_name_size=5)
  ggplot2::ggsave(paste0(omics, "_NRAS_DEG_vennDiagram_",Sys.Date(),".pdf"), width=5, height=5)
  
  ## run correlations
  corr.df <- as.data.frame(tidyr::expand_grid(x=names(degs), y=names(degs)))
  corr.df[,c("Pearson.est", "Pearson.p", "Spearman.est", "Spearman.p", "N")] <- NA
  DEG.df <- na.omit(data.table::rbindlist(degs, use.names = TRUE, idcol = "type"))
  if (omics == "Phospho") {
    Log2FC.df <- reshape2::dcast(DEG.df, SUB_SITE ~ type,value.var = "Log2FC", fill = NA)
    rownames(Log2FC.df) <- Log2FC.df$SUB_SITE
  } else {
    Log2FC.df <- reshape2::dcast(DEG.df, Gene ~ type,value.var = "Log2FC", fill = NA)
    rownames(Log2FC.df) <- Log2FC.df$Gene 
  }
  Log2FC.mat <- as.matrix(Log2FC.df[, 2:ncol(Log2FC.df)])
  for (i in names(degs)) {
    for (j in names(degs)) {
      temp.input <- na.omit(Log2FC.mat[,c(i, j)])
      if (nrow(temp.input) > 2) {
        temp.corr <- cor.test(temp.input[,1], temp.input[,2])
        corr.df[corr.df$x == i & corr.df$y == j,]$Pearson.p <- temp.corr$p.value
        corr.df[corr.df$x == i & corr.df$y == j,]$Pearson.est <- temp.corr$estimate
        temp.rho <- cor.test(temp.input[,1], temp.input[,2], method="spearman")
        corr.df[corr.df$x == i & corr.df$y == j,]$Spearman.p <- temp.rho$p.value
        corr.df[corr.df$x == i & corr.df$y == j,]$Spearman.est <- temp.rho$estimate
        corr.df[corr.df$x == i & corr.df$y == j,]$N <- nrow(temp.input)
        
        if (i != j) {
          stats_pearson <- substitute(
            r == est * "," ~ ~"p" ~ "=" ~ p,
            list(
              est = as.numeric(format(temp.corr$estimate, digits = 3)),
              p = format(temp.corr$p.value, digits = 3)
            )
          )
          colnames(temp.input) <- c("x", "y")
          ggplot2::ggplot(temp.input, aes(x=x, y=y)) + geom_point() + 
            theme_minimal(base_size=12) + ggplot2::geom_smooth(method = "lm", size = 1.5, 
                                                               linetype = "solid", color = "blue", se = TRUE, na.rm = TRUE) +
            ggplot2::geom_text(x = max(temp.input[,"x"]), y = min(temp.input[,"y"]),
                               vjust = "inward", hjust = "inward", colour = "blue", parse = TRUE,
                               label = as.character(as.expression(stats_pearson)), size = 8) + 
            ggtitle(paste0("Log2 Fold-change (N = ",nrow(temp.input),")")) +
            geom_hline(yintercept=0) + geom_vline(xintercept=0) + 
            labs(x=i, y=j, title = paste0("Log2 Fold-change (N = ",nrow(temp.input), ")"))
          ggplot2::ggsave(paste0(omics,"_",i, "_corrWith_", j, "_NRAS_diffexp_Log2FC.pdf"), width=5, height=5) 
          
          stats_spearman <- substitute(
            rho == est * "," ~ ~"p" ~ "=" ~ p,
            list(
              est = as.numeric(format(temp.rho$estimate, digits = 3)),
              p = format(temp.rho$p.value, digits = 3)
            )
          )
          colnames(temp.input) <- c("x", "y")
          ggplot2::ggplot(temp.input, aes(x=rank(x), y=rank(y))) + geom_point() + 
            theme_minimal(base_size=12) + ggplot2::geom_smooth(method = "lm", size = 1.5, 
                                                               linetype = "solid", color = "blue", se = TRUE, na.rm = TRUE) +
            ggplot2::geom_text(x = max(temp.input[,"x"]), y = min(temp.input[,"y"]),
                               vjust = "inward", hjust = "inward", colour = "blue", parse = TRUE,
                               label = as.character(as.expression(stats_spearman)), size = 8) + 
            ggtitle(paste0("Log2 Fold-change (N = ",nrow(temp.input),")")) +
            geom_hline(yintercept=0) + geom_vline(xintercept=0) +
            labs(x=i, y=j, title = paste0("Log2 Fold-change Rank (N = ",nrow(temp.input), ")"))
          ggplot2::ggsave(paste0(omics,"_",i, "_spearmanCorrWith_", j, "_NRAS_diffexp_Log2FC.pdf"), width=5, height=5) 
          
        }
      }
    }
  }
  write.csv(corr.df, paste0(omics, "_NRAS_diffexp_correlations_",Sys.Date(),".csv"), row.names=FALSE)
  corr.df$minusLogP <- NA
  corr.df[corr.df$Pearson.p != 0,]$minusLogP <- -log10(corr.df[corr.df$Pearson.p != 0,]$Pearson.p)
  
  corr.df2 <- corr.df
  corr.df2$x <- factor(corr.df2$x, levels=names(degs))
  corr.df2[,1] <- as.integer(corr.df2[,1])
  corr.df2$y <- factor(corr.df2$y, levels=names(degs))
  corr.df2[,2] <- as.integer(corr.df2[,2])
  corr.mat <- Matrix::sparseMatrix(i = corr.df2[,1], 
                                   j = corr.df2[,2], 
                                   x = corr.df2[,3])
  corr.mat <- as.matrix(corr.mat)
  colnames(corr.mat) <- names(degs)
  rownames(corr.mat) <- names(degs)
  corr.mat.plot <- ggcorrplot::ggcorrplot(corr.mat)
  ggplot2::ggsave(paste0(omics, "_Gilt_w_NRAS_ctl_deg_Log2FC_correlations_", Sys.Date(), ".pdf"), corr.mat.plot, width=5, height=5)

  corr.mat <- Matrix::sparseMatrix(i = corr.df2[,1], 
                                   j = corr.df2[,2], 
                                   x = corr.df2[,5])
  corr.mat <- as.matrix(corr.mat)
  colnames(corr.mat) <- names(degs)
  rownames(corr.mat) <- names(degs)
  corr.mat <- corr.mat[colSums(corr.mat) != 0, rowSums(corr.mat) != 0]
  corr.mat.plot <- ggcorrplot::ggcorrplot(corr.mat)
  ggplot2::ggsave(paste0(omics,"_Gilt_w_NRAS_ctl_deg_Log2FC_spearmanCorrelations_", Sys.Date(), ".pdf"), corr.mat.plot, width=5, height=5)
}

for (omics in names(all.inputs)) {
  inputs <- all.inputs[[omics]]
  setwd("~/OneDrive - PNNL/Documents/GitHub/Exp21_NRAS_ASO_treated_patients/proteomics/")
  dir.create("analysis")
  setwd("analysis")
  dir.create("Monocyte_vs_progenitor_signatures_20250530")
  setwd("Monocyte_vs_progenitor_signatures_20250530")
  degs <- list()
  for (i in names(inputs)) {
    # read in global diff exp
    plot.df <- read.csv(synapser::synGet(inputs[[i]])$path)
    sig.plot.df <- na.omit(plot.df[plot.df$adj.P.Val <= 0.001,])
    if (nrow(sig.plot.df) > 0) {
      degs[[i]] <- sig.plot.df
    }
    
    if (omics != "Phospho") {
      # label for van Galen mono or prog markers
      plot.df$Marker <- NA
      if (any(plot.df$Gene %in% vg.prog$Gene)) {
        plot.df[plot.df$Gene %in% vg.prog$Gene,]$Marker <- "Progenitor-like"
        nProgMarkers <- length(which(plot.df$Gene %in% vg.prog$Gene))
        nSigProgMarkers <- nrow(plot.df[plot.df$Gene %in% vg.prog$Gene &
                                          plot.df$adj.P.Val <= 0.001,])
      }
      if (any(plot.df$Gene %in% vg.mono$Gene)) {
        plot.df[plot.df$Gene %in% vg.mono$Gene,]$Marker <- "Monocyte-like"
        nMonoMarkers <- length(which(plot.df$Gene %in% vg.mono$Gene))
        nSigMonoMarkers <- nrow(plot.df[plot.df$Gene %in% vg.mono$Gene &
                                          plot.df$adj.P.Val <= 0.001,])
      }
      nMarkers <- nrow(plot.df[!is.na(plot.df$Marker),])
      nSigMarkers <- nrow(plot.df[!is.na(plot.df$Marker) & plot.df$adj.P.Val <= 0.001,])
      temp.title <- paste0(nSigMonoMarkers, " / ", nMonoMarkers, " monocyte-like genes and ",
                           nSigProgMarkers, " / ", nProgMarkers, " progenitor-like genes ",
                           "differentially expressed (adjusted p <= 0.001)")
      
      # # plot all features
      # limit.x <- ceiling(max(abs(plot.df$Log2FC)))
      # limit.y <- ceiling(max(abs(-log10(plot.df$P.Value))))
      # title.x <- "Log2 Fold Change (NRAS ASO / CTRL ASO)"
      # if (grepl(" KO",i)) {
      #   title.x <- "Log2 Fold Change (NRAS single guide / non-target guide)"
      # }
      # ggplot2::ggplot(plot.df, aes(x=Log2FC, y = -log10(P.Value), color = Marker)) +
      #   ggplot2::geom_point(size = 4) + theme_classic(base_size = 12) +
      #   ggrepel::geom_text_repel(
      #     data = subset(plot.df, !is.na(Marker)),
      #     mapping = aes(label = Gene, size = I(4)), nudge_y = 0.25
      #   ) + ggtitle(temp.title) +
      #   ggplot2::scale_color_manual(
      #     values = c(scales::hue_pal()(2), "grey"), name = "van Galen gene set",
      #     breaks = c("Progenitor-like", "Monocyte-like", "Other")
      #   ) +
      #   ggplot2::xlim(-limit.x, limit.x) +
      #   ggplot2::ylim(0, limit.y) +
      #   ggplot2::xlab(title.x) +
      #   ggplot2::ylab("-Log10 P-value") +
      #   ggplot2::geom_vline(xintercept = 0, linetype = "solid", color = "black",
      #                       linewidth = 0.5) +
      #   ggplot2::geom_hline(yintercept = -log10(0.001), linetype = "dashed", color = "black",
      #                       linewidth = 0.5)
      # ggplot2::ggsave(paste0(i, "_vanGalenMarkers_diffexp_q0.001_volcanoPlot.pdf"), width=5, height=3)
      
      # plot only markers
      plot.df <- na.omit(plot.df)
      limit.x <- ceiling(max(abs(plot.df$Log2FC)))
      limit.y <- ceiling(max(abs(-log10(plot.df$P.Value))))
      title.x <- "Log2 Fold Change (NRAS ASO / CTRL ASO)"
      if (grepl(" KO",i)) {
        title.x <- "Log2 Fold Change (NRAS single guide / non-target guide)"
      }
      ggplot2::ggplot(plot.df, aes(x=Log2FC, y = -log10(P.Value), color = Marker)) +
        ggplot2::geom_point(size = 4) + theme_classic(base_size = 12) +
        ggrepel::geom_text_repel(
          data = subset(plot.df, adj.P.Val <= 0.001),
          mapping = aes(label = Gene, size = I(4)), nudge_y = 0.25
        ) + ggtitle(temp.title) +
        ggplot2::scale_color_manual(
          values = c(scales::hue_pal()(2), "grey"), name = "van Galen gene set",
          breaks = c("Progenitor-like", "Monocyte-like", "Other")
        ) +
        ggplot2::xlim(-limit.x, limit.x) +
        ggplot2::ylim(0, limit.y) +
        ggplot2::xlab(title.x) +
        ggplot2::ylab("-Log10 P-value") +
        ggplot2::geom_vline(xintercept = 0, linetype = "solid", color = "black",
                            linewidth = 0.5) +
        ggplot2::geom_hline(yintercept = -log10(0.001), linetype = "dashed", color = "black",
                            linewidth = 0.5)
      ggplot2::ggsave(paste0(omics, "_", i, "_vanGalenMarkersOnly_diffexp_q0.001_volcanoPlot.pdf"), width=10, height=7) 
    }
  }
  
  # compare exps
  venn.list <- list()
  for (i in c("Patient ASO", "Cell Line ASO", "Cell Line KO")) {
    if(omics == "Phospho") {
      venn.list[[i]] <- degs[[i]]$SUB_SITE
    } else {
      venn.list[[i]] <- degs[[i]]$Gene
    }
  }
  ggvenn::ggvenn(venn.list, show_percentage = FALSE, text_size = 5, set_name_size=5)
  ggplot2::ggsave(paste0(omics, "_NRAS_DEG_vennDiagram_q0.001_",Sys.Date(),".pdf"), width=5, height=5)
  
  ## run correlations
  corr.df <- as.data.frame(tidyr::expand_grid(x=names(degs), y=names(degs)))
  corr.df[,c("Pearson.est", "Pearson.p", "Spearman.est", "Spearman.p", "N")] <- NA
  DEG.df <- na.omit(data.table::rbindlist(degs, use.names = TRUE, idcol = "type"))
  if (omics == "Phospho") {
    Log2FC.df <- reshape2::dcast(DEG.df, SUB_SITE ~ type,value.var = "Log2FC", fill = NA)
    rownames(Log2FC.df) <- Log2FC.df$SUB_SITE
  } else {
    Log2FC.df <- reshape2::dcast(DEG.df, Gene ~ type,value.var = "Log2FC", fill = NA)
    rownames(Log2FC.df) <- Log2FC.df$Gene 
  }
  Log2FC.mat <- as.matrix(Log2FC.df[, 2:ncol(Log2FC.df)])
  for (i in names(degs)) {
    for (j in names(degs)) {
      temp.input <- na.omit(Log2FC.mat[,c(i, j)])
      if (nrow(temp.input) > 2) {
        temp.corr <- cor.test(temp.input[,1], temp.input[,2])
        corr.df[corr.df$x == i & corr.df$y == j,]$Pearson.p <- temp.corr$p.value
        corr.df[corr.df$x == i & corr.df$y == j,]$Pearson.est <- temp.corr$estimate
        temp.rho <- cor.test(temp.input[,1], temp.input[,2], method="spearman")
        corr.df[corr.df$x == i & corr.df$y == j,]$Spearman.p <- temp.rho$p.value
        corr.df[corr.df$x == i & corr.df$y == j,]$Spearman.est <- temp.rho$estimate
        corr.df[corr.df$x == i & corr.df$y == j,]$N <- nrow(temp.input)
        
        if (i != j) {
          stats_pearson <- substitute(
            r == est * "," ~ ~"p" ~ "=" ~ p,
            list(
              est = as.numeric(format(temp.corr$estimate, digits = 3)),
              p = format(temp.corr$p.value, digits = 3)
            )
          )
          colnames(temp.input) <- c("x", "y")
          ggplot2::ggplot(temp.input, aes(x=x, y=y)) + geom_point() + 
            theme_minimal(base_size=12) + ggplot2::geom_smooth(method = "lm", size = 1.5, 
                                                               linetype = "solid", color = "blue", se = TRUE, na.rm = TRUE) +
            ggplot2::geom_text(x = max(temp.input[,"x"]), y = min(temp.input[,"y"]),
                               vjust = "inward", hjust = "inward", colour = "blue", parse = TRUE,
                               label = as.character(as.expression(stats_pearson)), size = 8) + 
            ggtitle(paste0("Log2 Fold-change (N = ",nrow(temp.input),")")) +
            geom_hline(yintercept=0) + geom_vline(xintercept=0) + 
            labs(x=i, y=j, title = paste0("Log2 Fold-change (N = ",nrow(temp.input), ")"))
          ggplot2::ggsave(paste0(omics,"_",i, "_corrWith_", j, "_NRAS_diffexp_q0.001_Log2FC.pdf"), width=5, height=5) 
          
          stats_spearman <- substitute(
            rho == est * "," ~ ~"p" ~ "=" ~ p,
            list(
              est = as.numeric(format(temp.rho$estimate, digits = 3)),
              p = format(temp.rho$p.value, digits = 3)
            )
          )
          colnames(temp.input) <- c("x", "y")
          ggplot2::ggplot(temp.input, aes(x=rank(x), y=rank(y))) + geom_point() + 
            theme_minimal(base_size=12) + ggplot2::geom_smooth(method = "lm", size = 1.5, 
                                                               linetype = "solid", color = "blue", se = TRUE, na.rm = TRUE) +
            ggplot2::geom_text(x = max(temp.input[,"x"]), y = min(temp.input[,"y"]),
                               vjust = "inward", hjust = "inward", colour = "blue", parse = TRUE,
                               label = as.character(as.expression(stats_spearman)), size = 8) + 
            ggtitle(paste0("Log2 Fold-change (N = ",nrow(temp.input),")")) +
            geom_hline(yintercept=0) + geom_vline(xintercept=0) +
            labs(x=i, y=j, title = paste0("Log2 Fold-change Rank (N = ",nrow(temp.input), ")"))
          ggplot2::ggsave(paste0(omics,"_",i, "_spearmanCorrWith_", j, "_NRAS_diffexp_q0.001_Log2FC.pdf"), width=5, height=5) 
          
        }
      }
    }
  }
  write.csv(corr.df, paste0(omics, "_NRAS_diffexp_correlations_q0.001_",Sys.Date(),".csv"), row.names=FALSE)
  corr.df$minusLogP <- NA
  corr.df[corr.df$Pearson.p != 0,]$minusLogP <- -log10(corr.df[corr.df$Pearson.p != 0,]$Pearson.p)
  
  corr.df2 <- corr.df
  corr.df2$x <- factor(corr.df2$x, levels=names(degs))
  corr.df2[,1] <- as.integer(corr.df2[,1])
  corr.df2$y <- factor(corr.df2$y, levels=names(degs))
  corr.df2[,2] <- as.integer(corr.df2[,2])
  corr.mat <- Matrix::sparseMatrix(i = corr.df2[,1], 
                                   j = corr.df2[,2], 
                                   x = corr.df2[,3])
  corr.mat <- as.matrix(corr.mat)
  colnames(corr.mat) <- names(degs)
  rownames(corr.mat) <- names(degs)
  corr.mat.plot <- ggcorrplot::ggcorrplot(corr.mat)
  ggplot2::ggsave(paste0(omics, "_Gilt_w_NRAS_ctl_deg_Log2FC_correlations_", Sys.Date(), ".pdf"), corr.mat.plot, width=5, height=5)
  
  corr.mat <- Matrix::sparseMatrix(i = corr.df2[,1], 
                                   j = corr.df2[,2], 
                                   x = corr.df2[,5])
  corr.mat <- as.matrix(corr.mat)
  colnames(corr.mat) <- names(degs)
  rownames(corr.mat) <- names(degs)
  corr.mat <- corr.mat[colSums(corr.mat) != 0, rowSums(corr.mat) != 0]
  corr.mat.plot <- ggcorrplot::ggcorrplot(corr.mat)
  ggplot2::ggsave(paste0(omics,"_Gilt_w_NRAS_ctl_deg_Log2FC_spearmanCorrelations_", Sys.Date(), ".pdf"), corr.mat.plot, width=5, height=5)
}
# `geom_smooth()` using formula = 'y ~ x'
# `geom_smooth()` using formula = 'y ~ x'
# `geom_smooth()` using formula = 'y ~ x'
# `geom_smooth()` using formula = 'y ~ x'
# Error in `[<-.data.frame`(`*tmp*`, corr.df$Pearson.p != 0, , value = list( : 
#                                                                              missing values are not allowed in subscripted assignments of data frames

# #### repeat for NRAS + Gilt vs. NRAS ####
# global.inputs <- list("Patient ASO" = "syn64742147", # none sig
#                "Cell Line ASO" = "syn64742334",
#                "Cell Line ASO 2HR" = "syn64741625",
#                "Cell Line ASO 6HR" = "syn64741801",
#                "Cell Line ASO 24HR" = "syn64741981")
# phospho.inputs <- list("Patient ASO" = "syn64742251", # none sig
#                       "Cell Line ASO" = "syn64742434",
#                       "Cell Line ASO 2HR" = "syn64741730",
#                       "Cell Line ASO 6HR" = "syn64741904",
#                       "Cell Line ASO 24HR" = "syn64742091")
# all.inputs <- list("Global" = global.inputs, "Phospho" = phospho.inputs)
# library(ggplot2)
# for (omics in names(all.inputs)) {
#   inputs <- all.inputs[[omics]]
#   setwd("~/OneDrive - PNNL/Documents/GitHub/Exp21_NRAS_ASO_treated_patients/proteomics/")
#   dir.create("analysis")
#   setwd("analysis")
#   dir.create("Monocyte_vs_progenitor_signatures_20250530")
#   setwd("Monocyte_vs_progenitor_signatures_20250530")
#   degs <- list()
#   for (i in names(inputs)) {
#     # read in global diff exp
#     plot.df <- read.csv(synapser::synGet(inputs[[i]])$path)
#     sig.plot.df <- na.omit(plot.df[plot.df$adj.P.Val <= 0.05,])
#     if (nrow(sig.plot.df) > 0) {
#       degs[[i]] <- sig.plot.df
#     }
#   }
#   library(ggplot2)
#   #compiled <- compile_mDEG(degs)
#   
#   # compare exps
#   venn.list <- list()
#   for (i in names(degs)) {
#     if(omics == "Phospho") {
#       venn.list[[i]] <- degs[[i]]$SUB_SITE
#     } else {
#       venn.list[[i]] <- degs[[i]]$Gene
#     }
#   }
#   ggvenn::ggvenn(venn.list, show_percentage = FALSE, text_size = 5, set_name_size=5)
#   ggplot2::ggsave(paste0(omics,"_Gilt_w_NRAS_ctl_DEG_vennDiagram_",Sys.Date(),".pdf"), width=5, height=5)
#   
#   ## run correlations
#   corr.df <- as.data.frame(tidyr::expand_grid(x=names(degs), y=names(degs)))
#   corr.df[,c("Pearson.est", "Pearson.p", "Spearman.est", "Spearman.p", "N")] <- NA
#   DEG.df <- na.omit(data.table::rbindlist(degs, use.names = TRUE, idcol = "type"))
#   if (omics == "Phospho") {
#     Log2FC.df <- reshape2::dcast(DEG.df, SUB_SITE ~ type,value.var = "Log2FC", fill = NA)
#     rownames(Log2FC.df) <- Log2FC.df$SUB_SITE
#   } else {
#     Log2FC.df <- reshape2::dcast(DEG.df, Gene ~ type,value.var = "Log2FC", fill = NA)
#     rownames(Log2FC.df) <- Log2FC.df$Gene 
#   }
#   Log2FC.mat <- as.matrix(Log2FC.df[, 2:ncol(Log2FC.df)])
#   for (i in names(degs)) {
#     for (j in names(degs)) {
#       temp.input <- na.omit(Log2FC.mat[,c(i, j)])
#       if (nrow(temp.input) > 2) {
#         temp.corr <- cor.test(temp.input[,1], temp.input[,2])
#         corr.df[corr.df$x == i & corr.df$y == j,]$Pearson.p <- temp.corr$p.value
#         corr.df[corr.df$x == i & corr.df$y == j,]$Pearson.est <- temp.corr$estimate
#         temp.rho <- cor.test(temp.input[,1], temp.input[,2], method="spearman")
#         corr.df[corr.df$x == i & corr.df$y == j,]$Spearman.p <- temp.rho$p.value
#         corr.df[corr.df$x == i & corr.df$y == j,]$Spearman.est <- temp.rho$estimate
#         corr.df[corr.df$x == i & corr.df$y == j,]$N <- nrow(temp.input)
#         
#         if (i != j) {
#           stats_pearson <- substitute(
#             r == est * "," ~ ~"p" ~ "=" ~ p,
#             list(
#               est = as.numeric(format(temp.corr$estimate, digits = 3)),
#               p = format(temp.corr$p.value, digits = 3)
#             )
#           )
#           colnames(temp.input) <- c("x", "y")
#           ggplot2::ggplot(temp.input, aes(x=x, y=y)) + geom_point() + 
#             theme_minimal(base_size=12) + ggplot2::geom_smooth(method = "lm", size = 1.5, 
#                                                                linetype = "solid", color = "blue", se = TRUE, na.rm = TRUE) +
#             ggplot2::geom_text(x = max(temp.input[,"x"]), y = min(temp.input[,"y"]),
#                                vjust = "inward", hjust = "inward", colour = "blue", parse = TRUE,
#                                label = as.character(as.expression(stats_pearson)), size = 8) + 
#             ggtitle(paste0("Log2 Fold-change (N = ",nrow(temp.input),")")) +
#             geom_hline(yintercept=0) + geom_vline(xintercept=0) + 
#             labs(x=i, y=j, title = paste0("Log2 Fold-change (N = ",nrow(temp.input), ")"))
#           ggplot2::ggsave(paste0(omics,"_",i, "_corrWith_", j, "_Gilt_w_NRAS_ctl_diffexp_Log2FC.pdf"), width=5, height=5) 
#           
#           stats_spearman <- substitute(
#             rho == est * "," ~ ~"p" ~ "=" ~ p,
#             list(
#               est = as.numeric(format(temp.rho$estimate, digits = 3)),
#               p = format(temp.rho$p.value, digits = 3)
#             )
#           )
#           colnames(temp.input) <- c("x", "y")
#           ggplot2::ggplot(temp.input, aes(x=rank(x), y=rank(y))) + geom_point() + 
#             theme_minimal(base_size=12) + ggplot2::geom_smooth(method = "lm", size = 1.5, 
#                                                                linetype = "solid", color = "blue", se = TRUE, na.rm = TRUE) +
#             ggplot2::geom_text(x = max(temp.input[,"x"]), y = min(temp.input[,"y"]),
#                                vjust = "inward", hjust = "inward", colour = "blue", parse = TRUE,
#                                label = as.character(as.expression(stats_spearman)), size = 8) + 
#             ggtitle(paste0("Log2 Fold-change (N = ",nrow(temp.input),")")) +
#             geom_hline(yintercept=0) + geom_vline(xintercept=0) +
#             labs(x=i, y=j, title = paste0("Log2 Fold-change Rank (N = ",nrow(temp.input), ")"))
#           ggplot2::ggsave(paste0(omics,"_",i, "_spearmanCorrWith_", j, "_Gilt_w_NRAS_ctl_diffexp_Log2FC.pdf"), width=5, height=5) 
#           
#         }
#       }
#     }
#   }
#   write.csv(corr.df, paste0(omics,"_Gilt_w_NRAS_ctl_diffexp_correlations_",Sys.Date(),".csv"), row.names=FALSE)
#   corr.df <- na.omit(corr.df)
#   corr.df$minusLogP <- NA
#   corr.df[corr.df$Pearson.p != 0,]$minusLogP <- -log10(corr.df[corr.df$Pearson.p != 0,]$Pearson.p)
#   
#   corr.df2 <- corr.df
#   corr.df2$x <- factor(corr.df2$x, levels=names(degs))
#   corr.df2[,1] <- as.integer(corr.df2[,1])
#   corr.df2$y <- factor(corr.df2$y, levels=names(degs))
#   corr.df2[,2] <- as.integer(corr.df2[,2])
#   corr.mat <- Matrix::sparseMatrix(i = corr.df2[,1], 
#                                    j = corr.df2[,2], 
#                                    x = corr.df2[,3])
#   corr.mat <- as.matrix(corr.mat)
#   colnames(corr.mat) <- names(degs)
#   rownames(corr.mat) <- names(degs)
#   corr.mat <- corr.mat[colSums(corr.mat) != 0, rowSums(corr.mat) != 0]
#   corr.mat.plot <- ggcorrplot::ggcorrplot(corr.mat)
#   ggplot2::ggsave(paste0(omics,"_Gilt_w_NRAS_ctl_deg_Log2FC_correlations_", Sys.Date(), ".pdf"), corr.mat.plot, width=5, height=5)
#   
#   corr.mat <- Matrix::sparseMatrix(i = corr.df2[,1], 
#                                    j = corr.df2[,2], 
#                                    x = corr.df2[,5])
#   corr.mat <- as.matrix(corr.mat)
#   colnames(corr.mat) <- names(degs)
#   rownames(corr.mat) <- names(degs)
#   corr.mat <- corr.mat[colSums(corr.mat) != 0, rowSums(corr.mat) != 0]
#   corr.mat.plot <- ggcorrplot::ggcorrplot(corr.mat)
#   ggplot2::ggsave(paste0(omics,"_Gilt_w_NRAS_ctl_deg_Log2FC_spearmanCorrelations_", Sys.Date(), ".pdf"), corr.mat.plot, width=5, height=5)
# }

#### run deconvolution for mono vs. prog ####
setwd("~/OneDrive - PNNL/Documents/GitHub/Exp21_NRAS_ASO_treated_patients/proteomics/")
dir.create("analysis")
setwd("analysis")
base.path <- "~/OneDrive - PNNL/Documents/GitHub/Exp21_NRAS_ASO_treated_patients/proteomics/analysis/"

# pt <- list("Meta" = "syn64732333", # exp 21
#            "Global" = "syn52751324",
#            "Phospho" = "syn52751606")
# cell <- list("Meta" = "syn52653527", # exp 20
#              "Global" = "syn52623904",
#              "Phospho" = "syn52623908")
# cellKO <- list("Meta" = "syn64731726", # exp 22
#                "Global" = "syn53212316",
#                "Phospho" = "syn53212321")
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
#synIDs <- c("syn53180709","syn52256487","syn53212344")
synIDs <- c("syn66726942","syn66726775","syn66705127")
names(synIDs) <- names(inputs)
exp.p.df <- data.frame()
library(ggplot2)
for (i in names(inputs)) {
  cat("Running",i,"\n")
  
  ### read in data
  meta.df <- readxl::read_excel(synapser::synGet(inputs[[i]]$Meta)$path, 
                                sheet = meta.sheet[[i]])
  global.df <- read.table(synapser::synGet(inputs[[i]]$Global)$path, sep = "\t")
  
  # remove samples that Sunil said to exclude
  if (grepl("patient", i, ignore.case = TRUE)) {
    excluded.patients <- c("20-00576", "20-00432", "18-00101", "17-01083", "16-00815", "16-01253")
    meta.df <- meta.df[!(meta.df$PatientID %in% excluded.patients),] 
  }
  
  # add format data so that samples are along rows
  global.df <- as.data.frame(t(global.df))
  global.df$Sample <- rownames(global.df)
  global.df <- global.df[,c("Sample", colnames(global.df)[1:(ncol(global.df)-1)])]
  if (grepl("patient", i, ignore.case = TRUE)) {
    global.df <- global.df[!(global.df$Sample %in% excluded.patients),]
  }
  global.df100 <- global.df[,colSums(is.na(global.df)) == 0]
  
  setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Exp21_NRAS_ASO_treated_patients/proteomics/analysis/")
  dir.create("Monocyte_vs_progenitor_signatures_20250530")
  setwd("Monocyte_vs_progenitor_signatures_20250530")
  
  ### evaluate signature
  wv.df <- compareSigs(global.df100, sigs) 
  write.csv(wv.df, paste0(i,"_wv.csv"), row.names = FALSE)

  ### is monocytic score lower in treated vs. untreated samples?
  ## prep contrasts
  temp.meta <- list()
  if (grepl("patient", i, ignore.case = TRUE)) { # exp 21
    contrasts <- c("NRAS", "NRAS_w_Gilt_ctl", "Gilt", "Gilt_w_NRAS_ctl")
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
  } else if (grepl("ASO", i, ignore.case = TRUE)) { # exp 20
    contrasts <- c("NRAS", "NRAS_w_Gilt_ctl", "Gilt", "Gilt_w_NRAS_ctl")
    # NRAS ASO vs. NEG ASO
    aso.treatments <- c("NRAS ASO", "CTRL ASO")
    meta.df.nras <- as.data.frame(meta.df[meta.df$Group %in% aso.treatments,])
    meta.df.nras$NRAS <- "Untreated"
    meta.df.nras[grepl("NRAS ASO", meta.df.nras$Sample),]$NRAS <- "Treated"
    rownames(meta.df.nras) <- meta.df.nras$MeasurementName
    temp.meta[[contrasts[1]]] <- meta.df.nras
    
    # Gilt + NRAS vs. Gilt
    nras.gilt.treatments <- c("NRAS ASO + GILT", "CTRL ASO + GILT")
    meta.df.nras.gilt <- as.data.frame(meta.df[meta.df$Group %in% nras.gilt.treatments,])
    meta.df.nras.gilt$NRAS_w_Gilt_ctl <- "Untreated"
    meta.df.nras.gilt[meta.df.nras.gilt$Group == "NRAS ASO + GILT",]$NRAS_w_Gilt_ctl <- "Treated"
    rownames(meta.df.nras.gilt) <- meta.df.nras.gilt$MeasurementName
    temp.meta[[contrasts[2]]] <- meta.df.nras.gilt
    
    # Gilt + NRAS vs. Gilt
    gilt.treatments <- c("CTRL ASO + GILT", "CTRL ASO")
    meta.df.gilt <- as.data.frame(meta.df[meta.df$Group %in% gilt.treatments,])
    meta.df.gilt$Gilt <- "Untreated"
    meta.df.gilt[meta.df.gilt$Group == "CTRL ASO + GILT",]$Gilt <- "Treated"
    rownames(meta.df.gilt) <- meta.df.gilt$MeasurementName
    temp.meta[[contrasts[3]]] <- meta.df.gilt
    
    # Gilt + NRAS vs. Gilt
    nras.gilt.treatments <- c("NRAS ASO + GILT", "NRAS ASO")
    meta.df.nras.gilt <- as.data.frame(meta.df[meta.df$Group %in% nras.gilt.treatments,])
    meta.df.nras.gilt$Gilt_w_NRAS_ctl <- "Untreated"
    meta.df.nras.gilt[meta.df.nras.gilt$Group == "NRAS ASO + GILT",]$Gilt_w_NRAS_ctl <- "Treated"
    rownames(meta.df.nras.gilt) <- meta.df.nras.gilt$MeasurementName
    temp.meta[[contrasts[4]]] <- meta.df.nras.gilt
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
  }
  
  Signature <- c("Overall",na.omit(unique(wv.df$Signature)))
  contrast.p.df <- data.frame()
  ## for each contrast:
  for (j in contrasts) {
    setwd(base.path)
    # t-test
    p.df <- data.frame(Signature, p=NA, Paired=FALSE)
    treated.samples <- rownames(temp.meta[[j]][temp.meta[[j]][,j] == "Treated",])
    untreated.samples <- rownames(temp.meta[[j]][temp.meta[[j]][,j] == "Untreated",])
    paired <- all(na.omit(wv.df[wv.df$Sample %in% treated.samples,]$Sample) == na.omit(wv.df[wv.df$Sample %in% untreated.samples,]$Sample))
    temp.wv <- list("Treated" = na.omit(wv.df[wv.df$Sample %in% treated.samples,]$WV),
                    "Untreated" = na.omit(wv.df[wv.df$Sample %in% untreated.samples,]$WV))
    p.df[p.df$Signature=="Overall",]$p <- t.test(temp.wv$Treated, temp.wv$Untreated, "less", paired=paired)$p.value
    p.df[p.df$Signature=="Overall","Paired"] <- paired
    
    # repeat t-test for each signature
    other.sigs <- na.omit(unique(wv.df$Signature))
    for (k in other.sigs) {
      temp.wv.df <- wv.df[wv.df$Signature==k,]
      paired <- all(na.omit(temp.wv.df[temp.wv.df$Sample %in% treated.samples,]$Sample) == na.omit(temp.wv.df[temp.wv.df$Sample %in% untreated.samples,]$Sample))
      temp.wv <- list("Treated" = na.omit(temp.wv.df[temp.wv.df$Sample %in% treated.samples,]$WV),
                      "Untreated" = na.omit(temp.wv.df[temp.wv.df$Sample %in% untreated.samples,]$WV))
      p.df[p.df$Signature==k,]$p <- t.test(temp.wv$Treated, temp.wv$Untreated, "less", paired=paired)$p.value
      p.df[p.df$Signature==k,"Paired"] <- paired
    }
    
    # plot
    p.df$Significance <- "p > 0.05"
    p.df$p <- as.numeric(p.df$p)
    if (any(p.df$p <= 0.05)) {
      p.df[p.df$p <= 0.05,]$Significance <- "p <= 0.05" 
    }
    p.df$Significance <- factor(p.df$Significance,
                                levels=c("p <= 0.05", "p > 0.05"))
    p.df$minusLogP<- 1E-4
    if (any(p.df$p != 0)) {
      p.df[p.df$p != 0,]$minusLogP <- -log(p.df$p, base=10) 
    }
    sigOrder <- p.df[order(p.df$p),]$Signature
    ggplot(p.df, aes(x=Signature, y=minusLogP, 
                     fill = Signature, alpha = 0.5)) + 
      geom_col() + theme_classic(base_size = 12) + ylab("-Log(P-value)") + 
      ggplot2::scale_x_discrete(limits = sigOrder) +
      scale_fill_manual(values=RColorBrewer::brewer.pal(length(Signature), "Set2"), 
                        breaks=c("Sorted","Lasry","Triana","van Galen", "Overall"))+
      ggtitle(paste("T-test: Monocyte scores are lower in",j, "treated samples"))
    ggsave(paste0(i,"_", j, "_pValue_WV.pdf"), width = 5, height = 5)
    
    p.df$Contrast <- j
    contrast.p.df <- rbind(contrast.p.df, p.df)
  }
  contrast.p.df$Input <- i
  exp.p.df <- rbind(exp.p.df, contrast.p.df)
  
  ## repeat with timepoint groupings if exp 20
  if (i == "Cell Line ASO") {
    timepoints <- na.omit(unique(meta.df$Time[meta.df$Time != "baseline"]))
    for (k in timepoints) {
      temp.meta.df <- meta.df[meta.df$Time == k,]
      temp.meta <- list()
      contrasts <- c("NRAS", "NRAS_w_Gilt_ctl", "Gilt", "Gilt_w_NRAS_ctl")
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
      
      for (j in contrasts) {
        setwd(base.path)
        # t-test
        p.df <- data.frame(Signature, p=NA, Paired=FALSE)
        treated.samples <- rownames(temp.meta[[j]][temp.meta[[j]][,j] == "Treated",])
        untreated.samples <- rownames(temp.meta[[j]][temp.meta[[j]][,j] == "Untreated",])
        paired <- all(na.omit(wv.df[wv.df$Sample %in% treated.samples,]$Sample) == na.omit(wv.df[wv.df$Sample %in% untreated.samples,]$Sample))
        temp.wv <- list("Treated" = na.omit(wv.df[wv.df$Sample %in% treated.samples,]$WV),
                        "Untreated" = na.omit(wv.df[wv.df$Sample %in% untreated.samples,]$WV))
        p.df[p.df$Signature=="Overall",]$p <- t.test(temp.wv$Treated, temp.wv$Untreated, "less", paired=paired)$p.value
        p.df[p.df$Signature=="Overall","Paired"] <- paired
        
        # repeat t-test for each signature
        other.sigs <- na.omit(unique(wv.df$Signature))
        for (n in other.sigs) {
          temp.wv.df <- wv.df[wv.df$Signature==n,]
          paired <- all(na.omit(temp.wv.df[temp.wv.df$Sample %in% treated.samples,]$Sample) == na.omit(temp.wv.df[temp.wv.df$Sample %in% untreated.samples,]$Sample))
          temp.wv <- list("Treated" = na.omit(temp.wv.df[temp.wv.df$Sample %in% treated.samples,]$WV),
                          "Untreated" = na.omit(temp.wv.df[temp.wv.df$Sample %in% untreated.samples,]$WV))
          p.df[p.df$Signature==n,]$p <- t.test(temp.wv$Treated, temp.wv$Untreated, "less", paired=paired)$p.value
          p.df[p.df$Signature==n,"Paired"] <- paired
        }
        
        # plot
        p.df$Significance <- "p > 0.05"
        p.df$p <- as.numeric(p.df$p)
        if (any(p.df$p <= 0.05)) {
          p.df[p.df$p <= 0.05,]$Significance <- "p <= 0.05" 
        }
        p.df$Significance <- factor(p.df$Significance,
                                    levels=c("p <= 0.05", "p > 0.05"))
        p.df$minusLogP<- 1E-4
        if (any(p.df$p != 0)) {
          p.df[p.df$p != 0,]$minusLogP <- -log(p.df$p, base=10) 
        }
        sigOrder <- p.df[order(p.df$p),]$Signature
        ggplot(p.df, aes(x=Signature, y=minusLogP, 
                         fill = Signature, alpha = 0.5)) + 
          geom_col() + theme_classic(base_size = 12) + ylab("-Log(P-value)") + 
          ggplot2::scale_x_discrete(limits = sigOrder) +
          scale_fill_manual(values=RColorBrewer::brewer.pal(length(Signature), "Set2"), 
                            breaks=c("Sorted","Lasry","Triana","van Galen", "Overall"))+
          ggtitle(paste("T-test: Monocyte scores are lower in",j, "treated samples"))
        ggsave(paste0(i,"_", j, "_", k, "_pValue_WV.pdf"), width = 5, height = 5)
        
        p.df$Contrast <- j
        p.df$Input <- paste(i, k)
        contrast.p.df <- rbind(contrast.p.df, p.df)
      }
      exp.p.df <- rbind(exp.p.df, contrast.p.df)
    }
  }
}
write.csv(exp.p.df, "NRAS_pValue_WV.csv", row.names = FALSE)

# plot for each contrast
contrasts <- c("NRAS", "NRAS_w_Gilt_ctl", "Gilt", "Gilt_w_NRAS_ctl")
for (j in contrasts) {
  p.df <- exp.p.df[exp.p.df$Contrast == j & !grepl(" HR",exp.p.df$Input),]
  p.df <- p.df[order(p.df$p),]
  p.df$Signature <- factor(p.df$Signature, levels=unique(sigOrder))
  ggplot2::ggplot(p.df, aes(x=Signature, y=minusLogP)) + geom_violin(alpha=0) +
    geom_point(aes(color=Input)) + geom_boxplot(width=0.2, alpha = 0) + 
    labs(y="-Log(P-value)") + theme_classic(base_size = 12) +
    ggtitle(paste("T-test: Monocyte scores are lower in",j, "treated samples"))
  ggsave(paste0(j,"_pValue_WV.pdf"), width = 5, height = 5) 
}

# also make plots just for exp 20 timecourse
exp.p.df <- read.csv("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Exp21_NRAS_ASO_treated_patients/proteomics/analysis/Monocyte_vs_progenitor_signatures_20250530/NRAS_pValue_WV.csv")
p.df <- exp.p.df[exp.p.df$Contrast == "NRAS" & grepl(" HR",exp.p.df$Input),]
p.df <- p.df[order(p.df$p),]
p.df$Timepoint <- sub("Cell Line ASO*.", "", p.df$Input)
sigOrder <- p.df[order(p.df$p),]$Signature
p.df$Signature <- factor(p.df$Signature, levels=unique(sigOrder))
ggplot2::ggplot(p.df, aes(x=Signature, y=minusLogP)) + geom_violin(alpha=0) +
  geom_point(aes(color=Timepoint)) + geom_boxplot(width=0.2, alpha = 0) + 
  labs(y="-Log(P-value)") + theme_classic(base_size = 12) +
  ggtitle(paste("T-test: Monocyte scores are lower in","NRAS", "treated samples"))
ggsave(paste0("NRAS_timecourse","_pValue_WV.pdf"), width = 5, height = 5) 

#### try deconvolution with just van Galen markers ####
### get van Galen markers
synapser::synLogin()
vg <- readxl::read_excel(synapser::synGet("syn61831560")$path)
vg.tumor <- vg[,8:(ncol(vg)-1)]
colnames(vg.tumor) <- vg.tumor[1,]
vg.tumor <- vg.tumor[2:(nrow(vg.tumor)-2),]

# prepare progenitor markers
vg.prog <- as.data.frame(vg.tumor[,"Progenitor-like"])
colnames(vg.prog)[1] <- "Gene"
vg.prog$weight <- -1

# prepare monocytic markers
vg.mono <- as.data.frame(vg.tumor[,"Monocyte-like"])

colnames(vg.mono)[1] <- "Gene"
vg.mono$weight <- 1

# combine prog and mono markers
vg.prog.mono <- rbind(vg.prog, vg.mono)
colnames(vg.prog.mono)[2] <- "Log2FC"
sigs <- list("van Galen" = vg.prog.mono)

### weighted voting deconvolution
setwd("~/OneDrive - PNNL/Documents/GitHub/Exp21_NRAS_ASO_treated_patients/proteomics/")
dir.create("analysis")
setwd("analysis")
base.path <- "~/OneDrive - PNNL/Documents/GitHub/Exp21_NRAS_ASO_treated_patients/proteomics/analysis/"

# pt <- list("Meta" = "syn64732333", # exp 21
#            "Global" = "syn52751324",
#            "Phospho" = "syn52751606")
# cell <- list("Meta" = "syn52653527", # exp 20
#              "Global" = "syn52623904",
#              "Phospho" = "syn52623908")
# cellKO <- list("Meta" = "syn64731726", # exp 22
#                "Global" = "syn53212316",
#                "Phospho" = "syn53212321")

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
#synIDs <- c("syn53180709","syn52256487","syn53212344")
synIDs <- c("syn66726942","syn66726775","syn66705127")
names(synIDs) <- names(inputs)
exp.p.df <- data.frame()
library(ggplot2)
for (i in names(inputs)) {
  cat("Running",i,"\n")
  
  ### read in data
  meta.df <- readxl::read_excel(synapser::synGet(inputs[[i]]$Meta)$path, 
                                sheet = meta.sheet[[i]])
  global.df <- read.table(synapser::synGet(inputs[[i]]$Global)$path, sep = "\t")
  
  # remove samples that Sunil said to exclude
  if (grepl("patient", i, ignore.case = TRUE)) {
    excluded.patients <- c("20-00576", "20-00432", "18-00101", "17-01083", "16-00815", "16-01253")
    meta.df <- meta.df[!(meta.df$PatientID %in% excluded.patients),] 
  }
  
  # add format data so that samples are along rows
  global.df <- as.data.frame(t(global.df))
  global.df$Sample <- rownames(global.df)
  global.df <- global.df[,c("Sample", colnames(global.df)[1:(ncol(global.df)-1)])]
  if (grepl("patient", i, ignore.case = TRUE)) {
    global.df <- global.df[!(global.df$Sample %in% excluded.patients),]
  }
  global.df100 <- global.df[,colSums(is.na(global.df)) == 0]
  
  setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Exp21_NRAS_ASO_treated_patients/proteomics/analysis/")
  dir.create("Monocyte_vs_progenitor_signatures_20250530")
  setwd("Monocyte_vs_progenitor_signatures_20250530")
  dir.create("vanGalenBinary")
  setwd("vanGalenBinary")
  
  ### evaluate signature
  wv.df <- compareSigs(global.df100, sigs) 
  write.csv(wv.df, paste0(i,"_wv.csv"), row.names = FALSE)
  
  ### is monocytic score lower in treated vs. untreated samples?
  ## prep contrasts
  temp.meta <- list()
  if (grepl("patient", i, ignore.case = TRUE)) { # exp 21
    contrasts <- c("NRAS", "NRAS_w_Gilt_ctl", "Gilt", "Gilt_w_NRAS_ctl")
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
  } else if (grepl("ASO", i, ignore.case = TRUE)) { # exp 20
    contrasts <- c("NRAS", "NRAS_w_Gilt_ctl", "Gilt", "Gilt_w_NRAS_ctl")
    # NRAS ASO vs. NEG ASO
    aso.treatments <- c("NRAS ASO", "CTRL ASO")
    meta.df.nras <- as.data.frame(meta.df[meta.df$Group %in% aso.treatments,])
    meta.df.nras$NRAS <- "Untreated"
    meta.df.nras[grepl("NRAS ASO", meta.df.nras$Sample),]$NRAS <- "Treated"
    rownames(meta.df.nras) <- meta.df.nras$MeasurementName
    temp.meta[[contrasts[1]]] <- meta.df.nras
    
    # Gilt + NRAS vs. Gilt
    nras.gilt.treatments <- c("NRAS ASO + GILT", "CTRL ASO + GILT")
    meta.df.nras.gilt <- as.data.frame(meta.df[meta.df$Group %in% nras.gilt.treatments,])
    meta.df.nras.gilt$NRAS_w_Gilt_ctl <- "Untreated"
    meta.df.nras.gilt[meta.df.nras.gilt$Group == "NRAS ASO + GILT",]$NRAS_w_Gilt_ctl <- "Treated"
    rownames(meta.df.nras.gilt) <- meta.df.nras.gilt$MeasurementName
    temp.meta[[contrasts[2]]] <- meta.df.nras.gilt
    
    # Gilt + NRAS vs. Gilt
    gilt.treatments <- c("CTRL ASO + GILT", "CTRL ASO")
    meta.df.gilt <- as.data.frame(meta.df[meta.df$Group %in% gilt.treatments,])
    meta.df.gilt$Gilt <- "Untreated"
    meta.df.gilt[meta.df.gilt$Group == "CTRL ASO + GILT",]$Gilt <- "Treated"
    rownames(meta.df.gilt) <- meta.df.gilt$MeasurementName
    temp.meta[[contrasts[3]]] <- meta.df.gilt
    
    # Gilt + NRAS vs. Gilt
    nras.gilt.treatments <- c("NRAS ASO + GILT", "NRAS ASO")
    meta.df.nras.gilt <- as.data.frame(meta.df[meta.df$Group %in% nras.gilt.treatments,])
    meta.df.nras.gilt$Gilt_w_NRAS_ctl <- "Untreated"
    meta.df.nras.gilt[meta.df.nras.gilt$Group == "NRAS ASO + GILT",]$Gilt_w_NRAS_ctl <- "Treated"
    rownames(meta.df.nras.gilt) <- meta.df.nras.gilt$MeasurementName
    temp.meta[[contrasts[4]]] <- meta.df.nras.gilt
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
  }
  
  Signature <- c("Overall",na.omit(unique(wv.df$Signature)))
  contrast.p.df <- data.frame()
  ## for each contrast:
  for (j in contrasts) {
    setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Exp21_NRAS_ASO_treated_patients/proteomics/analysis/")
    dir.create("Monocyte_vs_progenitor_signatures_20250530")
    setwd("Monocyte_vs_progenitor_signatures_20250530")
    dir.create("vanGalenBinary")
    setwd("vanGalenBinary")
    # t-test
    p.df <- data.frame(Signature, p=NA, Paired=FALSE)
    treated.samples <- rownames(temp.meta[[j]][temp.meta[[j]][,j] == "Treated",])
    untreated.samples <- rownames(temp.meta[[j]][temp.meta[[j]][,j] == "Untreated",])
    paired <- all(na.omit(wv.df[wv.df$Sample %in% treated.samples,]$Sample) == na.omit(wv.df[wv.df$Sample %in% untreated.samples,]$Sample))
    temp.wv <- list("Treated" = na.omit(wv.df[wv.df$Sample %in% treated.samples,]$WV),
                    "Untreated" = na.omit(wv.df[wv.df$Sample %in% untreated.samples,]$WV))
    p.df[p.df$Signature=="Overall",]$p <- t.test(temp.wv$Treated, temp.wv$Untreated, "less", paired=paired)$p.value
    p.df[p.df$Signature=="Overall","Paired"] <- paired
    wv.df$Treatment <- NA
    if (grepl("_w_", j)) {
      treatmentDescr <- sub("_w_", " + ", j)
      treatmentDescr <- sub("_ctl","",treatmentDescr)
      ctlDescr <- strsplit(treatmentDescr, " [+] ", )[[1]][[2]]
    } else {
      treatmentDescr <- j
      ctlDescr <- "Control"
    }
    wv.df[wv.df$Sample %in% treated.samples,]$Treatment <- treatmentDescr
    wv.df[wv.df$Sample %in% untreated.samples,]$Treatment <- ctlDescr
    if (i == "Cell Line ASO") {
      rownames(meta.df) <- meta.df$MeasurementName
      wv.df$Time <- NA
      timepoints <- na.omit(unique(meta.df$Time[meta.df$Time != "baseline"]))
      for (q in timepoints) {
        time.meta <- na.omit(meta.df[meta.df$Time == q,])
        wv.df[wv.df$Sample %in% time.meta$MeasurementName,]$Time <- q
      }
    }
    
    # repeat t-test for each signature
    other.sigs <- na.omit(unique(wv.df$Signature))
    for (k in other.sigs) {
      temp.wv.df <- wv.df[wv.df$Signature==k,]
      paired <- all(na.omit(temp.wv.df[temp.wv.df$Sample %in% treated.samples,]$Sample) == na.omit(temp.wv.df[temp.wv.df$Sample %in% untreated.samples,]$Sample))
      temp.wv <- list("Treated" = na.omit(temp.wv.df[temp.wv.df$Sample %in% treated.samples,]$WV),
                      "Untreated" = na.omit(temp.wv.df[temp.wv.df$Sample %in% untreated.samples,]$WV))
      p.df[p.df$Signature==k,]$p <- t.test(temp.wv$Treated, temp.wv$Untreated, "less", paired=paired)$p.value
      p.df[p.df$Signature==k,"Paired"] <- paired
      if (i == "Cell Line ASO") {
        ggplot2::ggplot(na.omit(temp.wv.df), aes(x=Treatment, y = WV)) + geom_violin(alpha=0) + 
          geom_point(aes(color=Time)) + geom_boxplot(width=0.2, alpha=0) + 
          labs(y="Monocytic vs. Progenitor-like Score") + 
          theme_classic(base_size = 12) + 
          ggtitle(paste0(i, " Treated Samples are\nLess Monocytic (p = ", round(p.df[p.df$Signature==k,]$p,2),")"))
      } else {
        ggplot2::ggplot(na.omit(temp.wv.df), aes(x=Treatment, y = WV)) + geom_violin(alpha=0) + 
          geom_point() + geom_boxplot(width=0.2, alpha=0) + 
          labs(y="Monocytic vs. Progenitor-like Score") + 
          theme_classic(base_size = 12) + 
          ggtitle(paste0(i, " Treated Samples are\nLess Monocytic (p = ", round(p.df[p.df$Signature==k,]$p,2),")"))
      }
      ggsave(paste0(i, "_", j,"_",k,"_WV.pdf"), width = 5, height = 5)
    }
    if (i == "Cell Line ASO") {
      ggplot2::ggplot(na.omit(wv.df), aes(x=Treatment, y = WV)) + geom_violin(alpha=0) + 
        geom_point(aes(color = Signature, shape = Time)) + geom_boxplot(width=0.2, alpha=0) + 
        labs(y="Monocytic vs. Progenitor-like Score") + 
        theme_classic(base_size = 12) + 
        ggtitle(paste0(i, " Treated Samples are\nLess Monocytic (p = ", round(p.df[p.df$Signature=="Overall",]$p,2),")")) 
    } else {
      ggplot2::ggplot(na.omit(wv.df), aes(x=Treatment, y = WV)) + geom_violin(alpha=0) + 
        geom_point(aes(color = Signature)) + geom_boxplot(width=0.2, alpha=0) + 
        labs(y="Monocytic vs. Progenitor-like Score") + 
        theme_classic(base_size = 12) + 
        ggtitle(paste0(i, " Treated Samples are\nLess Monocytic (p = ", round(p.df[p.df$Signature=="Overall",]$p,2),")"))
    }
    ggsave(paste0(i, "_", j,"_overall_WV.pdf"), width = 5, height = 5)
    
    # plot
    p.df$Significance <- "p > 0.05"
    p.df$p <- as.numeric(p.df$p)
    if (any(p.df$p <= 0.05)) {
      p.df[p.df$p <= 0.05,]$Significance <- "p <= 0.05" 
    }
    p.df$Significance <- factor(p.df$Significance,
                                levels=c("p <= 0.05", "p > 0.05"))
    p.df$minusLogP<- 1E-4
    if (any(p.df$p != 0)) {
      p.df[p.df$p != 0,]$minusLogP <- -log(p.df$p, base=10) 
    }
    sigOrder <- p.df[order(p.df$p),]$Signature
    ggplot(p.df, aes(x=Signature, y=minusLogP, 
                     fill = Signature, alpha = 0.5)) + 
      geom_col() + theme_classic(base_size = 12) + ylab("-Log(P-value)") + 
      ggplot2::scale_x_discrete(limits = sigOrder) +
      scale_fill_manual(values=RColorBrewer::brewer.pal(length(Signature), "Set2"), 
                        breaks=c("Sorted","Lasry","Triana","van Galen", "Overall"))+
      ggtitle(paste("T-test: Monocyte scores are lower in",j, "treated samples"))
    ggsave(paste0(i,"_", j, "_pValue_WV.pdf"), width = 5, height = 5)
    
    p.df$Contrast <- j
    contrast.p.df <- rbind(contrast.p.df, p.df)
  }
  contrast.p.df$Input <- i
  exp.p.df <- rbind(exp.p.df, contrast.p.df)
  
  ## repeat with timepoint groupings if exp 20
  if (i == "Cell Line ASO") {
    timepoints <- na.omit(unique(meta.df$Time[meta.df$Time != "baseline"]))
    for (k in timepoints) {
      temp.meta.df <- meta.df[meta.df$Time == k,]
      temp.meta <- list()
      contrasts <- c("NRAS", "NRAS_w_Gilt_ctl", "Gilt", "Gilt_w_NRAS_ctl")
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
      
      for (j in contrasts) {
        setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Exp21_NRAS_ASO_treated_patients/proteomics/analysis/")
        dir.create("Monocyte_vs_progenitor_signatures_20250530")
        setwd("Monocyte_vs_progenitor_signatures_20250530")
        dir.create("vanGalenBinary")
        setwd("vanGalenBinary")
        # t-test
        p.df <- data.frame(Signature, p=NA, Paired=FALSE)
        treated.samples <- rownames(temp.meta[[j]][temp.meta[[j]][,j] == "Treated",])
        untreated.samples <- rownames(temp.meta[[j]][temp.meta[[j]][,j] == "Untreated",])
        paired <- all(na.omit(wv.df[wv.df$Sample %in% treated.samples,]$Sample) == na.omit(wv.df[wv.df$Sample %in% untreated.samples,]$Sample))
        temp.wv <- list("Treated" = na.omit(wv.df[wv.df$Sample %in% treated.samples,]$WV),
                        "Untreated" = na.omit(wv.df[wv.df$Sample %in% untreated.samples,]$WV))
        p.df[p.df$Signature=="Overall",]$p <- t.test(temp.wv$Treated, temp.wv$Untreated, "less", paired=paired)$p.value
        p.df[p.df$Signature=="Overall","Paired"] <- paired
        
        # repeat t-test for each signature
        other.sigs <- na.omit(unique(wv.df$Signature))
        for (n in other.sigs) {
          temp.wv.df <- wv.df[wv.df$Signature==n,]
          paired <- all(na.omit(temp.wv.df[temp.wv.df$Sample %in% treated.samples,]$Sample) == na.omit(temp.wv.df[temp.wv.df$Sample %in% untreated.samples,]$Sample))
          temp.wv <- list("Treated" = na.omit(temp.wv.df[temp.wv.df$Sample %in% treated.samples,]$WV),
                          "Untreated" = na.omit(temp.wv.df[temp.wv.df$Sample %in% untreated.samples,]$WV))
          p.df[p.df$Signature==n,]$p <- t.test(temp.wv$Treated, temp.wv$Untreated, "less", paired=paired)$p.value
          p.df[p.df$Signature==n,"Paired"] <- paired
          ggplot2::ggplot(na.omit(temp.wv.df), aes(x=Treatment, y = WV)) + geom_violin(alpha=0) + 
            geom_point() + geom_boxplot(width=0.2, alpha=0) + 
            labs(y="Monocytic vs. Progenitor-like Score") + 
            theme_classic(base_size = 12) + 
            ggtitle(paste0(i, " Treated Samples are\nLess Monocytic (p = ", round(p.df[p.df$Signature==k,]$p,2),")"))
          ggsave(paste0(i, "_", j,"_", k, "_",n,"_WV.pdf"), width = 5, height = 5)
        }
        ggplot2::ggplot(na.omit(wv.df), aes(x=Treatment, y = WV)) + geom_violin(alpha=0) + 
          geom_point(aes(color = Signature)) + geom_boxplot(width=0.2, alpha=0) + 
          labs(y="Monocytic vs. Progenitor-like Score") + 
          theme_classic(base_size = 12) + 
          ggtitle(paste0(i, " Treated Samples are\nLess Monocytic (p = ", round(p.df[p.df$Signature=="Overall",]$p,2),")"))
        ggsave(paste0(i, "_", j,"_", k, "_overall_WV.pdf"), width = 5, height = 5)
        
        # plot
        p.df$Significance <- "p > 0.05"
        p.df$p <- as.numeric(p.df$p)
        if (any(p.df$p <= 0.05)) {
          p.df[p.df$p <= 0.05,]$Significance <- "p <= 0.05" 
        }
        p.df$Significance <- factor(p.df$Significance,
                                    levels=c("p <= 0.05", "p > 0.05"))
        p.df$minusLogP<- 1E-4
        if (any(p.df$p != 0)) {
          p.df[p.df$p != 0,]$minusLogP <- -log(p.df$p, base=10) 
        }
        sigOrder <- p.df[order(p.df$p),]$Signature
        ggplot(p.df, aes(x=Signature, y=minusLogP, 
                         fill = Signature, alpha = 0.5)) + 
          geom_col() + theme_classic(base_size = 12) + ylab("-Log(P-value)") + 
          ggplot2::scale_x_discrete(limits = sigOrder) +
          scale_fill_manual(values=RColorBrewer::brewer.pal(length(Signature), "Set2"), 
                            breaks=c("Sorted","Lasry","Triana","van Galen", "Overall"))+
          ggtitle(paste("T-test: Monocyte scores are lower in",j, "treated samples"))
        ggsave(paste0(i,"_", j, "_", k, "_pValue_WV.pdf"), width = 5, height = 5)
        
        p.df$Contrast <- j
        p.df$Input <- paste(i, k)
        contrast.p.df <- rbind(contrast.p.df, p.df)
      }
      exp.p.df <- rbind(exp.p.df, contrast.p.df)
    }
  }
}
write.csv(exp.p.df, "NRAS_pValue_WV_vanGalenBinary.csv", row.names = FALSE)

# plot for each contrast
contrasts <- c("NRAS", "NRAS_w_Gilt_ctl", "Gilt", "Gilt_w_NRAS_ctl")
for (j in contrasts) {
  p.df <- exp.p.df[exp.p.df$Contrast == j & !grepl(" HR",exp.p.df$Input),]
  p.df <- p.df[order(p.df$p),]
  p.df$Signature <- factor(p.df$Signature, levels=unique(sigOrder))
  ggplot2::ggplot(p.df, aes(x=Signature, y=minusLogP)) + geom_violin(alpha=0) +
    geom_point(aes(color=Input)) + geom_boxplot(width=0.2, alpha = 0) + 
    labs(y="-Log(P-value)") + theme_classic(base_size = 12) +
    ggtitle(paste("T-test: Monocyte scores are lower in",j, "treated samples"))
  ggsave(paste0(j,"_pValue_WV.pdf"), width = 5, height = 5) 
}

# also make plots just for exp 20 timecourse
exp.p.df <- read.csv("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Exp21_NRAS_ASO_treated_patients/proteomics/analysis/Monocyte_vs_progenitor_signatures_20250530/vanGalenBinary/NRAS_pValue_WV_vanGalenBinary.csv")
p.df <- exp.p.df[exp.p.df$Contrast == "NRAS" & grepl(" HR",exp.p.df$Input),]
p.df <- p.df[order(p.df$p),]
p.df$Timepoint <- sub("Cell Line ASO*.", "", p.df$Input)
sigOrder <- p.df[order(p.df$p),]$Signature
p.df$Signature <- factor(p.df$Signature, levels=unique(sigOrder))
ggplot2::ggplot(p.df, aes(x=Signature, y=minusLogP)) + geom_violin(alpha=0) +
  geom_point(aes(color=Timepoint)) + geom_boxplot(width=0.2, alpha = 0) + 
  labs(y="-Log(P-value)") + theme_classic(base_size = 12) +
  ggtitle(paste("T-test: Monocyte scores are lower in","NRAS", "treated samples"))
ggsave(paste0("NRAS_timecourse","_pValue_WV.pdf"), width = 5, height = 5) 


