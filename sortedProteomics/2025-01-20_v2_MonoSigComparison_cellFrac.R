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
synapser::synLogin()

source("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/helperScripts/circBar.R")
source("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/helperScripts/panSEA_helper_20240913.R")
evalOneMonoSig <- function(global.df100, frac.df, temp.sig, BeatAML, type="rna", gmt) {
  temp.sig <- na.omit(temp.sig)
  if (nrow(temp.sig) > 0) {
    ## sorted proteomics
    # perform weighted voting on sorted proteomics
    temp.sig <- temp.sig[,c("Gene", colnames(temp.sig)[1])]
    temp.sig2 <- temp.sig[temp.sig$Gene %in% colnames(global.df100)[2:ncol(global.df100)],]
    if (nrow(temp.sig2) > 0) {
      global.df100 <- global.df100[,c("Sample",temp.sig2$Gene)]
      if (nrow(global.df100) > 0) {
        sorted.wv <- panSEA::WV(global.df100, temp.sig2)
        temp.wv <- sorted.wv$scores 
        
        # evaluate accuracy based on t-test
        temp.test <- stats::t.test(temp.wv[grepl("CD14", temp.wv$Sample),]$WV,
                                   temp.wv[!grepl("CD14", temp.wv$Sample),]$WV, 
                                   "greater")
        temp.test.df <- as.data.frame(unlist(temp.test))
        colnames(temp.test.df)[1] <- "value"
        temp.test.df$variable <- rownames(temp.test.df)
        temp.test.df$N <- nrow(temp.wv)
        temp.test.df$N_CD14_Pos <- nrow(temp.wv[grepl("CD14", temp.wv$Sample),]$WV)
        temp.test.df$N_CD34_Pos <- nrow(temp.wv[grepl("CD34", temp.wv$Sample),]$WV)
        } else {
          frac.corr.result <- data.frame()
        }
      } else {
        temp.wv <- data.frame()
        temp.test.df <- data.frame()
      }
    } else {
      temp.wv <- data.frame()
      temp.test.df <- data.frame()
    }
    
    ## Beat AML
    # perform DMEA on Beat AML global proteomics
    expr <- BeatAML[[type]]
    temp.sig2 <- temp.sig[temp.sig$Gene %in% colnames(expr)[2:ncol(expr)],]
    nSamples <- length(expr$Barcode.ID[expr$Barcode.ID %in% BeatAML$drug$Barcode.ID])
    if (nrow(temp.sig2) > 2 & nSamples > 2) {
      expr <- expr[,c("Barcode.ID", colnames(expr)[colnames(expr) %in% temp.sig2$Gene])]
      if (ncol(expr) > 3 & nrow(expr) > 2) {
        DMEA.result <- panSEA::mDMEA(BeatAML$drug, gmt, list(expr), 
                                     list(temp.sig2), types=type,
                                     sample.names="Barcode.ID",
                                     weight.values=colnames(temp.sig2)[2], 
                                     scatter.plots = FALSE)
        DMEA.result.Spearman <- panSEA::mDMEA(BeatAML$drug, gmt, list(expr), 
                                     list(temp.sig2), types=type,
                                     sample.names="Barcode.ID",
                                     weight.values=colnames(temp.sig2)[2], 
                                     scatter.plots = FALSE, 
                                     rank.metric="Spearman.est")
        corr.df <- DMEA.result$all.results[[1]]$corr.result
        
        # compare to known cell fractions
        temp.wv2 <- DMEA.result$all.results[[1]]$WV.scores
        temp.wv.frac <- merge(frac.df, temp.wv2, by="Barcode.ID")
        if (nrow(temp.wv.frac) > 2) {
          frac.corr <- cor.test(temp.wv.frac$WV, temp.wv.frac[,2], 
                                method = "pearson")
          N <- nrow(temp.wv.frac)
          Pearson.est <- frac.corr$estimate
          Pearson.p <-frac.corr$p.value
          frac.corr <- cor.test(temp.wv.frac$WV, temp.wv.frac[,2], 
                                method = "spearman")
          Spearman.est <- frac.corr$estimate
          Spearman.p <-frac.corr$p.value
          frac.corr.result <- data.frame(Pearson.est, Pearson.p, 
                                         Spearman.est, Spearman.p, N)
      } else {
        corr.df <- data.frame()
        DMEA.result <- list()
        DMEA.result.Spearman <- list()
        frac.corr.result <- data.frame()
      }
    } else {
      corr.df <- data.frame()
      DMEA.result <- list()
      DMEA.result.Spearman <- list()
      frac.corr.result <- data.frame()
    }
  } else {
    temp.wv <- data.frame()
    temp.test.df <- data.frame()
    corr.df <- data.frame()
    frac.corr.result <- data.frame()
    DMEA.result <- list()
    DMEA.result.Spearman <- list()
  }
  
  return(list(wv = temp.wv, test = temp.test.df, drug.corr = corr.df, 
              frac.corr = frac.corr.result, 
              DMEA = DMEA.result, DMEA.Spearman=DMEA.result.Spearman))
}

evalMonoSig <- function(global.df100, frac.df, sig.matrix, BeatAML, 
                        types=rep("rna",ncol(sig.matrix)), gmt) {
  # evaluate each signature
  wv.df <- data.frame()
  test.df <- data.frame()
  drug.corr.df <- data.frame()
  frac.corr.df <- data.frame()
  sig.matrix$Gene <- rownames(sig.matrix)
  DMEA.results <- list()
  DMEA.results.Spearman <- list()
  for (i in 1:(ncol(sig.matrix)-1)) {
    cat("evaluating",names(sig.matrix)[i],"as",types[i],"\n")
    temp.result <- evalOneMonoSig(global.df100, frac.df, 
                                  sig.matrix[,c(i,ncol(sig.matrix))],
                                  BeatAML, types[i], gmt)
    DMEA.results[[names(sig.matrix)[i]]] <- temp.result$DMEA
    DMEA.results.Spearman[[names(sig.matrix)[i]]] <- temp.result$DMEA.Spearman
    temp.wv.df <- temp.result$wv
    temp.test.df <- temp.result$test
    temp.drug.corr.df <- temp.result$drug.corr
    temp.frac.corr.df <- temp.result$frac.corr
    temp.wv.df$Signature <- names(sig.matrix)[i]
    temp.test.df$Signature <- names(sig.matrix)[i]
    temp.drug.corr.df$Signature <- names(sig.matrix)[i]
    temp.frac.corr.df$Signature <- names(sig.matrix)[i]
    
    wv.df <- rbind(wv.df, temp.wv.df)
    test.df <- rbind(test.df, temp.test.df)
    drug.corr.df <- rbind(drug.corr.df, temp.drug.corr.df)
    frac.corr.df <- rbind(frac.corr.df, temp.frac.corr.df)
  }
  p.df <- test.df[test.df$variable == "p.value",]
  
  return(list(wv = wv.df, p = p.df, drug.corr = drug.corr.df, 
              frac.corr = frac.corr.df, DMEA = DMEA.results,
              DMEA.Spearman = DMEA.results.Spearman))
}

compareSigs <- function(global.df100, frac.df, sigs, 
                        value.var = "Log2FC", BeatAML, types=rep("rna",length(sigs)), gmt, 
                        fillVals = RColorBrewer::brewer.pal(length(sigs), "Set2")) {
  # combine signatures into matrix
  filtered.sigs.df <- data.table::rbindlist(sigs, use.names = TRUE, idcol = "Signature")
  sig.matrix <- reshape2::dcast(filtered.sigs.df, Gene ~ Signature, mean,
                                value.var = value.var)
  rownames(sig.matrix) <- sig.matrix$Gene
  sig.matrix$Gene <- NULL
  sig.matrix <- sig.matrix[,names(sigs)]
  
  # test signature matrix
  sigResults <- evalMonoSig(global.df100, frac.df, sig.matrix, BeatAML,types,gmt)
  wv.df <- sigResults$wv
  p.df <- sigResults$p
  drug.corr.df <- sigResults$drug.corr
  frac.corr.df <- sigResults$frac.corr
  
  # plot results
  p.df$Significance <- "p > 0.05"
  p.df$value <- as.numeric(p.df$value)
  if (any(p.df$value <= 0.05)) {
    p.df[p.df$value <= 0.05,]$Significance <- "p <= 0.05" 
  }
  p.df$Significance <- factor(p.df$Significance,
                                      levels=c("p <= 0.05", "p > 0.05"))
  p.df$minusLogFDR <- 1E-4
  if (any(p.df$value != 0)) {
    p.df[p.df$value != 0,]$minusLogFDR <- -log(p.df$value, base=10) 
  }
  sigOrder <- p.df[order(p.df$value),]$Signature
  ggplot(p.df, aes(x=Signature, y=-log(value, base=10), 
                   fill = Signature, alpha = 0.5)) + 
    geom_col() + theme_classic(base_size = 12) + ylab("-Log(P-value)") + 
    ggplot2::scale_x_discrete(limits = sigOrder) +
    scale_fill_manual(values=fillVals, 
                      breaks=c("Sorted","Lasry","Triana","van Galen"))+
    ggtitle("T-test: Monocyte scores are higher in CD14+ samples")
  ggsave("pValue_DIA_WV_signatureFill.pdf", width = 5, height = 5)
  
  p.df$`-Log(P-value)` <- -log(p.df$value, base=10)
  plot.df <- p.df
  plot.df$alpha <- 0.5
  circBar(plot.df, x="Signature", y = "-Log(P-value)", fill = "Signature", 
          alpha = "alpha", ymin = 0, ymax = 10, alpha_range=0.5, 
          y_ticks=seq(0,10,length.out=6),
          fillVals=fillVals,
          title="Monocytic signatures distinguish CD14+ and CD34+ samples",
          fname="pValue_DIA_WV_signatureFill_circBarPlot.pdf")
  
  rank.metrics <- c("Pearson.est", "Spearman.est")
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
      #sigOrder <- av.df[order(av.df$rank, decreasing=TRUE),]$Signature
      ggplot(na.omit(plot.df), aes(x=Signature, y=rank, fill = Signature, alpha=0.5)) + 
        geom_col() + theme_minimal(base_size = 12) + ylab(paste0(descr," Correlation Estimate")) + 
        facet_wrap(~ `Drug Treatment`) +
        ggplot2::scale_x_discrete(limits = sigOrder) +
        theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
        scale_fill_manual(values=fillVals, 
                          breaks=c("Sorted","Lasry","Triana","van Galen"))+
        ggtitle("Monocytic signatures predict drug sensitivity")
      ggsave(paste0("drugCorr_DIA_WV_signatureFill_",descr,".pdf"), width = 5, height = 5)
      
      for (j in doi.names) {
        plot.df <- drug.corr.df[drug.corr.df$`Drug Treatment` == j,]
        plot.df$rank <- plot.df[,i]
        sigOrder <- na.omit(unique(plot.df[order(plot.df$rank, decreasing=TRUE),]$Signature))
        plot.df$alpha <- 0.5
        circBar(na.omit(plot.df), x="Signature", y = "rank", fill = "Signature", 
                alpha = "alpha", ymin = 0, ymax = 1, alpha_range=0.5, 
                ytick_yScale = 2/3, ytick_yShift = 0, fillVals=fillVals,
                title=paste("Monocytic signatures predict", j, "sensitivity"),
                fname=paste0(j,"_Corr_DIA_WV_signatureFill_circBarPlot_",descr,".pdf"))
        ggplot(plot.df, aes(x=Signature, y=rank, fill = Signature, alpha=0.5)) + 
          geom_col() + theme_classic(base_size = 12) + 
          ylab(paste(descr, "Correlation Estimate")) + 
          ggplot2::scale_x_discrete(limits = sigOrder) +
          scale_fill_manual(values=fillVals, 
                            breaks=c("Sorted","Lasry","Triana","van Galen"))+
          ggtitle(paste("Monocytic signatures predict", j, "sensitivity"))
        ggsave(paste0(j,"_Corr_DIA_WV_signatureFill_barPlot_",descr,".pdf"), width = 5, height = 5)
      }
    }
    
    plot.df <- frac.corr.df
    plot.df$rank <- plot.df[,i]
    sigOrder <- plot.df[order(plot.df$rank, decreasing=TRUE),]$Signature
    ggplot(plot.df, aes(x=Signature, y=rank, fill = Signature, alpha=0.5)) + 
      geom_col() + theme_classic(base_size = 12) + ylab(paste(descr, "Correlation Estimate")) + 
      ggplot2::scale_x_discrete(limits = sigOrder) +
      scale_fill_manual(values=fillVals, 
                        breaks=c("Sorted","Lasry","Triana","van Galen"))+
      ggtitle("Monocytic signatures predict monocyte fraction")
    ggsave(paste0("fracCorr_DIA_WV_",descr,"_signatureFill.pdf"), width = 5, height = 5)
    plot.df$alpha <- 0.5
    circBar(plot.df, x="Signature", y = "rank", fill = "Signature", 
            alpha = "alpha", ymin = 0, ymax = 1, alpha_range=0.5, 
            ytick_yScale = 2/3, ytick_yShift = 0, fillVals=fillVals,
            title="Monocytic signatures predict monocyte fraction",
            fname=paste0("fracCorr_DIA_WV_signatureFill_circBarPlot_",descr,".pdf"))
  }
  
  return(list(wv = wv.df, p = p.df, drug.corr = drug.corr.df, 
              frac.corr = frac.corr.df, DMEA = sigResults$DMEA,
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

evalOneSigVen <- function(global.df100, temp.sig, BeatAML, type="global", gmt) {
  temp.sig <- na.omit(temp.sig)
  
  ## Beat AML
  # perform DMEA on Beat AML global proteomics
  expr <- BeatAML[[type]]
  temp.sig2 <- temp.sig[temp.sig$Gene %in% colnames(expr)[2:ncol(expr)],]
  nSamples <- length(expr$Barcode.ID[expr$Barcode.ID %in% BeatAML$drug$Barcode.ID])
  if (nrow(temp.sig2) > 2 & nSamples > 2) {
    expr <- expr[,c("Barcode.ID", colnames(expr)[colnames(expr) %in% temp.sig2$Gene])]
    if (ncol(expr) > 3 & nrow(expr) > 2) {
      DMEA.result <- panSEA::mDMEA(BeatAML$drug, gmt, list(expr), 
                                   list(temp.sig2), types=type,
                                   sample.names="Barcode.ID",
                                   weight.values=colnames(temp.sig2)[2], 
                                   scatter.plots = FALSE)
      corr.df <- DMEA.result$all.results[[1]]$corr.result
    } else {
      corr.df <- data.frame()
    }
  } else {
    corr.df <- data.frame()
  }
  
  return(corr.df)
}


optSig <- function(global.df100, temp.sig, BeatAML, type="global", gmt) {
  genes <- temp.sig$Gene
  # leave out each gene and then predict ven AUC
  sensPred <- data.frame()
  for (i in genes) {
    temp.sig2 <- temp.sig[temp.sig$Gene != i,]
    temp.sens <- evalOneSigVen(global.df100, temp.sig2, BeatAML, type, gmt)
    temp.sens$GeneLeftOut <- i
    sensPred <- rbind(sensPred, temp.sens)
  }
  write.csv(sensPred, "drugSensPrediction_sortedBead_geneLOO.csv", row.names = FALSE)
  venPred <- sensPred[sensPred$Drug == "Venetoclax",] # all significant
  venPred$delta <- venPred$Pearson.est - 0.704114664863842
  write.csv(venPred, "venSensPrediction_sortedBead_geneLOO.csv", row.names = FALSE)
  
  # filter for genes where ven correlation drops without the gene (delta < 0)
  venImpr <- venPred[venPred$delta<0,] # 1280 / 2597
  
  # start with gene w biggest impact (i.e., most negative delta) and increase # of genes until ven Pearson q < 0.05
  useful <- venImpr[order(venImpr$delta),]$Gene # 1280
  
  # # check that selected genes are fully covered in ven samples
  # # filter for genes with 3+ data points in Beat AML in samples with ven AUC
  # venAUC <- BeatAML$drug[!is.na(BeatAML$drug$Venetoclax),c("Barcode.ID","Venetoclax")]
  # venAUCprot <- BeatAML$global[BeatAML$global$Barcode.ID %in% venAUC$Barcode.ID,] # 9413 gene symbols
  # colSumsVen <- colSums(is.na(venAUCprot))
  # venAUCprot2 <- venAUCprot[,colSumsVen==0] # 6945 gene symbols including CD14, CD34
  # useful2 <- useful[useful %in% colnames(venAUCprot2)[2:ncol(venAUCprot2)]] # 1049; 3 most useful genes kept from before this filter
  
  # nValVen <- 122 - colSumsVen
  # any(nValVen < 3) # TRUE
  # nValVen <- as.data.frame(nValVen)
  # nValVen$Gene <- rownames(nValVen)
  # nValVen[nValVen$Gene ]
  # nValVen2 <- as.data.frame(nValVen[nValVen$nValVen == 122,])
  
  usedGenes <- c()
  minSensPred <- data.frame()
  for (i in useful) {
    usedGenes <- c(usedGenes, i)
    temp.sig2 <- temp.sig[temp.sig$Gene %in% usedGenes,]
    temp.sens <- evalOneSigVen(global.df100, temp.sig2, BeatAML, type, gmt)
    if (is.data.frame(temp.sens) & nrow(temp.sens) > 0){
      temp.sens$N_genes <- length(usedGenes)
      temp.sens$Genes <- paste0(usedGenes, collapse=", ")
      minSensPred <- rbind(minSensPred, temp.sens)
      if (temp.sens[temp.sens$Drug=="Venetoclax",]$Pearson.q < 0.05) {
        write.csv(temp.sens, "drugSensPrediction_sortedBead_minGenes.csv", row.names = FALSE)
        write.csv(temp.sig2, "drugSensPrediction_sortedBead_minGeneSignature.csv", row.names = FALSE)
        break
      } 
    }
  }
  # no DMEA results with KCTD12 or S100A8; first DMEA result is with SAMHD1 which is where Ven q<0.05
  # any(is.na(venImpr$delta)) returns FALSE so DMEA seems to have worked for the full loop with gene LOO
  
  # try just ven so can look at KCTD12, S100A8
  venAUC <- BeatAML$drug[!is.na(BeatAML$drug$Venetoclax),c("Barcode.ID","Venetoclax")]
  venAUCprot <- BeatAML$global[BeatAML$global$Barcode.ID %in% venAUC$Barcode.ID,] # 9413 gene symbols
  
  # calculate WV with genes
  usedGenes <- c()
  minSensPred <- data.frame()
  for (i in useful[useful %in% colnames(venAUCprot)[2:ncol(venAUCprot)]]) { # 1187 gene symbols
    usedGenes <- c(usedGenes, i)
    temp.sig2 <- temp.sig[temp.sig$Gene %in% usedGenes,]
    temp.wv <- DMEA::WV(venAUCprot, temp.sig2)$scores
    venAUC.wv <- merge(venAUC, temp.wv, by="Barcode.ID")
    temp.corr <- cor.test(venAUC.wv$Venetoclax, venAUC.wv$WV)
    temp.sens <- data.frame(Drug="Venetoclax", Pearson.est=temp.corr$estimate, Pearson.p=temp.corr$p.value,
                            N_genes = length(usedGenes), Genes = paste0(usedGenes, collapse=", "))
    if (is.data.frame(temp.sens) & nrow(temp.sens) > 0){
      minSensPred <- rbind(minSensPred, temp.sens)
      if (temp.sens[temp.sens$Drug=="Venetoclax",]$Pearson.p < 0.05) {
        write.csv(temp.sens, "drugSensPrediction_sortedBead_minGenes_v2.csv", row.names = FALSE)
        write.csv(temp.sig2, "drugSensPrediction_sortedBead_minGeneSignature_v2.csv", row.names = FALSE)
        break
      } 
    }
  }
  
  # what if just try all genes one-by-one
  indVenPred <- data.frame()
  for (i in useful[useful %in% colnames(venAUCprot)[2:ncol(venAUCprot)]]) { # 1187 gene symbols
    temp.sig2 <- temp.sig[temp.sig$Gene == i,]
    temp.wv <- DMEA::WV(venAUCprot, temp.sig2)$scores
    venAUC.wv <- merge(venAUC, temp.wv, by="Barcode.ID")
    temp.corr <- cor.test(venAUC.wv$Venetoclax, venAUC.wv$WV)
    temp.sens <- data.frame(Drug="Venetoclax", Pearson.est=temp.corr$estimate, Pearson.p=temp.corr$p.value,
                            N_genes = 1, Genes = i)
    if (is.data.frame(temp.sens) & nrow(temp.sens) > 0){
      indVenPred <- rbind(indVenPred, temp.sens)
    }
  }
  write.csv(indVenPred, "venSensPrediction_sortedBead_1Gene.csv", row.names = FALSE)
  
  indVenPred <- data.frame()
  for (i in temp.sig$Gene[temp.sig$Gene %in% colnames(venAUCprot)[2:ncol(venAUCprot)]]) { # 2504 gene symbols
    temp.sig2 <- temp.sig[temp.sig$Gene == i,]
    temp.wv <- DMEA::WV(venAUCprot, temp.sig2)$scores
    venAUC.wv <- merge(venAUC, temp.wv, by="Barcode.ID")
    if (nrow(venAUC.wv)>2) {
      temp.corr <- cor.test(venAUC.wv$Venetoclax, venAUC.wv$WV)
      temp.corr.sp <- cor.test(venAUC.wv$Venetoclax, venAUC.wv$WV, method="spearman")
      temp.sens <- data.frame(Drug="Venetoclax", Pearson.est=temp.corr$estimate, Pearson.p=temp.corr$p.value,
                              Spearman.est=temp.corr.sp$estimate, Spearman.p=temp.corr.sp$p.value,
                              N_genes = 1, Genes = i, N=nrow(venAUC.wv))
      indVenPred <- rbind(indVenPred, temp.sens)
    }
  }
  indVenPred$Pearson.q <- qvalue::qvalue(indVenPred$Pearson.p, pi0=1)$qvalue
  indVenPred$Spearman.q <- qvalue::qvalue(indVenPred$Spearman.p, pi0=1)$qvalue
  write.csv(indVenPred, "venSensPrediction_sortedBead_1GeneAll.csv", row.names = FALSE) # NCF2: 0.728081967	2.07E-21	0.722590303	0	1	NCF2	122	5.19E-18	0
  # 0.60161205	2.33E-13	0.59101979	0	1	PLBD1	122	8.85E-12	0 (rank 66 / 2504 based on Pearson.est)
  # 0.622730428	1.89E-14	0.646627299	0	1	CD14	122	1.10E-12	0 (rank 42)
  # 0.622703425	1.89E-14	0.639046861	0	1	BCL2	122	1.10E-12	0 (rank 43)
  # 0.067428548	0.460547261	0.049398422	0.588484512	1	CD34	122	0.501177897	0.621495242 (rank 2225)
  
  indVenPred <- read.csv("venSensPrediction_sortedBead_1GeneAll.csv")
  indVenPred2 <- data.frame()
  for (i in seq(1,100)) {
    topGenes <- indVenPred[indVenPred$Pearson.q <= 0.05,] %>% slice_max(Pearson.est, n=i)
    temp.sig2 <- temp.sig[temp.sig$Gene %in% topGenes$Genes,]
    temp.wv <- DMEA::WV(venAUCprot, temp.sig2)$scores
    venAUC.wv <- merge(venAUC, temp.wv, by="Barcode.ID")
    if (nrow(venAUC.wv)>2) {
      temp.corr <- cor.test(venAUC.wv$Venetoclax, venAUC.wv$WV)
      temp.corr.sp <- cor.test(venAUC.wv$Venetoclax, venAUC.wv$WV, method="spearman")
      temp.sens <- data.frame(Drug="Venetoclax", Pearson.est=temp.corr$estimate, Pearson.p=temp.corr$p.value,
                              Spearman.est=temp.corr.sp$estimate, Spearman.p=temp.corr.sp$p.value,
                              N_genes = length(topGenes$Genes), Genes = paste0(topGenes$Genes, collapse=", "), N=nrow(venAUC.wv))
      indVenPred2 <- rbind(indVenPred2, temp.sens)
    }
  }
  indVenPred2$Pearson.q <- qvalue::qvalue(indVenPred2$Pearson.p, pi0=1)$qvalue
  indVenPred2$Spearman.q <- qvalue::qvalue(indVenPred2$Spearman.p, pi0=1)$qvalue
  write.csv(indVenPred2, "venSensPrediction_sortedBead_1-100Genes.csv", row.names = FALSE)
  # Venetoclax 0.7918502 1.825468e-27 0.7834056 0 4 NCF2, FCGRT, KCTD12, CD93 122 6.183503e-26 0
  
  indVenPred3 <- data.frame()
  for (i in seq(1,100)) {
    topGenes <- indVenPred[indVenPred$Pearson.q <= 0.05,] %>% slice_max(Pearson.est, n=i)
    temp.sig2 <- temp.sig[temp.sig$Gene %in% topGenes$Genes,]
    temp.sig2$Log2FC <- 1 # give all genes equal weight
    temp.wv <- DMEA::WV(venAUCprot, temp.sig2)$scores
    venAUC.wv <- merge(venAUC, temp.wv, by="Barcode.ID")
    if (nrow(venAUC.wv)>2) {
      temp.corr <- cor.test(venAUC.wv$Venetoclax, venAUC.wv$WV)
      temp.corr.sp <- cor.test(venAUC.wv$Venetoclax, venAUC.wv$WV, method="spearman")
      temp.sens <- data.frame(Drug="Venetoclax", Pearson.est=temp.corr$estimate, Pearson.p=temp.corr$p.value,
                              Spearman.est=temp.corr.sp$estimate, Spearman.p=temp.corr.sp$p.value,
                              N_genes = length(topGenes$Genes), Genes = paste0(topGenes$Genes, collapse=", "), N=nrow(venAUC.wv))
      indVenPred3 <- rbind(indVenPred3, temp.sens)
    }
  }
  indVenPred3$Pearson.q <- qvalue::qvalue(indVenPred3$Pearson.p, pi0=1)$qvalue
  indVenPred3$Spearman.q <- qvalue::qvalue(indVenPred3$Spearman.p, pi0=1)$qvalue
  write.csv(indVenPred3, "venSensPrediction_sortedBead_1-100GenesUnweighted.csv", row.names = FALSE)
  # Venetoclax
  # 0.7953313
  # 7.423515e-28
  # 0.7847671
  # 0
  # 4
  # NCF2, FCGRT, KCTD12, CD93
  # 122
  # 2.387415e-26
  # 0
  
  # what if we rank genes by log2FC
  indVenPred2 <- data.frame()
  for (i in seq(1,100)) {
    temp.sig2 <- temp.sig %>% slice_max(abs(Log2FC), n=i)
    temp.wv <- DMEA::WV(venAUCprot, temp.sig2)$scores
    venAUC.wv <- merge(venAUC, temp.wv, by="Barcode.ID")
    if (nrow(venAUC.wv)>2) {
      temp.corr <- cor.test(venAUC.wv$Venetoclax, venAUC.wv$WV)
      temp.corr.sp <- cor.test(venAUC.wv$Venetoclax, venAUC.wv$WV, method="spearman")
      temp.sens <- data.frame(Drug="Venetoclax", Pearson.est=temp.corr$estimate, Pearson.p=temp.corr$p.value,
                              Spearman.est=temp.corr.sp$estimate, Spearman.p=temp.corr.sp$p.value,
                              N_genes = nrow(temp.sig2), Genes = paste0(temp.sig2$Gene, collapse=", "), N=nrow(venAUC.wv))
      indVenPred2 <- rbind(indVenPred2, temp.sens)
    }
  }
  indVenPred2$Pearson.q <- qvalue::qvalue(indVenPred2$Pearson.p, pi0=1)$qvalue
  indVenPred2$Spearman.q <- qvalue::qvalue(indVenPred2$Spearman.p, pi0=1)$qvalue
  write.csv(indVenPred2, "venSensPrediction_sortedBead_1-100GenesByLog2FC.csv", row.names = FALSE)
  
  indVenPred3 <- data.frame()
  for (i in seq(1,100)) {
    temp.sig2 <- temp.sig %>% slice_max(Log2FC, n=i)
    temp.sig2$Log2FC <- 1 # give all genes equal weight
    temp.wv <- DMEA::WV(venAUCprot, temp.sig2)$scores
    venAUC.wv <- merge(venAUC, temp.wv, by="Barcode.ID")
    if (nrow(venAUC.wv)>2) {
      temp.corr <- cor.test(venAUC.wv$Venetoclax, venAUC.wv$WV)
      temp.corr.sp <- cor.test(venAUC.wv$Venetoclax, venAUC.wv$WV, method="spearman")
      temp.sens <- data.frame(Drug="Venetoclax", Pearson.est=temp.corr$estimate, Pearson.p=temp.corr$p.value,
                              Spearman.est=temp.corr.sp$estimate, Spearman.p=temp.corr.sp$p.value,
                              N_genes = length(topGenes$Genes), Genes = paste0(topGenes$Genes, collapse=", "), N=nrow(venAUC.wv))
      indVenPred3 <- rbind(indVenPred3, temp.sens)
    }
  }
  indVenPred3$Pearson.q <- qvalue::qvalue(indVenPred3$Pearson.p, pi0=1)$qvalue
  indVenPred3$Spearman.q <- qvalue::qvalue(indVenPred3$Spearman.p, pi0=1)$qvalue
  write.csv(indVenPred3, "venSensPrediction_sortedBead_1-100GenesByPosLog2FCUnweighted.csv", row.names = FALSE)
  
  # what about aza+ven
  indPred <- data.frame()
  for (d in "Azacytidine - Venetoclax") {
    venAUC <- na.omit(BeatAML$drug[,c("Barcode.ID",d)])
    venAUCprot <- BeatAML$global[BeatAML$global$Barcode.ID %in% venAUC$Barcode.ID,] # 9413 gene symbols
    # make sure each protein (columns 2+) are quantified in 2+ samples (rows)
    venAUCprot3 <- venAUCprot[,colSums(is.na(venAUCprot)) < (nrow(venAUCprot)-2)]
    
    indVenPred <- data.frame()
    for (i in temp.sig$Gene[temp.sig$Gene %in% colnames(venAUCprot3)[2:ncol(venAUCprot3)]]) { # 2504 gene symbols
      temp.sig2 <- temp.sig[temp.sig$Gene == i,]
      temp.wv <- DMEA::WV(venAUCprot3, temp.sig2)$scores
      venAUC.wv <- merge(venAUC, temp.wv, by="Barcode.ID")
      if (nrow(venAUC.wv) > 2) {
        temp.corr <- cor.test(venAUC.wv[[d]], venAUC.wv$WV)
        temp.corr.sp <- cor.test(venAUC.wv[[d]], venAUC.wv$WV, method="spearman")
        temp.sens <- data.frame(Drug=d, Pearson.est=temp.corr$estimate, Pearson.p=temp.corr$p.value,
                                Spearman.est=temp.corr.sp$estimate, Spearman.p = temp.corr.sp$p.value,
                                N_genes = 1, Genes = i, N=nrow(venAUC.wv))
        indVenPred <- rbind(indVenPred, temp.sens) 
      }
    }
    indPred <- rbind(indPred, indVenPred)
  }
  indPred$Pearson.q <- qvalue::qvalue(indPred$Pearson.p, pi0=1)$qvalue
  indPred$Spearman.q <- qvalue::qvalue(indPred$Spearman.p, pi0=1)$qvalue
  write.csv(indPred, "azaVenPrediction_sortedBead_1GeneAll.csv", row.names = FALSE) # PLBD1 but N=15: r=0.9102487,p=2.474523e-06,q=0.006181358; rho=0.9035714,p=0,q=0
  # NCF2: r=0.8183966, p=0.0001916759, q=0.02769824; rho=0.7607143,p=0.00151044,q=0.08174293 (rank 16 / 2498 based on Pearson.est)
  # 0.7586739 1.041759e-03 0.7642857 1.392934e-03 1 CD14 15 0.040035613 0.08174293 (rank 65 / 2498 based on Pearson.est)
  
  
  # what about other drugs
  indPred <- data.frame()
  for (d in colnames(BeatAML$drug)[2:ncol(BeatAML$drug)]) {
    venAUC <- na.omit(BeatAML$drug[,c("Barcode.ID",d)])
    venAUCprot <- BeatAML$global[BeatAML$global$Barcode.ID %in% venAUC$Barcode.ID,] # 9413 gene symbols
    
    indVenPred <- data.frame()
    for (i in temp.sig$Gene[temp.sig$Gene %in% colnames(venAUCprot)[2:ncol(venAUCprot)]]) { # 2504 gene symbols
      temp.sig2 <- temp.sig[temp.sig$Gene == i,]
      temp.wv <- DMEA::WV(venAUCprot, temp.sig2)$scores
      venAUC.wv <- merge(venAUC, temp.wv, by="Barcode.ID")
      if (nrow(venAUC.wv) > 2) {
        temp.corr <- cor.test(venAUC.wv[[d]], venAUC.wv$WV)
        temp.corr.sp <- cor.test(venAUC.wv[[d]], venAUC.wv$WV, method="spearman")
        temp.sens <- data.frame(Drug=d, Pearson.est=temp.corr$estimate, Pearson.p=temp.corr$p.value,
                                Spearman.est=temp.corr.sp$estimate, Spearman.p = temp.corr.sp$p.value,
                                N_genes = 1, Genes = i)
        indVenPred <- rbind(indVenPred, temp.sens) 
      }
    }
    indPred <- rbind(indPred, indVenPred)
  }
  write.csv(indPred, "drugSensPrediction_sortedBead_1GeneAll.csv", row.names = FALSE)
  
  return(list(loo = sensPred, ind=indVenPred))
  #return(list(loo = sensPred, venLoo = venPred, min=minSensPred, minSig=temp.sig2))
}

optMonoSig <- function(frac.meta, temp.sig, BeatAML, type="global", gmt) {
  genes <- temp.sig$Gene
  
  # try just ven so can look at KCTD12, S100A8
  fracProt <- BeatAML$global[BeatAML$global$Barcode.ID %in% frac.meta$Barcode.ID,] # 9413 gene symbols
  
  # calculate WV with genes
  indPred <- data.frame()
  testGenes <- genes[genes %in% colnames(fracProt)[2:ncol(fracProt)]]
  testedGenes <- unique(indPred$Genes)
  for (i in testGenes[!(testGenes %in% testedGenes)]) { # 2504 gene symbols
    temp.sig2 <- temp.sig[temp.sig$Gene == i,]
    temp.wv <- DMEA::WV(fracProt, temp.sig2)$scores
    venAUC.wv <- merge(frac.meta, temp.wv, by="Barcode.ID")
    if (nrow(venAUC.wv) > 2) {
      temp.corr <- cor.test(venAUC.wv$`%.Monocytes.in.PB`, venAUC.wv$WV)
      temp.corr.sp <- cor.test(venAUC.wv$`%.Monocytes.in.PB`, venAUC.wv$WV, method="spearman")
      temp.sens <- data.frame(Drug="%.Monocytes.in.PB", Pearson.est=temp.corr$estimate, Pearson.p=temp.corr$p.value,
                              Spearman.est=temp.corr.sp$estimate, Spearman.p = temp.corr.sp$p.value,
                              N_genes = 1, Genes = i, N=nrow(venAUC.wv))
      indPred <- rbind(indPred, temp.sens) 
    }
  }
  indPred$Pearson.q <- qvalue::qvalue(indPred$Pearson.p, pi0=1)$qvalue
  indPred$Spearman.q <- qvalue::qvalue(indPred$Spearman.p, pi0=1)$qvalue
  write.csv(indPred, "monoFracPrediction_sortedBead_1GeneAll.csv", row.names = FALSE) # CD93: pearson r 0.6326378, p 3.618788e-12, q 9.057827e-09; spearman rho 0.5979159, p 1.003341e-10, q 8.371211e-08
  # 0.496935085	2.25E-07	0.535276572	1.62E-08	1	CD14	97	1.41E-05	1.39E-06 (rank 41 by Pearson.est)
  # 0.231769655	0.022355499	0.07091542	0.49003613	1	CD34	97	0.045940734	0.55175908 (rank 1221)
  # 0.351859294	0.000408995	0.321062822	0.001343852	1	BCL2	97	0.00215519	0.004220402 (rank 477)
  # 0.522081296	4.15E-08	0.556290752	3.30E-09	1	NCF2	97	5.62E-06	5.58E-07 (rank 19)
  # 0.47261114	1.02E-06	0.521189999	4.42E-08	1	PLBD1	97	3.50E-05	2.76E-06 (rank 74)
  # note: % monocytes in PB correlation with Ven AUC: r=0.408, p=0.000897, N=63; Aza+Ven AUC: r=0.584, p = 0.0594, N=11 from: "sortedAMLproteomics_20241111.pptx" in PTRC2 folder
  
  return(indPred)
  #return(list(loo = sensPred, venLoo = venPred, min=minSensPred, minSig=temp.sig2))
}


optVenSigs <- function(sigs, BeatAML, predAll, types=c("global","rna","rna","rna"), 
                       weighted=c(TRUE, FALSE), weights=c("Log2FC","Pearson.est"),
                       selectionMetrics = c("Pearson.est","Pearson.p")) {
  venAUC <- BeatAML$drug[!is.na(BeatAML$drug$Venetoclax),c("Barcode.ID","Venetoclax")]
  
  allPred <- data.frame()
  for (s in 1:length(sigs)) {
    venAUCprot <- BeatAML[[types[s]]][BeatAML[[types[s]]]$Barcode.ID %in% venAUC$Barcode.ID,] # 9413 gene symbols
    
    temp.sig <- sigs[[s]]
    indVenPred <- predAll[predAll$Genes %in% temp.sig$Gene & predAll$DataType == types[s],]
    cat("Trying",names(sigs)[s],"signature\n")
    # indVenPred <- data.frame()
    # for (i in temp.sig$Gene[temp.sig$Gene %in% colnames(venAUCprot)[2:ncol(venAUCprot)]]) { # 2504 gene symbols
    #   temp.sig2 <- temp.sig[temp.sig$Gene == i,]
    #   
    #   temp.wv <- DMEA::WV(venAUCprot, temp.sig2)$scores
    #   venAUC.wv <- merge(venAUC, temp.wv, by="Barcode.ID")
    #   if (nrow(venAUC.wv)>2) {
    #     temp.corr <- cor.test(venAUC.wv$Venetoclax, venAUC.wv$WV)
    #     temp.corr.sp <- cor.test(venAUC.wv$Venetoclax, venAUC.wv$WV, method="spearman")
    #     temp.sens <- data.frame(Drug="Venetoclax", Pearson.est=temp.corr$estimate, Pearson.p=temp.corr$p.value,
    #                             Spearman.est=temp.corr.sp$estimate, Spearman.p=temp.corr.sp$p.value,
    #                             N_genes = 1, Genes = i, N=nrow(venAUC.wv))
    #     indVenPred <- rbind(indVenPred, temp.sens)
    #   }
    # }
    # indVenPred$Pearson.q <- qvalue::qvalue(indVenPred$Pearson.p, pi0=1)$qvalue
    # indVenPred$Spearman.q <- qvalue::qvalue(indVenPred$Spearman.p, pi0=1)$qvalue
    # indVenPred$SelectionMetric <- NA
    # indVenPred$Weighted <- TRUE
    # indVenPred$Signature <- names(sigs)[s]
    # allPred <- rbind(allPred, indVenPred)
    # 
    # indVenPred0 <- data.frame()
    # for (i in temp.sig$Gene[temp.sig$Gene %in% colnames(venAUCprot)[2:ncol(venAUCprot)]]) { # 2504 gene symbols
    #   temp.sig2 <- temp.sig[temp.sig$Gene == i,]
    #   temp.sig2$Log2FC <- 1
    #   temp.wv <- DMEA::WV(venAUCprot, temp.sig2)$scores
    #   venAUC.wv <- merge(venAUC, temp.wv, by="Barcode.ID")
    #   if (nrow(venAUC.wv)>2) {
    #     temp.corr <- cor.test(venAUC.wv$Venetoclax, venAUC.wv$WV)
    #     temp.corr.sp <- cor.test(venAUC.wv$Venetoclax, venAUC.wv$WV, method="spearman")
    #     temp.sens <- data.frame(Drug="Venetoclax", Pearson.est=temp.corr$estimate, Pearson.p=temp.corr$p.value,
    #                             Spearman.est=temp.corr.sp$estimate, Spearman.p=temp.corr.sp$p.value,
    #                             N_genes = 1, Genes = i, N=nrow(venAUC.wv))
    #     indVenPred0 <- rbind(indVenPred0, temp.sens)
    #   }
    # }
    # indVenPred0$Pearson.q <- qvalue::qvalue(indVenPred0$Pearson.p, pi0=1)$qvalue
    # indVenPred0$Spearman.q <- qvalue::qvalue(indVenPred0$Spearman.p, pi0=1)$qvalue
    # indVenPred2$SelectionMetric <- NA
    # indVenPred0$Weighted <- FALSE
    # indVenPred0$Weight <- NA
    # indVenPred0$Signature <- names(sigs)[s]
    # allPred <- rbind(allPred, indVenPred0)
    # 
    for (w in weighted) {
      cat("Weighted:",w,"\n")
      for (weight in weights) {
        cat("Weight:",weight,"\n")
        for (metric in selectionMetrics) {
          cat("Selection metric:",metric,"\n")
          indVenPred2 <- data.frame()
          for (i in seq(1,100)) {
            cat(i," ")
            indVenPred$rank <- indVenPred[,metric]
            if (metric == "Pearson.p") {
              topGenes <- indVenPred[indVenPred$Pearson.q <= 0.05,] %>% slice_min(rank, n=i)
            } else {
              topGenes <- indVenPred[indVenPred$Pearson.q <= 0.05,] %>% slice_max(rank, n=i)
            }
            temp.sig2 <- temp.sig[temp.sig$Gene %in% topGenes$Genes,]
            if (!w) {
              temp.sig2[,weight] <- 1 # give all genes equal weight 
            } else if (weight == "Pearson.est") {
              temp.sig2 <- topGenes[,c("Genes",weight)]
            }
            
            # make sure all genes are in expression
            temp.sig2 <- as.data.frame(dplyr::distinct(temp.sig2[temp.sig2[,1] %in% colnames(venAUCprot)[2:ncol(venAUCprot)],]))
            
            # make sure there are no NAs in expression
            venAUCprot2 <- as.data.frame(na.omit(venAUCprot[,c("Barcode.ID",unique(temp.sig2[,1]))]))
            
            
            temp.wv <- DMEA::WV(venAUCprot2, temp.sig2, weight.values=weight)$scores
            venAUC.wv <- merge(venAUC, temp.wv, by="Barcode.ID")
            if (nrow(venAUC.wv)>2) {
              temp.corr <- cor.test(venAUC.wv$Venetoclax, venAUC.wv$WV)
              temp.corr.sp <- cor.test(venAUC.wv$Venetoclax, venAUC.wv$WV, method="spearman")
              temp.sens <- data.frame(Drug="Venetoclax", Pearson.est=temp.corr$estimate, Pearson.p=temp.corr$p.value,
                                      Spearman.est=temp.corr.sp$estimate, Spearman.p=temp.corr.sp$p.value,
                                      N_genes = length(unique(temp.sig2[,1])), Genes = paste0(unique(temp.sig2[,1]), collapse=", "),
                                      Gene_directions=paste0(sign(temp.sig2[,2]), collapse=", "), allGenesSameDir = length(unique(sign(temp.sig2[,2])))==1,
                                      allGenes_direction = ifelse(length(unique(sign(temp.sig2[,2])))==1, unique(sign(temp.sig2[,2])), NA),
                                      N_genes_up = length(sign(temp.sig2[,2])[sign(temp.sig2[,2])==1]), 
                                      N_genes_dn = length(sign(temp.sig2[,2])[sign(temp.sig2[,2]) == -1]),
                                      N=nrow(venAUC.wv))
              indVenPred2 <- rbind(indVenPred2, temp.sens)
            }
          }
          indVenPred2$Pearson.q <- qvalue::qvalue(indVenPred2$Pearson.p, pi0=1)$qvalue
          indVenPred2$Spearman.q <- qvalue::qvalue(indVenPred2$Spearman.p, pi0=1)$qvalue
          indVenPred2$SelectionMetric <- metric
          indVenPred2$Weighted <- w
          indVenPred2$Weight <- ifelse(w, weight, NA)
          indVenPred2$Signature <- names(sigs)[s]
          allPred <- rbind(allPred, indVenPred2)
        }
      }
    }
    #   
    # 
    # 
    # indVenPred2 <- data.frame()
    # for (i in seq(1,100)) {
    #   topGenes <- indVenPred[indVenPred$Pearson.q <= 0.05,] %>% slice_min(Pearson.p, n=i)
    #   temp.sig2 <- temp.sig[temp.sig$Gene %in% topGenes$Genes,]
    #   temp.wv <- DMEA::WV(venAUCprot, temp.sig2)$scores
    #   venAUC.wv <- merge(venAUC, temp.wv, by="Barcode.ID")
    #   if (nrow(venAUC.wv)>2) {
    #     temp.corr <- cor.test(venAUC.wv$Venetoclax, venAUC.wv$WV)
    #     temp.corr.sp <- cor.test(venAUC.wv$Venetoclax, venAUC.wv$WV, method="spearman")
    #     temp.sens <- data.frame(Drug="Venetoclax", Pearson.est=temp.corr$estimate, Pearson.p=temp.corr$p.value,
    #                             Spearman.est=temp.corr.sp$estimate, Spearman.p=temp.corr.sp$p.value,
    #                             N_genes = length(topGenes$Genes), Genes = paste0(topGenes$Genes, collapse=", "), N=nrow(venAUC.wv))
    #     indVenPred2 <- rbind(indVenPred2, temp.sens)
    #   }
    # }
    # indVenPred2$Pearson.q <- qvalue::qvalue(indVenPred2$Pearson.p, pi0=1)$qvalue
    # indVenPred2$Spearman.q <- qvalue::qvalue(indVenPred2$Spearman.p, pi0=1)$qvalue
    # indVenPred2$SelectionMetric <- "Pearson.p"
    # indVenPred2$Weighted <- TRUE
    # indVenPred2$Weight <- "Log2FC"
    # indVenPred2$Signature <- names(sigs)[s]
    # allPred <- rbind(allPred, indVenPred2)
    # 
    # indVenPred3 <- data.frame()
    # for (i in seq(1,100)) {
    #   topGenes <- indVenPred[indVenPred$Pearson.q <= 0.05,] %>% slice_max(Pearson.est, n=i)
    #   temp.sig2 <- temp.sig[temp.sig$Gene %in% topGenes$Genes,]
    #   temp.sig2$Log2FC <- 1 # give all genes equal weight
    #   temp.wv <- DMEA::WV(venAUCprot, temp.sig2)$scores
    #   venAUC.wv <- merge(venAUC, temp.wv, by="Barcode.ID")
    #   if (nrow(venAUC.wv)>2) {
    #     temp.corr <- cor.test(venAUC.wv$Venetoclax, venAUC.wv$WV)
    #     temp.corr.sp <- cor.test(venAUC.wv$Venetoclax, venAUC.wv$WV, method="spearman")
    #     temp.sens <- data.frame(Drug="Venetoclax", Pearson.est=temp.corr$estimate, Pearson.p=temp.corr$p.value,
    #                             Spearman.est=temp.corr.sp$estimate, Spearman.p=temp.corr.sp$p.value,
    #                             N_genes = length(topGenes$Genes), Genes = paste0(topGenes$Genes, collapse=", "), N=nrow(venAUC.wv))
    #     indVenPred3 <- rbind(indVenPred3, temp.sens)
    #   }
    # }
    # indVenPred3$Pearson.q <- qvalue::qvalue(indVenPred3$Pearson.p, pi0=1)$qvalue
    # indVenPred3$Spearman.q <- qvalue::qvalue(indVenPred3$Spearman.p, pi0=1)$qvalue
    # indVenPred2$SelectionMetric <- "Pearson.est"
    # indVenPred3$Weighted <- FALSE
    # indVenPred3$Weight <- NA
    # indVenPred3$Signature <- names(sigs)[s]
    # allPred <- rbind(allPred, indVenPred3)
    # 
    # indVenPred3 <- data.frame()
    # for (i in seq(1,100)) {
    #   topGenes <- indVenPred[indVenPred$Pearson.q <= 0.05,] %>% slice_min(Pearson.p, n=i)
    #   temp.sig2 <- temp.sig[temp.sig$Gene %in% topGenes$Genes,]
    #   temp.sig2$Log2FC <- 1 # give all genes equal weight
    #   temp.wv <- DMEA::WV(venAUCprot, temp.sig2)$scores
    #   venAUC.wv <- merge(venAUC, temp.wv, by="Barcode.ID")
    #   if (nrow(venAUC.wv)>2) {
    #     temp.corr <- cor.test(venAUC.wv$Venetoclax, venAUC.wv$WV)
    #     temp.corr.sp <- cor.test(venAUC.wv$Venetoclax, venAUC.wv$WV, method="spearman")
    #     temp.sens <- data.frame(Drug="Venetoclax", Pearson.est=temp.corr$estimate, Pearson.p=temp.corr$p.value,
    #                             Spearman.est=temp.corr.sp$estimate, Spearman.p=temp.corr.sp$p.value,
    #                             N_genes = length(topGenes$Genes), Genes = paste0(topGenes$Genes, collapse=", "), N=nrow(venAUC.wv))
    #     indVenPred3 <- rbind(indVenPred3, temp.sens)
    #   }
    # }
    # indVenPred3$Pearson.q <- qvalue::qvalue(indVenPred3$Pearson.p, pi0=1)$qvalue
    # indVenPred3$Spearman.q <- qvalue::qvalue(indVenPred3$Spearman.p, pi0=1)$qvalue
    # indVenPred2$SelectionMetric <- "Pearson.p"
    # indVenPred3$Weighted <- FALSE
    # indVenPred3$Weight <- NA
    # indVenPred3$Signature <- names(sigs)[s]
    # allPred <- rbind(allPred, indVenPred3)
  }
  return(allPred)
}

optVen <- function(BeatAML, types=c("global","rna")) {
  venAUC <- BeatAML$drug[!is.na(BeatAML$drug$Venetoclax),c("Barcode.ID","Venetoclax")]
  
  allPred <- data.frame()
  for (t in types) {
    venAUCprot <- BeatAML[[t]][BeatAML[[t]]$Barcode.ID %in% venAUC$Barcode.ID,] # 9413 gene symbols
    
    # keep genes with 3+ unique non-NA values
    keepGenes <- names(venAUCprot)[sapply(venAUCprot, function(x) length(na.omit(unique(x)))) >= 3]
    venAUCprot <- venAUCprot[,keepGenes]
    
    indVenPred <- data.frame()
    for (i in colnames(venAUCprot)[2:ncol(venAUCprot)]) { # 2504 gene symbols
      temp.sig2 <- data.frame(Gene=i, Log2FC=1)
      
      temp.wv <- DMEA::WV(venAUCprot, temp.sig2)$scores
      venAUC.wv <- merge(venAUC, temp.wv, by="Barcode.ID")
      if (nrow(venAUC.wv)>2) {
        temp.corr <- cor.test(venAUC.wv$Venetoclax, venAUC.wv$WV)
        temp.corr.sp <- cor.test(venAUC.wv$Venetoclax, venAUC.wv$WV, method="spearman")
        temp.sens <- data.frame(Drug="Venetoclax", Pearson.est=temp.corr$estimate, Pearson.p=temp.corr$p.value,
                                Spearman.est=temp.corr.sp$estimate, Spearman.p=temp.corr.sp$p.value,
                                N_genes = 1, Genes = i, N=nrow(venAUC.wv))
        indVenPred <- rbind(indVenPred, temp.sens)
      }
    }
    indVenPred$Pearson.q <- qvalue::qvalue(indVenPred$Pearson.p, pi0=1)$qvalue
    indVenPred$Spearman.q <- qvalue::qvalue(indVenPred$Spearman.p, pi0=1)$qvalue
    indVenPred$Weighted <- FALSE
    indVenPred$Signature <- NA
    indVenPred$DataType <- t
    allPred <- rbind(allPred, indVenPred)
  }
  return(allPred)
}


#### run predictions with all genes ####
setwd("~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/")
setwd("data")

#gmt.drug <- readRDS("gmt_BeatAML_drug_MOA_2024-02-22.rds")
# drug.info <- read.csv("~/OneDrive - PNNL/Documents/PTRC2/BeatAML_single_drug_moa_2025-01-20.csv",
#                       stringsAsFactors=FALSE, fileEncoding="latin1")
# gmt.drug <- DMEA::as_gmt(drug.info, sep=", ")
# saveRDS(gmt.drug, "gmt_BeatAML_drug_MOA_2025-01-20.rds")
gmt.drug <- readRDS("gmt_BeatAML_drug_MOA_2025-01-20.rds")

# load sorted proteomics
# synapser::synLogin()
# global.df <- read.csv(synapser::synGet("syn58895933")$path) # DIA
# global.df <- global.df[ , which(colMeans(!is.na(global.df)) >= 0.75)] # 36 out of 48 samples are kept
# outliers <- c("X00839_CD34plusFlow", "X00117_CD34plus", 
#               "X00432_CD14plus", "X00251_CD14plus", "X00105_CD14plusFlow")
# # syn.test <- synapser::synGet("syn58914135") # for some reason, can't access pre-filtered version ???
# rownames(global.df) <- global.df$Gene
# global.df$Gene <- NULL
# global.df <- as.data.frame(t(global.df))
# global.df$Sample <- rownames(global.df)
# global.df <- global.df[,c("Sample", colnames(global.df)[1:(ncol(global.df)-1)])]
# global.df <- global.df[!(global.df$Sample %in% outliers),]
# global.df100 <- global.df[,colSums(is.na(global.df)) == 0]
# global.df100 <- global.df100[!grepl("flow",global.df100$Sample, ignore.case=TRUE),] # 17 samples

# load sorted proteomics signature
# sig.paths <- list("Sorted" = "analysis/combined24-27/DIA_2batches_noOutliers_noMSC/Sort Type_Bead/CD14_Pos_vs_Neg/global/Differential_expression/Differential_expression_results.csv",
#                   "van Galen" = "data/externalSignatures/formatted/notFilteredForMalignant/Differential_expression_van_Galen_AML_D0_Mono-like_vs_Prog-like_noNA.csv",
#                   "Triana" = "data/externalSignatures/formatted/Triana_RNA_AML_100PercentCells_Classical-Monocytes_vs_HSCs-and-MPPs_differentialExpression.csv",
#                   "Lasry" = "data/externalSignatures/formatted/notFilteredForMalignant/Differential_expression_Lasry_AML_CD14PosMonocyte_vs_HSC_protein-coding.csv")


# sig.paths <- list("Sorted" = "analysis/combined24-27/DIA_2batches_noOutliers_noMSC/Sort Type_Bead/CD14_Pos_vs_Neg/global/Differential_expression/Differential_expression_results.csv",
#                   "van Galen" = "data/externalSignatures/formatted/notFilteredForMalignant/Differential_expression_van_Galen_AML_D0_Mono-like_vs_Prog-like_protein-coding.csv",
#                   "Triana" = "data/externalSignatures/formatted/Triana_RNA_AML_100PercentCells_Classical-Monocytes_vs_HSCs-and-MPPs_differentialExpression_protein-coding.csv",
#                   "Lasry" = "data/externalSignatures/formatted/notFilteredForMalignant/Differential_expression_Lasry_AML_CD14PosMonocyte_vs_MPP_protein-coding.csv")
#cd14.sig <- na.omit(read.csv(synapser::synGet("syn64543462")$path)) # DIA

# sig.paths <- list("Sorted" = "analysis/combined24-27/DIA_2batches_noOutliers_noMSC/Sort Type_Bead/CD14_Pos_vs_Neg/global/Differential_expression/Differential_expression_results.csv",
#                   "van Galen" = "data/externalSignatures/formatted/notFilteredForMalignant/Differential_expression_van_Galen_AML_D0_Mono-like_vs_Prog-like.csv",
#                   "Triana" = "data/externalSignatures/formatted/Triana_RNA_AML_100PercentCells_Classical-Monocytes_vs_HSCs-and-MPPs_differentialExpression.csv",
#                   "Lasry" = "data/externalSignatures/formatted/notFilteredForMalignant/Differential_expression_Lasry_AML_CD14PosMonocyte_vs_MPP.csv")

sig.paths <- list("Sorted" = "analysis/combined24-27/using_cellType-sortType-patient_factors/DIA_2batches_noOutliers_noMSC_Bead/no_filter/CD14_Pos_vs_Neg/global/Differential_expression/Differential_expression_results.csv",
                  "van Galen" = "data/externalSignatures/formatted/notFilteredForMalignant/Differential_expression_van_Galen_AML_D0_Mono-like_vs_Prog-like.csv",
                  "Triana" = "data/externalSignatures/formatted/Triana_RNA_AML_100PercentCells_Classical-Monocytes_vs_HSCs-and-MPPs_differentialExpression.csv",
                  "Lasry" = "data/externalSignatures/formatted/notFilteredForMalignant/Differential_expression_Lasry_AML_CD14PosMonocyte_vs_MPP.csv")

setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/")
# prot.coding <- read.table("data/externalSignatures/uniprotkb_taxonomy_id_9606_AND_existenc_2025_07_14.tsv", header=TRUE, fill=TRUE)
# prot.coding <- prot.coding[prot.coding$Gene != "",]
# prot.coding.genes <- unique(unlist(strsplit(prot.coding$Gene," "))) # 14535
# sorted <- read.csv(sig.paths$Sorted)
# prot.coding.genes <- unique(c(prot.coding.genes, sorted$Gene)) # 21314
# 
# prot.coding <- read.table("data/externalSignatures/uniprotkb_taxonomy_id_9606_2025_07_14.tsv", header=TRUE, fill=TRUE, sep="\t")
# prot.coding <- prot.coding[prot.coding$Gene.Names != "",]
# prot.coding.genes <- unique(unlist(strsplit(prot.coding$Gene.Names," "))) # 126950
# sorted <- read.csv(sig.paths$Sorted)
# prot.coding.genes <- unique(c(prot.coding.genes, sorted$Gene)) # 127396

sorted.patients <- unique(dia.wo.out$meta$patient)
BeatAML <- load_not_norm_BeatAML_for_DMEA3(exclude.samples=sorted.patients)

uniprot.fasta <- seqinr::read.fasta("data/externalSignatures/UP000005640_9606.fasta.gz")
gene.info <- seqinr::getAnnot(uniprot.fasta)
gene.info2 <- sub(".*GN=","",gene.info)
gene.info2 <- sub(" .*","", gene.info2)
prot.coding.genes <- unique(gene.info2) # 20639
prot.coding.genes <- prot.coding.genes[!grepl("_HUMAN", prot.coding.genes)] # 20402
sorted <- read.csv(sig.paths$Sorted)
prot.coding.genes2 <- unique(c(prot.coding.genes, sorted$Gene)) # 20406
prot.coding.genes3 <- unique(c(prot.coding.genes2, colnames(BeatAML$global)[2:ncol(BeatAML$global)])) # 20707
sorted$Gene[!(sorted$Gene %in% prot.coding.genes)] # "STK19"    "SLC22A18" "TMEM199"  "CCDC115" were not in my uniprot download
length(sorted$Gene[(sorted$Gene %in% prot.coding.genes)]) / length(sorted$Gene) # 99.94192%
colnames(BeatAML$global)[2:ncol(BeatAML$global)][!(colnames(BeatAML$global)[2:ncol(BeatAML$global)] %in% prot.coding.genes)] # 305 were not in my uniprot download
length(colnames(BeatAML$global)[2:ncol(BeatAML$global)][(colnames(BeatAML$global)[2:ncol(BeatAML$global)] %in% prot.coding.genes)]) / length(colnames(BeatAML$global)[2:ncol(BeatAML$global)]) # 96.7598%
saveRDS(prot.coding.genes3, "proteinCodingGenes.rds")
setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/")
prot.coding.genes3 <- readRDS("proteinCodingGenes.rds")

global.df <- dia.wo.out$global
#rownames(global.df) <- global.df$Gene
#global.df$Gene <- NULL
global.df <- as.data.frame(t(global.df))
global.df$Sample <- rownames(global.df)
global.df <- global.df[,c("Sample", colnames(global.df)[1:(ncol(global.df)-1)])]
global.df100 <- global.df[,colSums(is.na(global.df)) == 0]
global.df100 <- global.df100[!grepl("_f_",global.df100$Sample, ignore.case=TRUE),] # 39 samples out of 51, 4417 proteins out of 6888

# import signatures and filter
sigs <- list()
setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/")
for (i in names(sig.paths)) {
  sigs[[i]] <- read.csv(sig.paths[[i]])
  sigs[[i]] <- sigs[[i]][sigs[[i]]$Gene %in% prot.coding.genes3,] # added 20250714
  sigs[[i]] <- na.omit(sigs[[i]][sigs[[i]]$adj.P.Val <= 0.05,c("Gene","Log2FC")])
}
saveRDS(sigs,"mono_vs_prog_sigs_ProteinCoding.rds")
sigs$`van Galen`$Gene[sigs$`van Galen`$Gene %in% sigs$Sorted$Gene] # "S100A9" "KIT"    "CEBPD"  "S100A8"
c("S100A9","KIT","CEBPD","S100A8")[c("S100A9","KIT","CEBPD","S100A8") %in% sigs$Triana$Gene] # "S100A9" "KIT"    "S100A8" are the 3 genes in all 4 sigs
venn.list <- list()
colorOrder <- c("Sorted","Lasry","Triana","van Galen")
for (i in colorOrder) {
  venn.list[[i]] <- unique(sigs[[i]]$Gene)
}
library(ggvenn)
fillVals = RColorBrewer::brewer.pal(length(sigs), "Set2")
ggvenn::ggvenn(venn.list, show_percentage = FALSE, fill_color=fillVals, set_name_size=5, text_size=5)
ggsave(paste0("mono_vs_prog_signature_vennDiagram_",Sys.Date(),".pdf"),width=7, height=7)

# run correlations between signatures
sig.df <- data.table::rbindlist(sigs, use.names=TRUE, idcol="Signature")
sig.df <- reshape2::dcast(sig.df, Gene ~ Signature, value.var="Log2FC")
corr <- data.frame()
for (i in colorOrder) {
  otherSigs <- colorOrder[colorOrder != i]
  temp.input <- sig.df[,c("Gene",i,otherSigs)]
  temp.corr <- DMEA::rank_corr(temp.input,plots=FALSE)$result
  temp.corr[nrow(temp.corr)+1,] <- c(i,1,rep(NA, ncol(temp.corr)-2))
  temp.corr$Signature <- i
  corr <- rbind(corr, temp.corr)
}
write.csv(corr, paste0("mono_vs_prog_signature_correlations_withSelfCorr_",Sys.Date(),".csv"), row.names=FALSE) # all significant
#maxAbsEst <- max(abs(corr$Pearson.est))
corr$Pearson.est <- as.numeric(corr$Pearson.est)
corr$Signature <- as.character(corr$Signature)
corr$Drug <- as.character(corr$Drug)
sigOrder <- corr[corr$Drug == "Sorted",]
sigOrder <- sigOrder[order(sigOrder$Pearson.est, decreasing=TRUE),]$Signature
corr$Signature <- factor(corr$Signature, levels=sigOrder)
corr$Drug <- factor(corr$Drug, levels=sigOrder)
ggplot(corr, aes(x=Drug, y=Signature, fill=Pearson.est)) + geom_tile() + 
  scale_fill_gradient2(limits=c(-1,1), low="blue", mid="grey", high="red")+labs(fill="Pearson r")+
  theme_minimal(base_size=16) + theme(axis.text=element_text(vjust=1, hjust=1, angle=45), axis.title=element_blank())
ggsave(paste0("mono_vs_prog_signature_correlations_withSelfCorr_orderBySortedR_",Sys.Date(),".pdf"), width=4, height=2.5)

# load cell fraction data
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
frac.meta <- na.omit(num.meta[,c("labId", "%.Monocytes.in.PB")]) # 625 samples
colnames(frac.meta)[1] <- "Barcode.ID"
#frac.meta$labId <- sub(".*-","X",frac.meta$labId) # don't need to exclude sorted patients here because they are already excluded from global data

setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/")
dir.create(paste0("Monocyte_vs_progenitor_signatures_beadOnly_",Sys.Date()))
setwd(paste0("Monocyte_vs_progenitor_signatures_beadOnly_",Sys.Date()))

# sorted.patients <- c("18-00105", "21-00839", "22-00571", "22-00117", "16-01184",
#                      "19-00074", "18-00103", "21-00432", "17-01060", "22-00251")
sorted.patients <- unique(dia.wo.out$meta$patient)
BeatAML <- load_not_norm_BeatAML_for_DMEA3(exclude.samples=sorted.patients)

# evaluate signature
evalResults <- compareSigs(global.df100, frac.meta, sigs, BeatAML = BeatAML, types=c("global", "rna", "rna", "rna"), gmt = gmt.drug) 
write.csv(evalResults$wv, "wv.csv", row.names = FALSE)
write.csv(evalResults$p, "pValues.csv", row.names = FALSE)
write.csv(evalResults$drug.corr, "drugCorrelations.csv", row.names = FALSE)
write.csv(evalResults$frac.corr, "cellFractionCorrelations.csv", row.names = FALSE)
saveRDS(evalResults$DMEA, "DMEA.rds")
saveRDS(evalResults$DMEA.Spearman, "DMEA_Spearman.rds")
all.DMEA.files <- list()
for (i in names(sigs)) {
  DMEA.files <- list("DMEA_WV_results.csv" =
                       evalResults$DMEA[[i]]$all.results[[1]]$WV.scores,
                     "DMEA_unused_weights.csv" =
                       evalResults$DMEA[[i]]$all.results[[1]]$unused.weights,
                     "DMEA_results.csv" =
                       evalResults$DMEA[[i]]$all.results[[1]]$result,
                     "DMEA_results_Spearman.csv" =
                       evalResults$DMEA.Spearman[[i]]$all.results[[1]]$result,
                     "DMEA_correlation_results.csv" = 
                       evalResults$DMEA[[i]]$all.results[[1]]$corr.result,
                     "DMEA_volcano_plot.pdf" =
                       evalResults$DMEA[[i]]$all.results[[1]]$volcano.plot,
                     "DMEA_volcano_plot_Spearman.pdf" =
                       evalResults$DMEA.Spearman[[i]]$all.results[[1]]$volcano.plot) 
  all.DMEA.files[[i]] <- DMEA.files
}
save_to_synapse_v2(all.DMEA.files#, "syn64606612"
                   )

# redo plots
setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/")
setwd("Monocyte_vs_progenitor_signatures_beadOnly_2025-07-14")
drug.corr.df <- read.csv("drugCorrelations.csv")
frac.corr.df <- read.csv("cellFractionCorrelations.csv")

fillVals = RColorBrewer::brewer.pal(length(sigs), "Set2")
rank.metrics <- c("Pearson.est", "Spearman.est")
for (i in rank.metrics) {
  descr <- stringr::str_split_1(i, "[.]")[1]
  if (descr == "Pearson") {
    ylab <- "Pearson r"
  } else {
    ylab <- "Spearman rho"
  }
  
  # drug sensitivity
  for (j in c("Aza + Ven", "Ven", "Aza")) {
    plot.df <- drug.corr.df[drug.corr.df$`Drug.Treatment` == j,]
    plot.df$rank <- plot.df[,i]
    sigOrder <- na.omit(unique(plot.df[order(plot.df$rank, decreasing=TRUE),]$Signature))
    # ggplot(plot.df, aes(x=Signature, y=rank, fill = Signature, alpha=0.5)) + 
    #   geom_col() + theme_classic(base_size = 12) + 
    #   ylab(ylab) + 
    #   ggplot2::scale_x_discrete(limits = sigOrder) +
    #   scale_fill_manual(values=fillVals, 
    #                     breaks=c("Sorted","Lasry","Triana","van Galen"))+
    #   ggtitle(paste("Monocytic signatures predict", j, "sensitivity"))
    # ggsave(paste0(j,"_Corr_DIA_WV_signatureFill_barPlot_",descr,"_", Sys.Date(), ".pdf"), width = 5, height = 3)
    
    ggplot(plot.df, aes(x=Signature, y=rank, fill = Signature)) + 
      geom_col(alpha=0.5, show.legend=FALSE) + theme_classic(base_size = 12) + 
      ylab(ylab) + theme(axis.title.x=element_blank(), 
                         axis.title.y=element_text(size=16),
                         axis.text.x=element_text(size=16, angle=45, vjust=1,hjust=1)) +
      ggplot2::scale_x_discrete(limits = sigOrder) +
      scale_fill_manual(values=fillVals, 
                        breaks=c("Sorted","Lasry","Triana","van Galen"))#+
    #ggtitle(paste("Monocytic signatures predict", j, "sensitivity"))
    ggsave(paste0(j,"_Corr_DIA_WV_signatureFill_barPlot_",descr,"_", Sys.Date(), ".pdf"), width = 2, height = 2)
  }
  
  # mono fraction
  plot.df <- frac.corr.df
  plot.df$rank <- plot.df[,i]
  sigOrder <- plot.df[order(plot.df$rank, decreasing=TRUE),]$Signature
  ggplot(plot.df, aes(x=Signature, y=rank, fill = Signature)) + 
    geom_col(alpha=0.5, show.legend=FALSE) + theme_classic(base_size = 12) + ylab(ylab) + 
    theme(axis.title.x=element_blank(), 
          axis.title.y=element_text(size=16),
          axis.text.x=element_text(size=16, angle=45, vjust=1,hjust=1)) +
    ggplot2::scale_x_discrete(limits = sigOrder) +
    scale_fill_manual(values=fillVals, 
                      breaks=c("Sorted","Lasry","Triana","van Galen"))#+
    #ggtitle("Monocytic signatures predict monocyte fraction")
  ggsave(paste0("fracCorr_DIA_WV_",descr,"_signatureFill_", Sys.Date(),".pdf"), width = 2, height = 2)
}

#### find min number of important genes ####
setwd("~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/")
setwd("data")

gmt.drug <- readRDS("gmt_BeatAML_drug_MOA_2025-01-20.rds")

# load sorted proteomics signature
# sig.paths <- list("Sorted" = "analysis/combined24-27/DIA_2batches_noOutliers_noMSC/Sort Type_Bead/CD14_Pos_vs_Neg/global/Differential_expression/Differential_expression_results.csv",
#                   "van Galen" = "data/externalSignatures/formatted/notFilteredForMalignant/Differential_expression_van_Galen_AML_D0_Mono-like_vs_Prog-like_noNA.csv",
#                   "Triana" = "data/externalSignatures/formatted/Triana_RNA_AML_100PercentCells_Classical-Monocytes_vs_HSCs-and-MPPs_differentialExpression.csv",
#                   "Lasry" = "data/externalSignatures/formatted/notFilteredForMalignant/Differential_expression_Lasry_AML_CD14PosMonocyte_vs_HSC_protein-coding.csv")

# sig.paths <- list("Sorted" = "analysis/combined24-27/DIA_2batches_noOutliers_noMSC/Sort Type_Bead/CD14_Pos_vs_Neg/global/Differential_expression/Differential_expression_results.csv",
#                   "van Galen" = "data/externalSignatures/formatted/notFilteredForMalignant/Differential_expression_van_Galen_AML_D0_Mono-like_vs_Prog-like_protein-coding.csv",
#                   "Triana" = "data/externalSignatures/formatted/Triana_RNA_AML_100PercentCells_Classical-Monocytes_vs_HSCs-and-MPPs_differentialExpression_protein-coding.csv",
#                   "Lasry" = "data/externalSignatures/formatted/notFilteredForMalignant/Differential_expression_Lasry_AML_CD14PosMonocyte_vs_MPP_protein-coding.csv")


dia.wo.out <- readRDS("~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/analysis/DIA_2batches_noOutliers.rds")

global.df <- dia.wo.out$global
#rownames(global.df) <- global.df$Gene
#global.df$Gene <- NULL
global.df <- as.data.frame(t(global.df))
global.df$Sample <- rownames(global.df)
global.df <- global.df[,c("Sample", colnames(global.df)[1:(ncol(global.df)-1)])]
global.df100 <- global.df[,colSums(is.na(global.df)) == 0]
global.df100 <- global.df100[!grepl("_f_",global.df100$Sample, ignore.case=TRUE),] # 39 samples out of 51, 4417 proteins out of 6888

sorted.patients <- unique(dia.wo.out$meta$patient)
BeatAML <- load_not_norm_BeatAML_for_DMEA3(exclude.samples=sorted.patients)

# import signatures and filter
sigs <- list()
setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/")
for (i in names(sig.paths)) {
  sigs[[i]] <- read.csv(sig.paths[[i]])
  sigs[[i]] <- na.omit(sigs[[i]][sigs[[i]]$adj.P.Val <= 0.05,c("Gene","Log2FC")])
}

# 
# # optimize sorted sig
# optSorted <- optSig(global.df100, sigs[["Sorted"]], BeatAML=BeatAML, type="global", gmt=gmt.drug)
# #write.csv(optSorted$loo, "drugSensPrediction_sortedBead_geneLOO.csv", row.names = FALSE)
# #write.csv(optSorted$venLoo, "venSensPrediction_sortedBead_geneLOO.csv", row.names = FALSE)
# #write.csv(optSorted$min, "drugSensPrediction_sortedBead_minGenes.csv", row.names = FALSE)
# #write.csv(optSorted$minSig, "drugSensPrediction_sortedBead_minGeneSignature.csv", row.names = FALSE)
# 
# venn.list <- list()
# colorOrder <- c("Sorted","Lasry","Triana","van Galen")
# for (i in colorOrder) {
#   venn.list[[i]] <- unique(sigs[[i]]$Gene)
# }
# library(ggvenn)
# 
# ggvenn::ggvenn(venn.list, show_percentage = FALSE, fill_color=fillVals, set_name_size=5, text_size=5)
# ggsave("mono_vs_prog_signature_vennDiagram.pdf",width=7, height=7)
# 
# # run correlations between signatures
# sig.df <- data.table::rbindlist(sigs, use.names=TRUE, idcol="Signature")
# sig.df <- reshape2::dcast(sig.df, Gene ~ Signature, value.var="Log2FC")
# corr <- data.frame()
# for (i in colorOrder) {
#   otherSigs <- colorOrder[colorOrder != i]
#   temp.input <- sig.df[,c("Gene",i,otherSigs)]
#   temp.corr <- DMEA::rank_corr(temp.input,plots=FALSE)$result
#   temp.corr[nrow(temp.corr)+1,] <- c(i,1,rep(NA, ncol(temp.corr)-2))
#   temp.corr$Signature <- i
#   corr <- rbind(corr, temp.corr)
# }
# write.csv(corr, "mono_vs_prog_signature_correlations.csv", row.names=FALSE) # all significant
# write.csv(corr, "mono_vs_prog_signature_correlations_withSelfCorr.csv", row.names=FALSE) # all significant
# #maxAbsEst <- max(abs(corr$Pearson.est))
# corr$Pearson.est <- as.numeric(corr$Pearson.est)
# ggplot(corr, aes(x=Drug, y=Signature, fill=Pearson.est)) + geom_tile() + 
#   scale_fill_gradient2(limits=c(-1,1), low="blue", mid="grey", high="red")+labs(fill="Pearson r")+
#   theme_minimal(base_size=16) + theme(axis.text=element_text(vjust=1, hjust=1, angle=45), axis.title=element_blank())
# ggsave("mono_vs_prog_signature_correlations_withSelfCorr.pdf", width=4, height=2.5)
# 
# setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/")
# dir.create("Monocyte_vs_progenitor_signatures_beadOnly_2025-05-30")
# setwd("Monocyte_vs_progenitor_signatures_beadOnly_2025-05-30")
# 
# sorted.patients <- unique(dia.wo.out$meta$patient)
# BeatAML <- load_not_norm_BeatAML_for_DMEA3(exclude.samples=sorted.patients)
# 

# redo plots
setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/")
setwd("Monocyte_vs_progenitor_signatures_beadOnly_2025-07-14")
predAll <- optVen(BeatAML)
predAll <- read.csv("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/Monocyte_vs_progenitor_signatures_beadOnly_2025-07-09/venSensPredictions_noSignature_2025-07-09.csv")
predAll$ProteinCoding <- FALSE
# txs <- transcripts(edb, filter=GeneNameFilter(unique(predAll[predAll$N_genes==1,]$Genes)), columns = "tx_biotype")
# protein.coding.genes <- txs[txs$tx_biotype == "protein_coding",]$gene_name
predAll[predAll$Genes %in% prot.coding.genes3,]$ProteinCoding <- TRUE # 74910
write.csv(predAll, "venSensPredictions_noSignature_2025-07-14.csv", row.names=FALSE)
# best are:
# Venetoclax 0.8133223 3.292705e-26 0.8022521 0 11 LRRC25, HMOX1, LRP1, SLC15A3, LILRB2, LILRA6, COTL1, CHST15, RBM47, FCGRT, SGSH 106 2.750886e-25 0 FALSE Lasry

pred <- optVenSigs(sigs, BeatAML, predAll) # previously selected by Pearson.est
pred$DataType <- "rna"
pred[pred$Signature=="Sorted",]$DataType <- "global"
pred$ProteinCoding <- TRUE # was false before
# library(ensembldb)
# BiocManager::install("EnsDb.Hsapiens.v86")
# library(EnsDb.Hsapiens.v86)
# edb <- EnsDb.Hsapiens.v86
# ## Evaluate whether we have protein annotation available
# hasProteinData(edb)
# listTables(edb)
# txs <- transcripts(edb, filter=GeneNameFilter(unique(pred[pred$N_genes==1,]$Genes)), columns = "tx_biotype")
# protein.coding.genes <- txs[txs$tx_biotype == "protein_coding",]$gene_name # 17704

# # split gene vectors for each row
# geneList <- strsplit(pred$Genes, ", ")
# pred[all(geneList %in% protein.coding.genes),]$ProteinCoding <- TRUE # 339

write.csv(pred, "venSensPredictions_2025-09-04.csv", row.names=FALSE)
pred <- read.csv("venSensPredictions_2025-09-04.csv")
#topPred <- pred[which.max(pred$Pearson.est) | which.max(-log10(pred$Pearson.p)),]
# best are all non-weighted:
# RNA: Lasry: SLC7A7, LILRA6, LRP1, SGSH, LILRB1, CLEC7A, RBM47, FGR, LRRC25, SLC15A3, LILRB2, IQSEC1, CHST15, TNFRSF1B, HMOX1, CD1D (16, r=0.811)
# Protein: Sorted: NCF2, FCGRT, KCTD12, CD93 (4, r=0.795)
# RNA: Lasry: LRRC25 (1, r=0.760)
# Protein: Sorted: NCF2 (1, r=0.728)

# ggplot(pred[pred$Pearson.q <= 0.05,], aes(x=N_genes, y=Pearson.est, color=Signature, shape=Weighted, alpha=N)) +
#   geom_point() + geom_smooth(se=FALSE, linetype="dashed") + theme_classic() + scale_x_continuous(transform = "log10") +
#   geom_hline(yintercept =0, linetype="dashed", color="gray")+labs(x="# of Genes", y="Pearson Correlation Estimate") + 
#   ggrepel::geom_label_repel(data=rbind(pred[pred$Signature=="Sorted" & pred$Pearson.q<=0.05,] %>% slice_max(Pearson.est, n=1, with_ties=FALSE),
#                                        pred[pred$Signature=="Lasry" & pred$Pearson.q<=0.05,] %>% slice_max(Pearson.est, n=1, with_ties=FALSE)),
#                             aes(label=Genes))
# ggsave("Ngenes_allSigs_PearsonEst_wLabelandNalpha.pdf", width=5, height=3)
# 
# ggplot(pred[pred$N_genes==1 & pred$Pearson.q<=0.05,], 
#        aes(x=Pearson.est, y=-log10(Pearson.p), color=Signature, shape=Weighted, alpha=N)) +
#   geom_point() + theme_classic() +
#   geom_hline(yintercept =-log10(0.05), linetype="dashed", color="gray")+
#   labs(y="-Log(P-value)", x="Pearson Correlation Estimate") + 
#   ggrepel::geom_label_repel(data=rbind(pred[pred$N_genes==1 & pred$Pearson.q<=0.05,] %>% slice_min(Pearson.p, n=1, with_ties=FALSE),
#                                        pred[pred$N_genes==1 & pred$Pearson.q<=0.05,] %>% slice_max(Pearson.est, n=1, with_ties=FALSE)),
#                             aes(label=Genes))
# ggsave("1gene_allSigs_PearsonEst_withNalpha.pdf", width=4, height=3)

#singlePred <- rbind(pred[pred$N_genes==1 & !pred$Weighted,colnames(predAll)], predAll)
# singlePred2 <- plyr::ddply(singlePred, .(Genes,DataType,Pearson.est,Pearson.p,N), summarize,
#                            Signature=paste0(na.omit(unique(Signature)), collapse=", "),
#                            ProteinCoding = any(ProteinCoding))
# there are still genes in a signature also included as not being in a signature?
# also 301 global genes are labeled as not protein-coding
singlePred2 <- predAll
for (i in names(sigs)) {
  if (any(!is.na(singlePred2$Signature))) {
    singlePred2[singlePred2$Genes %in% sigs[[i]]$Gene & !is.na(singlePred2$Signature),]$Signature <- 
      paste0(singlePred2[singlePred2$Genes %in% sigs[[i]]$Gene & !is.na(singlePred2$Signature),]$Signature, ", ", i)
  }
  singlePred2[singlePred2$Genes %in% sigs[[i]]$Gene & is.na(singlePred2$Signature),]$Signature <- i
}
singlePred2$Significant <- FALSE
singlePred2[singlePred2$Pearson.p <= 0.05,]$Significant <- TRUE
singlePred2[is.na(singlePred2$Signature),]$Signature <- "None"
#singlePred2[singlePred2$Signature=="",]$Signature <- "None"
ggplot(singlePred2, # was 32,254 rows on 2025-07-09, now 35,870 rows on 2025-07-14
       aes(x=Pearson.est, y=-log10(Pearson.p), color=Signature, shape=DataType, alpha=N)) +
  geom_point() + theme_classic() + scale_color_manual(values=c("gray",scales::hue_pal()(length(unique(singlePred2$Signature))-1)),
                                                      breaks=c("None", unique(singlePred2[order(singlePred2$Pearson.p),]$Signature)[unique(singlePred2[order(singlePred2$Pearson.p),]$Signature) != "None"]))+
  geom_hline(yintercept =-log10(0.05), linetype="dashed", color="gray")+
  labs(y="-Log(P-value)", x="Pearson Correlation Estimate",shape="Data Type") + 
  scale_shape_manual(values=c(16,17),labels=c("Protein","RNA"))+#scale_alpha_continuous()+
  ggrepel::geom_label_repel(data=dplyr::distinct(rbind(singlePred2[singlePred2$DataType=="global",] %>% slice_min(Pearson.p, n=1, with_ties=FALSE),
                                       singlePred2[singlePred2$DataType=="global" & singlePred2$Signature!="",] %>% slice_min(Pearson.p, n=1, with_ties=FALSE),
                                       singlePred2[singlePred2$DataType=="rna",] %>% slice_min(Pearson.p, n=1, with_ties=FALSE),
                                       singlePred2[singlePred2$DataType=="rna" & singlePred2$Signature!="",] %>% slice_min(Pearson.p, n=1, with_ties=FALSE))),
                            aes(label=Genes))
ggsave(paste0("1gene_allSigs_PearsonEst_withNalpha_v2_",Sys.Date(),".pdf"), width=4, height=3)

singlePred3 <- predAll[predAll$Pearson.p<0.05 & predAll$N>=100,]
for (i in names(sigs)) {
  singlePred3[,i] <- FALSE
  singlePred3[singlePred3$Genes %in% sigs[[i]]$Gene,i] <- TRUE
}
singlePred3$Significant <- FALSE
singlePred3[singlePred3$Pearson.p <= 0.05,]$Significant <- TRUE
singlePred3$`None (Protein)` <- TRUE
singlePred3[singlePred3$DataType=="global" & (singlePred3$Sorted | singlePred3$Triana | singlePred3$Lasry | singlePred3$`van Galen`),]$`None (Protein)` <- FALSE
singlePred3$`None (RNA)` <- TRUE
singlePred3[singlePred3$DataType=="rna" & (singlePred3$Sorted | singlePred3$Triana | singlePred3$Lasry | singlePred3$`van Galen`),]$`None (RNA)` <- FALSE
long3 <- reshape2::melt(singlePred3[singlePred3$ProteinCoding,c("Genes","Pearson.est",names(sigs),"None (Protein)","None (RNA)")], 
                        id.vars = c("Genes", "Pearson.est"), variable.name="Signature",measured.vars=c(names(sigs),"None (Protein)","None (RNA)"))
long3 <- long3[long3$value,]

ggplot(na.omit(long3), 
       aes(x=reorder(Signature, -abs(Pearson.est), FUN=median), y=abs(Pearson.est), color=Signature)) + geom_violin(alpha=0) + 
  geom_point() + geom_boxplot(width=0.4, alpha = 0) + theme_classic() + #scale_color_manual(values=c("gray",scales::hue_pal()(length(unique(singlePred2$Signature))-1)))+
  scale_color_manual(values=c(RColorBrewer::brewer.pal(length(unique(pred$Signature)), "Set2"),"darkgrey","grey"),
                     breaks=c(names(sigs),"None (Protein)", "None (RNA)"))+
  labs(y="Absolute Pearson r", x="Signature") + theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))
ggsave(paste0("1gene_proteinCoding_allSigs_PearsonEst_boxPlot_minN100_",Sys.Date(),".pdf"), width=4, height=4)
ggsave(paste0("1gene_proteinCoding_allSigs_PearsonEst_boxPlot_minN100_width6_",Sys.Date(),".pdf"), width=6, height=4)
ggsave(paste0("1gene_proteinCoding_allSigs_PearsonEst_boxPlot_minN100_width5_",Sys.Date(),".pdf"), width=5, height=4)

singlePred3 <- predAll[predAll$Pearson.p<0.05,]
for (i in names(sigs)) {
  singlePred3[,i] <- FALSE
  singlePred3[singlePred3$Genes %in% sigs[[i]]$Gene,i] <- TRUE
}
singlePred3$Significant <- FALSE
singlePred3[singlePred3$Pearson.p <= 0.05,]$Significant <- TRUE
singlePred3$`None (Protein)` <- TRUE
singlePred3[singlePred3$DataType=="global" & (singlePred3$Sorted | singlePred3$Triana | singlePred3$Lasry | singlePred3$`van Galen`),]$`None (Protein)` <- FALSE
singlePred3$`None (RNA)` <- TRUE
singlePred3[singlePred3$DataType=="rna" & (singlePred3$Sorted | singlePred3$Triana | singlePred3$Lasry | singlePred3$`van Galen`),]$`None (RNA)` <- FALSE
long3 <- reshape2::melt(singlePred3[singlePred3$ProteinCoding,c("Genes","Pearson.p",names(sigs),"None (Protein)","None (RNA)")], 
                        id.vars = c("Genes", "Pearson.p"), variable.name="Signature",measured.vars=c(names(sigs),"None (Protein)","None (RNA)"))
long3 <- long3[long3$value,]

ggplot(na.omit(long3), 
       aes(x=reorder(Signature, -log10(Pearson.p), FUN=median), y=-log10(Pearson.p), color=Signature)) + geom_violin(alpha=0) + 
  geom_point() + geom_boxplot(width=0.4, alpha = 0) + theme_classic() + #scale_color_manual(values=c("gray",scales::hue_pal()(length(unique(singlePred2$Signature))-1)))+
  scale_color_manual(values=c(RColorBrewer::brewer.pal(length(unique(pred$Signature)), "Set2"),"darkgrey","grey"),
                     breaks=c(names(sigs),"None (Protein)", "None (RNA)"))+
  labs(y="-Log(Pearson p-value)", x="Signature") + theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))
ggsave(paste0("1gene_proteinCoding_allSigs_PearsonP_boxPlot_minN100_",Sys.Date(),".pdf"), width=4, height=4)
ggsave(paste0("1gene_proteinCoding_allSigs_PearsonP_boxPlot_minN100_width6_",Sys.Date(),".pdf"), width=6, height=4)
ggsave(paste0("1gene_proteinCoding_allSigs_PearsonP_boxPlot_minN100_width5_",Sys.Date(),".pdf"), width=5, height=4)

ggplot(singlePred2[singlePred2$ProteinCoding,], # 23,781 (73.73%); note: 9153/9411 (97.26%) global rows are protein-coding? maybe because of difference in database?
       aes(x=Pearson.est, y=-log10(Pearson.p), color=Signature, shape=DataType, alpha=N)) +
  geom_point() + theme_classic() + #scale_color_manual(values=c("gray",scales::hue_pal()(length(unique(singlePred2$Signature))-1)))+
  scale_color_manual(values=c("gray",scales::hue_pal()(length(unique(singlePred2$Signature))-1)),
                     breaks=c("None", unique(singlePred2[order(singlePred2$Pearson.p),]$Signature)[unique(singlePred2[order(singlePred2$Pearson.p),]$Signature) != "None"]))+
  geom_hline(yintercept =-log10(0.05), linetype="dashed", color="gray")+
  labs(y="-Log(P-value)", x="Pearson Correlation Estimate",shape="Data Type") + 
  scale_shape_manual(values=c(16,17),labels=c("Protein","RNA"))+#scale_alpha_continuous()+
  ggrepel::geom_label_repel(data=dplyr::distinct(rbind(singlePred2[singlePred2$DataType=="global" & singlePred2$ProteinCoding,] %>% slice_min(Pearson.p, n=1, with_ties=FALSE),
                                                       singlePred2[singlePred2$DataType=="global" & singlePred2$Signature!="None" & singlePred2$ProteinCoding,] %>% slice_min(Pearson.p, n=1, with_ties=FALSE),
                                                       singlePred2[singlePred2$DataType=="rna" & singlePred2$ProteinCoding,] %>% slice_min(Pearson.p, n=1, with_ties=FALSE),
                                                       singlePred2[singlePred2$DataType=="rna" & singlePred2$Signature!="None" & singlePred2$ProteinCoding,] %>% slice_min(Pearson.p, n=1, with_ties=FALSE))),
                            aes(label=Genes), show_guides=FALSE) + theme(legend.position="top") + guides(color=guide_legend(position="right"), alpha=guide_legend(position="bottom"))
ggsave(paste0("1gene_proteinCoding_allSigs_PearsonEst_withNalpha_v3_",Sys.Date(),".pdf"), width=4, height=4)
ggsave(paste0("1gene_proteinCoding_allSigs_PearsonEst_withNalpha_v3_wider_",Sys.Date(),".pdf"), width=5, height=4)

ggplot(singlePred2[singlePred2$ProteinCoding,], # 23,781 (73.73%); note: 9153/9411 (97.26%) global rows are protein-coding? maybe because of difference in database?
       aes(x=Pearson.est, y=-log10(Pearson.p), color=Signature, shape=DataType, alpha=N)) +
  geom_point() + theme_classic() + #scale_color_manual(values=c("gray",scales::hue_pal()(length(unique(singlePred2$Signature))-1)))+
  scale_color_manual(values=c("gray",scales::hue_pal()(length(unique(singlePred2$Signature))-1)),
                     breaks=c("None", unique(singlePred2[order(singlePred2$Pearson.p),]$Signature)[unique(singlePred2[order(singlePred2$Pearson.p),]$Signature) != "None"]))+
  geom_hline(yintercept =-log10(0.05), linetype="dashed", color="gray")+
  labs(y="-Log(P-value)", x="Pearson Correlation Estimate",shape="Data Type") + 
  scale_shape_manual(values=c(16,17),labels=c("Protein","RNA"))+#scale_alpha_continuous()+
  ggrepel::geom_label_repel(data=dplyr::distinct(rbind(singlePred2[singlePred2$DataType=="global" & singlePred2$ProteinCoding,] %>% slice_min(Pearson.p, n=1, with_ties=FALSE),
                                                       singlePred2[singlePred2$DataType=="global" & singlePred2$Signature!="None" & singlePred2$ProteinCoding,] %>% slice_min(Pearson.p, n=1, with_ties=FALSE),
                                                       singlePred2[singlePred2$DataType=="rna" & singlePred2$ProteinCoding,] %>% slice_min(Pearson.p, n=1, with_ties=FALSE),
                                                       singlePred2[singlePred2$DataType=="rna" & singlePred2$Signature!="None" & singlePred2$ProteinCoding,] %>% slice_min(Pearson.p, n=1, with_ties=FALSE))),
                            aes(label=Genes), show_guides=TRUE) + theme(legend.position="top") + guides(color=guide_legend(position="right"), alpha=guide_legend(position="bottom"))
ggsave(paste0("1gene_proteinCoding_allSigs_PearsonEst_withNalpha_",Sys.Date(),".pdf"), width=4, height=4)
ggsave(paste0("1gene_proteinCoding_allSigs_PearsonEst_withNalpha_wider_",Sys.Date(),".pdf"), width=5, height=4)

# do the protein-coding genes which belong to a signature have higher Pearson est or lower p than those which don't belong to a signature?
p.test <- t.test(singlePred2[singlePred2$ProteinCoding & singlePred2$Signature!="None",]$Pearson.p, # 5572 with mean 0.130
                     singlePred2[singlePred2$ProteinCoding & singlePred2$Signature=="None",]$Pearson.p, # 18334 with mean 0.226
                   "less")
p.test
# yes: p=6.458747e-131

est.test <- t.test(singlePred2[singlePred2$ProteinCoding & singlePred2$Signature!="None",]$Pearson.est, # 5572 with mean 0.0144
                   singlePred2[singlePred2$ProteinCoding & singlePred2$Signature=="None",]$Pearson.est, # 18334 with mean -0.0263
                   "greater")
# yes: p=1.493237e-16

est.test <- t.test(abs(singlePred2[singlePred2$ProteinCoding & singlePred2$Signature!="None",]$Pearson.est), # 5572 with mean 0.296
                   abs(singlePred2[singlePred2$ProteinCoding & singlePred2$Signature=="None",]$Pearson.est), # 18334 with mean 0.206
                   "greater")
# yes: p=6.533182e-250

# do proteins have more significant Pearson correlations than genes?
# sharedGenes <- plyr::ddply(singlePred2, .(Genes), summarize,
#                            protP = Pearson.p[DataType=="global"],
#                            rnaP = Pearson.p[DataType=="rna"])
sharedGenes <- reshape2::dcast(singlePred2, Genes~DataType, value.var="Pearson.p")
sharedGenes <- na.omit(sharedGenes) # 8861 gene symbols in both rna and global data
p.omics.test <- t.test(sharedGenes$global, sharedGenes$rna, "less", paired=TRUE)
# no, p = 1

# do proteins have less significant Pearson correlations than genes?
p.omics.test <- t.test(sharedGenes$global, sharedGenes$rna, "greater", paired=TRUE)
# yes, p = 4.58224e-20

p.omics.test <- t.test(singlePred2[singlePred2$ProteinCoding & singlePred2$DataType=="global",]$Pearson.p, 
                       singlePred2[singlePred2$ProteinCoding & singlePred2$DataType=="rna",]$Pearson.p, 
                       "less", paired=FALSE)
# no, p = 1

p.omics.test <- t.test(singlePred2[singlePred2$ProteinCoding & singlePred2$DataType=="global",]$Pearson.p, 
                       singlePred2[singlePred2$ProteinCoding & singlePred2$DataType=="rna",]$Pearson.p, 
                       "greater", paired=FALSE)
# yes, p = 3.062636e-11

sharedGenes <- reshape2::dcast(singlePred2, Genes~DataType, value.var="Pearson.est")
sharedGenes <- na.omit(sharedGenes) # 8861 gene symbols in both rna and global data
p.omics.test <- t.test(sharedGenes$global, sharedGenes$rna, "less", paired=TRUE)
p.omics.test$p.value
# no, p = 1

# do proteins have less significant Pearson correlations than genes?
p.omics.test <- t.test(sharedGenes$global, sharedGenes$rna, "greater", paired=TRUE)
p.omics.test$p.value
# yes, p = 8.946384e-72

p.omics.test <- t.test(singlePred2[singlePred2$ProteinCoding & singlePred2$DataType=="global",]$Pearson.est, 
                       singlePred2[singlePred2$ProteinCoding & singlePred2$DataType=="rna",]$Pearson.est, 
                       "less", paired=FALSE)
p.omics.test$p.value
# no, p = 1

p.omics.test <- t.test(singlePred2[singlePred2$ProteinCoding & singlePred2$DataType=="global",]$Pearson.est, 
                       singlePred2[singlePred2$ProteinCoding & singlePred2$DataType=="rna",]$Pearson.est, 
                       "greater", paired=FALSE)
p.omics.test$p.value
# yes, p = 3.401272e-44

fillVals = RColorBrewer::brewer.pal(length(unique(pred$Signature)), "Set2")
p1 <- ggplot(dplyr::distinct(pred[pred$SelectionMetric=="Pearson.p" & !pred$allGenesSameDir & pred$Weighted,]), aes(x=N_genes, y=Pearson.est, color=Signature, shape=Weight)) +
  geom_point(size=3) + #geom_smooth(se=FALSE, linetype="dashed", alpha=0.5) + 
  theme_classic() + scale_x_continuous(transform = "log10") +
  scale_color_manual(values=fillVals, breaks=names(sig.paths)) + xlab("# of Features") + ylab("Pearson r")
ggsave(paste0("Ngenes_allSigs_PearsonEst_PearsonPSelection_",Sys.Date(),".pdf"), p1, width=5, height=4) # was height 3

fillVals = RColorBrewer::brewer.pal(length(unique(pred$Signature)), "Set2")
p2 <- ggplot(dplyr::distinct(pred[pred$SelectionMetric=="Pearson.p" & !pred$allGenesSameDir & pred$Weighted,]), aes(x=N_genes, y=-log10(Pearson.p), color=Signature, shape=Weight)) +
  geom_point(size=3) + #geom_smooth(se=FALSE, linetype="dashed", alpha=0.5) + 
  theme_classic() + scale_x_continuous(transform = "log10") +
  scale_color_manual(values=fillVals, breaks=names(sig.paths)) + xlab("# of Features") + ylab("-Log(Pearson p)")
ggsave(paste0("Ngenes_allSigs_PearsonP_PearsonPSelection_",Sys.Date(),".pdf"), p2, width=5, height=4) # was height 3

library(patchwork)
p2/p1+plot_layout(guides="collect")
ggsave(paste0("Ngenes_allSigs_PearsonEstORp_PearsonPSelection_",Sys.Date(),".pdf"), width=5, height=4) # was height 3

pred$label <- pred$Genes
pred[pred$Genes=="LRRC25, HMOX1, LRP1, SLC15A3, LILRB2, LILRA6, CHST15, RBM47, SGSH, SLC7A7, TNFRSF1B, LILRB1, CD1D, FGR, IQSEC1, CLEC7A",]$label <- 
  "LRRC25, HMOX1, LRP1, SLC15A3,\nLILRB2, LILRA6, CHST15, RBM47,\nSGSH, SLC7A7, TNFRSF1B, LILRB1,\nCD1D, FGR, IQSEC1, CLEC7A"
pred[pred$Genes=="SLC7A7, LILRA6, LRP1, SGSH, LILRB1, CLEC7A, RBM47, FGR, LRRC25, SLC15A3, LILRB2, IQSEC1, CHST15, TNFRSF1B, HMOX1, CD1D",]$label <- 
  "SLC7A7, LILRA6, LRP1, SGSH,/nLILRB1, CLEC7A, RBM47, FGR,/nLRRC25, SLC15A3, LILRB2, IQSEC1,/nCHST15, TNFRSF1B, HMOX1, CD1D"
ggplot(pred, aes(x=N_genes, y=Pearson.est, color=Signature, shape=Weight, alpha=N)) +
  geom_point() + geom_smooth(se=FALSE, linetype="dashed", show_guides=FALSE) + theme_classic() + scale_x_continuous(transform = "log10") +
  scale_color_manual(values=fillVals, breaks=names(sig.paths)) + 
  geom_hline(yintercept =0, linetype="dashed", color="gray")+labs(x="# of Genes", y="Pearson Correlation Estimate") + 
  ggrepel::geom_label_repel(data=rbind(pred[pred$Signature=="Sorted" & pred$Pearson.q<=0.05,] %>% slice_max(Pearson.est, n=1, with_ties=FALSE),
                                       pred[pred$Signature=="Lasry" & pred$Pearson.q<=0.05,] %>% slice_max(Pearson.est, n=1, with_ties=FALSE)),
                            aes(label=label), box.padding = 1, show_guides=FALSE) + #theme(legend.position="top") + 
  guides(#color=guide_legend(position="right"), 
    alpha=guide_legend(position="bottom"))
ggsave(paste0("Ngenes_allSigs_PearsonEst_wLabelandNalpha_",Sys.Date(),".pdf"), width=5, height=4) # was height 3

#pred <- read.csv("venSensPredictions_2025-07-14.csv")
pred$label <- pred$Genes
pred[pred$Genes=="SLC7A7, LILRA6, LRP1, SGSH, LILRB1, CLEC7A, RBM47, FGR, LRRC25, SLC15A3, LILRB2, IQSEC1, CHST15, TNFRSF1B, HMOX1, CD1D",]$label <- 
  "SLC7A7, LILRA6, LRP1, SGSH,/nLILRB1, CLEC7A, RBM47, FGR,/nLRRC25, SLC15A3, LILRB2, IQSEC1,/nCHST15, TNFRSF1B, HMOX1, CD1D"
fillVals = RColorBrewer::brewer.pal(length(unique(pred$Signature)), "Set2")
topGenesPerSig <- plyr::ddply(pred, .(Signature), summarize,
                              topGeneByQ = Genes[which.min(Pearson.q)],
                              topGeneByEst = Genes[which.max(Pearson.est)])
topIndGenesPerSig <- plyr::ddply(pred[pred$N_genes==1,], .(Signature), summarize,
                              topGeneByQ = Genes[which.min(Pearson.q)],
                              topGeneByEst = Genes[which.max(Pearson.est)]) # same whether by Q or Est
topIndPred <- pred[pred$N_genes == 1 & pred$Genes %in% topIndGenesPerSig$topGeneByQ,]
predMin <- dplyr::distinct(pred[pred$N_genes > 1 | (pred$Genes == topIndGenesPerSig[topIndGenesPerSig$Signature=="Sorted",]$topGeneByEst & pred$Signature == "Sorted") |
                  (pred$Genes == topIndGenesPerSig[topIndGenesPerSig$Signature=="Lasry",]$topGeneByEst & pred$Signature == "Lasry") |
                  (pred$Genes == topIndGenesPerSig[topIndGenesPerSig$Signature=="Triana",]$topGeneByEst & pred$Signature == "Triana") |
                  (pred$Genes == topIndGenesPerSig[topIndGenesPerSig$Signature=="van Galen",]$topGeneByEst & pred$Signature == "van Galen"),])
ggplot(predMin, aes(x=N_genes, y=Pearson.est, color=Signature, shape=Weight, alpha=N)) +
  geom_point() + geom_smooth(se=FALSE, aes(linetype=Weight), show_guides=TRUE) + theme_classic() + scale_x_continuous(transform = "log10") +
  scale_color_manual(values=fillVals, breaks=names(sig.paths)) + scale_linetype_manual(values=c("dashed","dotted"))+
  geom_hline(yintercept =0, linetype="dashed", color="gray")+labs(x="# of Genes", y="Pearson Correlation Estimate") + 
  ggrepel::geom_label_repel(data=rbind(pred[pred$Signature=="Sorted" & pred$Pearson.q<=0.05,] %>% slice_max(Pearson.est, n=1, with_ties=FALSE),
                                       pred[pred$Signature=="Lasry" & pred$Pearson.q<=0.05,] %>% slice_max(Pearson.est, n=1, with_ties=FALSE)),
                            aes(label=label), box.padding = 1, show_guides=FALSE) + #theme(legend.position="top") + 
  guides(#color=guide_legend(position="right"), 
    alpha=guide_legend(position="bottom"))
ggsave(paste0("Ngenes_allSigs_PearsonEst_wLabelandNalpha_",Sys.Date(),".pdf"), width=5, height=4) # was height 3

ggplot(predMin, aes(x=N_genes, y=Pearson.est, color=Signature, shape=Weight)) +
  geom_point() + geom_smooth(se=FALSE, aes(linetype=Weight), show_guides=TRUE) + theme_classic() + scale_x_continuous(transform = "log10") +
  scale_color_manual(values=fillVals, breaks=names(sig.paths)) + scale_linetype_manual(values=c("dashed","dotted"))+
  geom_hline(yintercept =0, linetype="dashed", color="gray")+labs(x="# of Genes", y="Pearson Correlation Estimate") + 
  ggrepel::geom_label_repel(data=rbind(pred[pred$Signature=="Sorted" & pred$Pearson.q<=0.05,] %>% slice_max(Pearson.est, n=1, with_ties=FALSE),
                                       pred[pred$Signature=="Lasry" & pred$Pearson.q<=0.05,] %>% slice_max(Pearson.est, n=1, with_ties=FALSE)),
                            aes(label=label), box.padding = 1, show_guides=FALSE)
ggsave(paste0("Ngenes_allSigs_PearsonEst_wLabel_",Sys.Date(),".pdf"), width=5, height=4) # was height 3
write.csv(predMin, paste0("Ngenes_allSigs_PearsonEst_wLabel_",Sys.Date(),".csv"), row.names=FALSE)

predMinLong <- reshape2::melt(predMin[!predMin$allGenesSameDir,], 
id.vars=c("Drug","Genes", "Gene_directions", "allGenesSameDir", "SelectionMetric", 
"Weighted", "Weight", "Signature", "DataType", "ProteinCoding", "label","N_genes"))
p1 <- ggplot(predMinLong[predMinLong$variable == "Pearson.est" & predMinLong$Weighted,], aes(x=N_genes, y=value, color=Signature, shape=Weight)) +
  geom_point() + geom_smooth(se=FALSE, aes(linetype=Weight), show_guides=TRUE) + theme_classic() + scale_x_continuous(transform = "log10") +
  scale_color_manual(values=fillVals, breaks=names(sig.paths)) + scale_linetype_manual(values=c("dashed","dotted"))+
  #geom_hline(yintercept =0, linetype="dashed", color="gray")+
  labs(x="# of Genes", y = "Pearson r") #+ facet_grid(variable~.) 
p1
#+
  # ggrepel::geom_label_repel(data=rbind(pred[pred$Signature=="Sorted" & pred$Pearson.q<=0.05,] %>% slice_max(Pearson.est, n=1, with_ties=FALSE),
  #                                      pred[pred$Signature=="Lasry" & pred$Pearson.q<=0.05,] %>% slice_max(Pearson.est, n=1, with_ties=FALSE)),
  #                           aes(label=label), box.padding = 1, show_guides=FALSE)

p2 <- ggplot(predMinLong[predMinLong$variable == "Pearson.p" & predMinLong$Weighted,], aes(x=N_genes, y=-log10(value), color=Signature, shape=Weight)) +
  geom_point() + geom_smooth(se=FALSE, aes(linetype=Weight), show_guides=TRUE) + theme_classic() + scale_x_continuous(transform = "log10") +
  scale_color_manual(values=fillVals, breaks=names(sig.paths)) + scale_linetype_manual(values=c("dashed","dotted"))+
  #geom_hline(yintercept =0, linetype="dashed", color="gray")+
  labs(x="# of Genes", y="-Log(P-value)")
p2
# ggrepel::geom_label_repel(data=rbind(pred[pred$Signature=="Sorted" & pred$Pearson.q<=0.05,] %>% slice_max(Pearson.est, n=1, with_ties=FALSE),
#                                      pred[pred$Signature=="Lasry" & pred$Pearson.q<=0.05,] %>% slice_max(Pearson.est, n=1, with_ties=FALSE)),
#                           aes(label=label), box.padding = 1, show_guides=FALSE)
library(patchwork)
p2/p1+plot_layout(guides="collect")
ggsave(paste0("Ngenes_allSigs_PearsonEstORp_wLabel_",Sys.Date(),".pdf"), width=5, height=4) # was height 3
write.csv(predMin, paste0("Ngenes_allSigs_PearsonEst_wLabel_",Sys.Date(),".csv"), row.names=FALSE)

# topPred <- rbind(pred[pred$N_genes==1 & pred$Pearson.q<=0.05,] %>% slice_min(Pearson.p, n=1, with_ties=FALSE),
#                  pred[pred$N_genes==1 & pred$Pearson.q<=0.05,] %>% slice_max(Pearson.est, n=1, with_ties=FALSE),
#                  pred[pred$DataType=="global" & pred$Pearson.q<=0.05,] %>% slice_max(Pearson.est, n=1, with_ties=FALSE),
#                  pred[pred$DataType=="rna" & pred$Pearson.q<=0.05,] %>% slice_max(Pearson.est, n=1, with_ties=FALSE))
topSinglePred <- dplyr::distinct(rbind(singlePred2[singlePred2$DataType=="global" & singlePred2$ProteinCoding,] %>% slice_min(Pearson.p, n=1, with_ties=FALSE),
                                 singlePred2[singlePred2$DataType=="global" & singlePred2$Signature!="None" & singlePred2$ProteinCoding,] %>% slice_min(Pearson.p, n=1, with_ties=FALSE),
                                 singlePred2[singlePred2$DataType=="rna" & singlePred2$ProteinCoding,] %>% slice_min(Pearson.p, n=1, with_ties=FALSE),
                                 singlePred2[singlePred2$DataType=="rna" & singlePred2$Signature!="None" & singlePred2$ProteinCoding,] %>% slice_min(Pearson.p, n=1, with_ties=FALSE)))
write.csv(topSinglePred, "topVenSensPredictions_1gene_2025-07-15.csv", row.names=FALSE)

topPred <- dplyr::distinct(rbind(pred[pred$DataType=="global" & pred$N_genes == 1,] %>% slice_min(Pearson.p, n=1, with_ties=FALSE),
                                       pred[pred$DataType=="global" & pred$N_genes > 1,] %>% slice_min(Pearson.p, n=1, with_ties=FALSE),
                                       pred[pred$DataType=="rna" & pred$N_genes==1,] %>% slice_min(Pearson.p, n=1, with_ties=FALSE),
                                       pred[pred$DataType=="rna" & pred$N_genes>1,] %>% slice_min(Pearson.p, n=1, with_ties=FALSE)))
write.csv(topPred, "topVenSensPredictions_2025-07-15.csv", row.names=FALSE)


topMixPred <- dplyr::distinct(rbind(pred[pred$DataType=="global" & pred$N_genes > 1 & !pred$allGenesSameDir,] %>% slice_min(Pearson.p, n=1, with_ties=FALSE),
                                 pred[pred$DataType=="rna" & pred$N_genes>1 & !pred$allGenesSameDir,] %>% slice_min(Pearson.p, n=1, with_ties=FALSE),
                                 pred[pred$DataType=="global" & pred$N_genes > 1 & !pred$allGenesSameDir,] %>% slice_max(Pearson.est, n=1, with_ties=FALSE),
                                 pred[pred$DataType=="rna" & pred$N_genes>1 & !pred$allGenesSameDir,] %>% slice_max(Pearson.est, n=1, with_ties=FALSE)))
write.csv(topMixPred, "topMixVenSensPredictions_2025-07-15.csv", row.names=FALSE)

# what are the signatures associated with these top mixes?
sigs2 <- list()
for (i in 1:nrow(topMixPred)) {
  temp.sig <- sigs[[topMixPred$Signature[i]]]
  indVenPred <- predAll[predAll$Genes %in% temp.sig$Gene & predAll$DataType == topMixPred$DataType[i],]
  selGenes <- strsplit(topMixPred$Genes[i], ", ")[[1]]
  topGenes <- indVenPred[indVenPred$Pearson.q <= 0.05 & indVenPred$Genes %in% selGenes,c("Genes",topMixPred$Weight[i])]
  temp.name <- paste0(topMixPred$Signature[i], ": ", topMixPred$N_genes[i], 
                      ifelse(topMixPred$Signature[i]=="Sorted", " proteins", " genes"))
  colnames(topGenes)[1] <- "Gene"
  sigs2[[temp.name]] <- topGenes
}
saveRDS(sigs2, "topMixVenSensSigs_2025-07-15.rds")
