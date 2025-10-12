# PTRC2 sorted AML proteomics poster: figure 4
# Author: Belinda B. Garana, belinda.garana@pnnl.gov
remove(list=ls())
library(plyr);library(dplyr);library(ggplot2);library(synapser)
library(patchwork);library(msigdbr)
setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/PTRC2/sortedAMLproteomics/SACB_2025/poster/Figure_4")

synapser::synLogin()

### redo just using significantly corr features
source("~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/panSEA_helper_20240913.R")

#### compile target ranks across all analyses ####
# get DEGs
corr.result <- read.csv(synapser::synGet("syn64543462")$path)
corr.result$Rank <- rank(corr.result$Log2FC)
bestRank <- max(corr.result$Rank)
bestRankCorr <- max(corr.result$Rank)
corr.result$BestRank <- FALSE
corr.result[corr.result$Rank == bestRank,]$BestRank <- TRUE
corr.result$percentile <- percent_rank(corr.result$Log2FC)

# get enriched pathways and their members
all.gsea <- read.csv(synapser::synGet("syn64543470")$path)
# GSEA
hallmark <- msigdbr::msigdbr(species="Homo sapiens",category="H")
hallmark <- dplyr::distinct(hallmark[hallmark$gs_name %in% all.gsea$Feature_set,c("gs_name","gene_symbol")]) # 2358 rows
colnames(hallmark) <- c("Feature_set","Gene")
hallmark <- plyr::ddply(hallmark, .(Feature_set), summarize,
                        Gene = paste0(unique(Gene), collapse = " "))
all.gsea2 <- merge(hallmark, all.gsea, by="Feature_set", all.y = TRUE)
all.gsea2$Rank <- rank(-all.gsea2$NES)
bestRank <- max(all.gsea2$Rank)
bestRankGSEA <- max(all.gsea2$Rank)
all.gsea2$BestRank <- FALSE
all.gsea2[all.gsea2$Rank == bestRank,]$BestRank <- TRUE
all.gsea2$percentile <- percent_rank(all.gsea2$NES)

# get correlated drugs and their targets
drug.corr <- read.csv(synapser::synGet("syn64606618")$path)
drug.corr$Significant <- FALSE
drug.corr[drug.corr$Pearson.q <= 0.05,]$Significant <- TRUE
drug.corr <- drug.corr[drug.corr$Significant,] # 264
beatAML.drug.info <- read.csv("~/OneDrive - PNNL/Documents/PTRC2/BeatAML_single_drug_moa_2025-01-20.csv",
                              stringsAsFactors=FALSE, fileEncoding="latin1")
colnames(beatAML.drug.info)[4] <- "broad_id"
drug.info <- read.csv("https://raw.githubusercontent.com/BelindaBGarana/DMEA/refs/heads/shiny-app/Inputs/PRISM_secondary-screen-replicate-treatment-info.csv")
drug.info <- merge(beatAML.drug.info, drug.info, by="broad_id")
drug.info <- dplyr::distinct(drug.info[,c("Drug","moa.x","moa.y","target")]) # 142
colnames(drug.info)[2] <- "moa"
drug.corr.wInfo <- merge(drug.info, drug.corr, by="Drug") # 45
drug.corr.wInfo <- drug.corr.wInfo[!is.na(drug.corr.wInfo$target),] # 207
drug.targets <- unique(unlist(strsplit(drug.corr.wInfo$target, ", "))) # 291
drug.corr.wInfo$Rank <- rank(-drug.corr.wInfo$Pearson.est)
bestRank <- max(drug.corr.wInfo$Rank)
bestRankDrug <- max(drug.corr.wInfo$Rank)
drug.corr.wInfo$BestRank <- FALSE
drug.corr.wInfo[drug.corr.wInfo$Rank == bestRank,]$BestRank <- TRUE
drug.corr.wInfo$percentile <- percent_rank(-drug.corr.wInfo$Pearson.est)

# # get enriched drug MOAs and their members 
moa.results <- read.csv(synapser::synGet("syn64606616")$path)
moa.results$Significant <- FALSE
moa.results[moa.results$p_value <= 0.05 & moa.results$FDR_q_value <= 0.25,]$Significant <- TRUE
moa.results <- moa.results[moa.results$Significant,] # 19
colnames(moa.results)[3] <- "moa"
moa.info <- dplyr::distinct(drug.info[drug.info$moa %in% moa.results$moa,c("moa","target")]) # 60 MOAs
moa.info <- plyr::ddply(moa.info, .(moa), summarize,
                        Gene = paste0(unique(target), collapse = " "))
moa.results2 <- merge(moa.results, moa.info, by="moa", all.x = TRUE)
moa.results2$Rank <- rank(-moa.results2$NES)
bestRank <- max(moa.results2$Rank)
bestRankMOA <- max(moa.results2$Rank)
moa.results2$BestRank <- FALSE
moa.results2[moa.results2$Rank == bestRank,]$BestRank <- TRUE
moa.results2$percentile <- percent_rank(-moa.results2$NES)

# get selected network nodes and their centrality
net.centr <- read.csv(file.path("~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/analysis/",
                                paste0("PCSF_Protein_", "2025-01-15"),
                                "centrality_optimalRun.csv"))
net.centr$Rank <- rank(net.centr$eigen_centrality)
bestRank <- max(net.centr$Rank)
bestRankNet <- max(net.centr$Rank)
net.centr$BestRank <- FALSE
net.centr[net.centr$Rank == bestRank,]$BestRank <- TRUE
net.centr$percentile <- percent_rank(net.centr$eigen_centrality)

all.features <- unique(c(corr.result$Gene, 
                         unlist(strsplit(all.gsea2$Gene, ", ")), 
                         unlist(strsplit(moa.results2$Drug_set, ", ")), 
                         drug.targets, net.centr$name)) # 7460
other.features <- unique(c(corr.result$Gene, 
                           unlist(strsplit(all.gsea2$Gene, ", ")), 
                           net.centr$name)) # 7427

Target <- drug.targets[drug.targets %in% other.features] # 58; identify drug targets also implicated in at least 1 other analysis
excluded.targets <- drug.targets[!(drug.targets %in% other.features)] # 32; for example, ADRA1A is first and uniprot protein name is ADA1A
# so I think these are gene symbols
toxic <- c("CD14_Pos","CD34_Pos","none")
for (j in toxic) {
  if (j == "CD14_Pos") {
    toxCorr <- corr.result[corr.result$Log2FC > 0,]
    toxGSEA <- all.gsea2[all.gsea2$NES > 0,]
    toxDrug <- drug.corr.wInfo[drug.corr.wInfo$Pearson.est < 0,]
    toxMOA <- moa.results2[moa.results2$NES < 0,]
    toxNet <- net.centr[grepl("positive",net.centr$runID),]
  } else if (j == "CD34_Pos") {
    toxCorr <- corr.result[corr.result$Log2FC < 0,]
    toxGSEA <- all.gsea2[all.gsea2$NES < 0,]
    toxDrug <- drug.corr.wInfo[drug.corr.wInfo$Pearson.est > 0,]
    toxMOA <- moa.results2[moa.results2$NES > 0,]
    toxNet <- net.centr[grepl("negative",net.centr$runID),]
  } else {
    toxCorr <- corr.result
    toxGSEA <- all.gsea2
    toxDrug <- drug.corr.wInfo
    toxMOA <- moa.results2
    toxNet <- net.centr
  }
  percScoresInfo <- data.frame(Target)
  percScoresInfo[,c("Correlation", "GSEA", "Network", 
            "Drug", "DMEA", "MeanPerc", "N","N_best")] <- NA
  percScoresInfo$Best <- ""
  percScores <- percScoresInfo
  nScores <- percScoresInfo
  for (i in 1:length(Target)) {
    nRankVals <- 0
    nBestRanks <- 0
    
    # correlations
    if (Target[i] %in% unique(toxCorr$Gene)) {
      tempCorr <- na.omit(toxCorr[toxCorr$Gene == Target[i],])
      tempTypes <- paste0(unique(tempCorr$type), collapse = " ")
      tempScores <- round(mean(tempCorr$Log2FC), digits=2)
      tempPerc <- round(mean(tempCorr$percentile), digits=2)
      
      # keep score
      nRankVals <- nRankVals + nrow(tempCorr)
      if (any(tempCorr$Rank == bestRankCorr)) {
        nBestRanks <- nBestRanks + length(which(tempCorr$Rank == bestRankCorr))
        percScoresInfo$Best[i] <- paste(percScoresInfo$Best[i], "Correlation")
      }
      
      # update data frame
      percScoresInfo$Correlation[i] <- paste0(tempTypes, " | ", tempScores, " | ", tempPerc)
      percScores$Correlation[i] <- mean(tempCorr$percentile)
      nScores$Correlation[i] <- nrow(tempCorr)
    }
    
    # GSEA
    if (Target[i] %in% unique(unlist(strsplit(toxGSEA$Gene, " ")))) {
      tempGSEA <- na.omit(toxGSEA[grepl(Target[i], toxGSEA$Gene),])
      tempTypes <- paste0(unique(tempGSEA$type), collapse = " ")
      tempSets <- paste0(tempGSEA$Feature_set, collapse = " ")
      tempScores <- round(mean(tempGSEA$NES), digits=2)
      tempPerc <- round(mean(tempGSEA$percentile), digits=2)
      
      # keep score
      nRankVals <- nRankVals + nrow(tempGSEA)
      if (any(tempGSEA$Rank == bestRankGSEA)) {
        nBestRanks <- nBestRanks + length(which(tempGSEA$Rank == bestRankGSEA))
        percScoresInfo$Best[i] <- paste(percScoresInfo$Best[i], "GSEA")
      }
      
      # update data frame
      percScoresInfo$GSEA[i] <- paste0(tempTypes, " | ", tempSets, " | ", tempScores, " | ", tempPerc)
      percScores$GSEA[i] <- mean(tempGSEA$percentile)
      nScores$GSEA[i] <- nrow(tempGSEA)
    }
    
    # Network
    if (Target[i] %in% toxNet$name) {
      tempGSEA <- na.omit(toxNet[toxNet$name == Target[i],])
      tempScores <- round(mean(tempGSEA$eigen_centrality),2)
      tempPerc <- round(mean(tempGSEA$percentile),2)
      
      # keep score
      nRankVals <- nRankVals + nrow(tempGSEA)
      if (any(tempGSEA$Rank == bestRankNet)) {
        nBestRanks <- nBestRanks + length(which(tempGSEA$Rank == bestRankNet))
        percScoresInfo$Best[i] <- paste(percScoresInfo$Best[i], "Network")
      }
      
      # update data frame
      percScoresInfo$Network[i] <- paste0(tempScores, " | ", tempPerc)
      percScores$Network[i] <- mean(tempGSEA$percentile)
      nScores$Network[i] <- nrow(tempGSEA)
    }
    
    # Drug correlation
    if (Target[i] %in% toxDrug$target) {
      tempGSEA <- na.omit(toxDrug[grepl(Target[i], toxDrug$target),])
      tempDrugs <- paste0(tempGSEA$name, collapse = " ")
      tempScores <- round(mean(tempGSEA$Pearson.est),2)
      tempPerc <- round(mean(tempGSEA$percentile),2)
      
      # keep score
      nRankVals <- nRankVals + nrow(tempGSEA)
      if (any(tempGSEA$Rank == bestRankDrug)) {
        nBestRanks <- nBestRanks + length(which(tempGSEA$Rank == bestRankDrug))
        percScoresInfo$Best[i] <- paste(percScoresInfo$Best[i], "Drug")
      }
      
      # update data frame
      percScoresInfo$Drug[i] <- paste0(tempDrugs, " | ", tempScores, " | ", tempPerc)
      percScores$Drug[i] <- mean(tempGSEA$percentile)
      nScores$Drug[i] <- nrow(tempGSEA)
    }
    
    # DMEA
    if (Target[i] %in% unique(unlist(strsplit(toxMOA$Drug_set, " ")))) {
      tempGSEA <- na.omit(toxMOA[grepl(Target[i], toxMOA$Drug_set),])
      tempTypes <- paste0(unique(tempGSEA$type), collapse = " ")
      tempSets <- paste0(tempGSEA$moa, collapse = " ")
      tempScores <- round(mean(tempGSEA$NES),2)
      tempPerc <- round(mean(tempGSEA$percentile),2)
      
      # keep score
      nRankVals <- nRankVals + nrow(tempGSEA)
      if (any(tempGSEA$Rank == bestRankMOA)) {
        nBestRanks <- nBestRanks + length(which(tempGSEA$Rank == bestRankMOA))
        percScoresInfo$Best[i] <- paste(percScoresInfo$Best[i], "DMEA")
      }
      
      # update data frame
      percScoresInfo$DMEA[i] <- paste0(tempTypes, " | ", tempSets, " | ", tempScores, " | ", tempPerc)
      percScores$DMEA[i] <- mean(tempGSEA$percentile)
      nScores$DMEA[i] <- nrow(tempGSEA)
    }
    
    # update data frame
    percScores$N[i] <- nRankVals
    percScores$N_best[i] <- nBestRanks
    percScores$MeanPerc[i] <- mean(as.numeric(percScores[i,2:6]), na.rm = TRUE)
  }
  percScoresInfo$N_best <- percScores$N_best
  percScoresInfo$N <- percScores$N
  percScoresInfo$MeanPerc <- percScores$MeanPerc
  percScoresInfo$Best <- substring(percScoresInfo$Best, 2)
  nScores$N_best <- percScoresInfo$N_best
  nScores$N <- percScoresInfo$N
  nScores$MeanPerc <- percScoresInfo$MeanPerc
  nScores$Best <- percScoresInfo$Best
  percScores$Best <- percScoresInfo$Best
  percScores$MeanPercTimesN <- percScores$MeanPerc*percScores$N
  percScoresInfo$MeanPercTimesN <- percScores$MeanPercTimesN
  
  percScores$N_analyses <- rowSums(!is.na(percScores[,2:6]))
  percScores$MeanPercTimesN_analyses <- percScores$N_analyses * percScores$MeanPerc
  percScoresInfo$MeanPercTimesN_analyses <- percScores$MeanPercTimesN_analyses
  write.csv(nScores, paste0(j,"_N_Percentile_scores_",Sys.Date(),".csv"), row.names = FALSE)
  #nScores <- read.csv(paste0(j,"_N_Percentile_scores_",Sys.Date(),".csv"))
  write.csv(percScores, paste0(j,"_Percentile_scores_",Sys.Date(),".csv"), row.names = FALSE)
  write.csv(percScoresInfo, paste0(j, "_Percentile_scores_info_",Sys.Date(),".csv"), row.names = FALSE)
  
  filtered.percScores <- percScores[percScores$N_analyses >= 3,]
  if (nrow(filtered.percScores) > 0) {
    filtered.percScores <- filtered.percScores[order(-filtered.percScores$MeanPercTimesN_analyses),]
    targetOrder <- filtered.percScores$Target
    write.csv(filtered.percScores, paste0(j,"_Percentile_scores_min3analyses_",Sys.Date(),".csv"), row.names = FALSE)

    filtered.percScores$N <- NULL
    filtered.percScores$N_best <- NULL
    filtered.percScores$N_analyses <- NULL
    filtered.percScores$MeanPerc <- NULL
    filtered.percScores$MeanPercTimesN <- NULL
    filtered.percScores$MeanPercTimesN_analyses <- NULL
    nScores2 <- reshape2::melt(nScores[nScores$Target %in% filtered.percScores$Target,
                                       colnames(filtered.percScores)], 
                               id = c("Target","Best"), variable.name = "Analysis", value.name="N")
    filtered.percScores <- reshape2::melt(filtered.percScores, id = c("Target", "Best"), 
                                          variable.name = "Analysis", value.name = "Mean Percentile")
    
    plot.df <- merge(filtered.percScores, nScores2, by=c("Target","Best","Analysis"))
    plot.df$`Mean Percentile` <- as.numeric(plot.df$`Mean Percentile`)
    plot.df$N <- as.numeric(plot.df$N)
    plot.df <- na.omit(plot.df)
    plot.df$`Top Score` <- FALSE
    if (any(plot.df$Best == plot.df$Analysis)) {
      plot.df[plot.df$Best == plot.df$Analysis,]$`Top Score` <- TRUE 
    }
    # create dot plot to represent filtered.percScores
    dot.plot <- ggplot2::ggplot(plot.df, aes(y = Target, x = Analysis, color = 100*`Mean Percentile`, size=N)) + 
      ggplot2::geom_point() + scale_color_gradient2(low="grey",high="red", limits=c(0, 100)) +
      theme_classic() + ggplot2::labs(color = "Mean Percentile", size = "N") + 
      ggplot2::scale_y_discrete(limits = targetOrder) +
      geom_point(data = subset(plot.df, `Top Score`), col = "black", stroke = 1.5, shape = 21) +
      theme(axis.text.y = element_text(angle = 45, vjust=1, hjust=1),
            axis.text.x = element_text(angle=90,vjust=0.5),
            axis.title.x=element_text(angle=180)
      ) + theme(#legend.position="top", 
        legend.direction="vertical") + 
      coord_flip()
    dot.plot
    saveRDS(dot.plot, paste0(j,"_Target_min3analyses_dotPlot_",Sys.Date(),".rds"))
    ggsave(paste0(j,"_Target_min3analyses_dotPlot_wider_",Sys.Date(),".pdf"),dot.plot,width=11, height=4)
    
    dot.plot <- ggplot2::ggplot(plot.df, aes(y = Target, x = Analysis, color = 100*`Mean Percentile`, size=N)) + 
      ggplot2::geom_point() + scale_color_gradient2(low="grey",high="red", limits=c(0, 100)) +
      theme_classic() + ggplot2::labs(color = "Mean Percentile", size = "N") + 
      ggplot2::scale_y_discrete(limits = targetOrder) +
      geom_point(data = subset(plot.df, `Top Score`), col = "black", stroke = 1.5, shape = 21) +
      theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1),
            #axis.text.x = element_text(angle=90,vjust=0.5),
            #axis.title.x=element_text(angle=180)
      ) + theme(#legend.position="top", 
        legend.direction="vertical") + 
      coord_flip()
    dot.plot
    saveRDS(dot.plot, paste0(j,"_Target_min3analyses_dotPlot_forHorizontal_",Sys.Date(),".rds"))
    ggsave(paste0(j,"_Target_min3analyses_dotPlot_forHorizontal_",Sys.Date(),".pdf"),dot.plot,width=4, height=2)
    ggsave(paste0(j,"_Target_min3analyses_dotPlot_forHorizontal_wider_",Sys.Date(),".pdf"),dot.plot,width=11, height=4)
  }
}
dot.plot.amp <- readRDS(paste0("CD14_Pos","_Target_min3analyses_dotPlot_forHorizontal_",Sys.Date(),".rds"))
dot.plot.amp
ggsave(paste0("CD14_Pos","_Target_min3analyses_dotPlot_forHorizontal",Sys.Date(),".pdf"),dot.plot,width=5, height=4)
ggsave(paste0("CD14_Pos","_Target_min3analyses_dotPlot_v2_",Sys.Date(),".pdf"),dot.plot,width=10.5, height=2.5)
