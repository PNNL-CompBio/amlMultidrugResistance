# plots for SACB poster re: sorted AML proteomics
library(plyr);library(dplyr);library(synapser);library(ggplot2)
setwd("~/OneDrive - PNNL/Documents/PTRC2/sortedAMLproteomics/SACB_2025/poster/")

n.top <- 5
synapser::synLogin()
#### 1. diffexp bar plot: top 20 ####
# load input
diffexp <- na.omit(read.csv(synapser::synGet("syn69336657")$path)) # was syn64543462 before considering cell type, sort type, patient factors in differential expression
sig.diffexp <- diffexp[diffexp$adj.P.Val <= 0.05,]
top.pos.diffexp <- sig.diffexp %>% slice_max(Log2FC, n = n.top)
top.neg.diffexp <- sig.diffexp %>% slice_min(Log2FC, n = n.top)
top.diffexp <- rbind(top.pos.diffexp, top.neg.diffexp)
top.diffexp$Significance <- "Upregulated"
top.diffexp[top.diffexp$Log2FC < 0,]$Significance <- "Downregulated"
top.diffexp$Significance <- factor(top.diffexp$Significance, levels=c("Upregulated","Downregulated"))

# vertical bar plot of top 5 sig results
nSig <- nrow(sig.diffexp)
nTotal <- nrow(diffexp)
title <- paste0("Differentially expressed genes (",nSig,"/",nTotal, " with adjusted p <= 0.05)")
#setOrder <- top.diffexp[order(top.diffexp$Log2FC*-log10(top.diffexp$adj.P.Val)),]$Gene
setOrder <- top.diffexp[order(top.diffexp$Log2FC),]$Gene

# bold CD14 and CD34
geneFace <- rep("plain", nrow(top.diffexp))
names(geneFace) <- setOrder
geneFace[c("CD14","CD34")] <- "bold"

diffexp.bar <- ggplot(top.diffexp, aes(x=Gene, y=Log2FC, fill=Significance, alpha=0.5)) + geom_col()+
  ggplot2::scale_x_discrete(limits = setOrder) +
  scale_fill_manual(values=c("red","blue"), breaks = c("Upregulated", "Downregulated"))+
  theme(axis.text.x=element_text(face=geneFace)) +
  theme_classic(base_size = 12) + ggtitle(title) + coord_flip()
diffexp.bar
ggsave(paste0("diffexp_top_",n.top,"_sig_absLog2FC_barPlot_",Sys.Date(),".pdf"), diffexp.bar, width=7, height=7)

title <- paste0("Differentially expressed genes\n(",nSig,"/",nTotal, " with adjusted p <= 0.05)")
top.diffexp$`-Log(FDR)` <- -log(top.diffexp$adj.P.Val, base=10)
absMaxLog2FC <- max(abs(top.diffexp$Log2FC))
maxLogFDR <- ceiling(max(top.diffexp$`-Log(FDR`))
top.diffexp$Significant <- TRUE
diffexp.dot <- ggplot(top.diffexp, aes(x=Gene, y="CD14+ vs. CD34+", color=Log2FC, size=`-Log(FDR)`)) + geom_point()+
  ggplot2::scale_x_discrete(limits = setOrder) +
  theme(axis.text.x=element_text(face=geneFace)) +
  scale_size(limits=c(1,maxLogFDR), range = c(0.5,7)) +
  scale_color_gradient2(low="blue",high="red", mid="grey", limits=c(-absMaxLog2FC, absMaxLog2FC)) +
  geom_point(data = subset(top.diffexp, Significant), col = "black", stroke = 1.5, shape = 21) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text = element_text(size=16)) + 
  theme_classic(base_size = 12) + ggtitle(title) + coord_flip()
diffexp.dot
ggsave(paste0("diffexp_top_",n.top,"_sig_absLog2FC_dotPlot_",Sys.Date(),".pdf"), diffexp.dot, width=4, height=4)

n.top <- 10
#### 2. GSEA bar plot: Hallmark ####
# load input
gsea <- read.csv(synapser::synGet("syn69336665")$path) # was syn64543470 before considering cell type, sort type, patient factors in differential expression
gsea$Significance <- "FDR > 0.25"
gsea[gsea$FDR_q_value <= 0.25 & gsea$p_value <= 0.05,]$Significance <- "FDR <= 0.25"
gsea$Significance <- factor(gsea$Significance, levels = c("FDR <= 0.25", "FDR > 0.25"))
gsea$`Gene Set` <- sub("HALLMARK_","",gsea$Feature_set)
top.gsea <- gsea[gsea$Significance == "FDR <= 0.25",] %>% slice_max(abs(NES), n = n.top)

# vertical bar plot of top sig results
nSig <- nrow(gsea[gsea$Significance == "FDR <= 0.25",]) # 24
nTotal <- nrow(gsea)
title <- paste0("Hallmark Gene Sets (",nSig,"/",nTotal, " Enriched)")
setOrder <- top.gsea[order(top.gsea$NES),]$`Gene Set`
gsea.bar <- ggplot(top.gsea, aes(x=`Gene Set`, y=NES, fill=Significance)) + geom_col()+
  ggplot2::scale_x_discrete(limits = setOrder) +
  scale_fill_manual(values=c("red","grey"), breaks = c("FDR <= 0.25", "FDR > 0.25"))+
  theme_classic(base_size = 12) + ggtitle(title) + coord_flip()
gsea.bar
ggsave(paste0("GSEA_Hallmark_top_",n.top,"_sets_barPlot_",Sys.Date(),".pdf"), gsea.bar, width=7, height=7)

title <- paste0("Hallmark Gene Sets\n(",nSig,"/",nTotal, " Enriched)")
top.gsea$`-Log(FDR)` <- -log(top.gsea$FDR_q_value, base=10)
top.gsea[top.gsea$FDR_q_value==0,]$`-Log(FDR)` <- -log(1E-4, base=10)
absMaxNES <- max(abs(top.gsea$NES))
maxLogFDR <- ceiling(max(top.gsea$`-Log(FDR`))
top.gsea$Significant <- TRUE
gsea.dot <- ggplot(top.gsea, aes(x=`Gene Set`, y="CD14+ vs. CD34+", color=NES, size=`-Log(FDR)`)) + geom_point()+
  ggplot2::scale_x_discrete(limits = setOrder) +
  scale_size(limits=c(2,maxLogFDR), range = c(1,5)) +
  scale_color_gradient2(low="blue",high="red", mid="grey", limits=c(-absMaxNES, absMaxNES)) +
  geom_point(data = subset(top.gsea, Significant), col = "black", stroke = 1.5, shape = 21) + 
  theme_classic(base_size = 12) + ggtitle(title) + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), 
        axis.text = element_text(size=16)) + coord_flip()
gsea.dot
ggsave(paste0("gsea_top_",n.top,"_sig_absNES_",Sys.Date(),".pdf"), gsea.dot, width=5.5, height=4) # height 2 for n.top = 5, 5 for n.top = 20

title <- paste0("Hallmark Gene Sets\n(",nSig,"/",nTotal, " Enriched)")
gsea$`-Log(FDR)` <- -log(gsea$FDR_q_value, base=10)
gsea[gsea$FDR_q_value==0,]$`-Log(FDR)` <- -log(1E-4, base=10)
absMaxNES <- max(abs(gsea$NES))
maxLogFDR <- ceiling(max(gsea$`-Log(FDR`))
gsea$Significant <- FALSE
gsea[gsea$FDR_q_value <= 0.25 & gsea$p_value <= 0.05,]$Significant <- TRUE
sig.gsea <- gsea[gsea$Significant,]
setOrder <- sig.gsea[order(sig.gsea$NES),]$`Gene Set`
gsea.dot <- ggplot(sig.gsea, aes(x=`Gene Set`, y="CD14+ vs. CD34+", color=NES, size=`-Log(FDR)`)) + geom_point()+
  ggplot2::scale_x_discrete(limits = setOrder) +
  scale_size(limits=c(0.75,maxLogFDR), range = c(0.5,5)) +
  scale_color_gradient2(low="blue",high="red", mid="grey", limits=c(-absMaxNES, absMaxNES)) +
  geom_point(data = sig.gsea, col = "black", stroke = 1.5, shape = 21) + 
  theme_classic(base_size = 12) + ggtitle(title) + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), 
        axis.text = element_text(size=16)) + coord_flip()
gsea.dot
ggsave(paste0("gsea_sig_",Sys.Date(),".pdf"), gsea.dot, width=5.5, height=6) # height 2 for n.top = 5, 5 for n.top = 20


setOrder <- gsea[order(gsea$NES),]$`Gene Set`
gsea.dot <- ggplot(gsea, aes(x=`Gene Set`, y="CD14+ vs. CD34+", color=NES, size=`-Log(FDR)`)) + geom_point()+
  ggplot2::scale_x_discrete(limits = setOrder) +
  scale_size(limits=c(0,maxLogFDR), range = c(0.5,5)) +
  scale_color_gradient2(low="blue",high="red", mid="grey", limits=c(-absMaxNES, absMaxNES)) +
  geom_point(data = subset(gsea, Significant), col = "black", stroke = 1.5, shape = 21) + 
  theme_classic(base_size = 12) + ggtitle(title) + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), 
        axis.text = element_text(size=16)) + coord_flip()
gsea.dot
ggsave(paste0("gsea_",Sys.Date(),".pdf"), gsea.dot, width=6.25, height=12) # height 2 for n.top = 5, 5 for n.top = 20


gsea <- read.csv(synapser::synGet("syn69336665")$path) # was syn64543470 before considering cell type, sort type, patient factors in differential expression
top.gsea <- gsea %>% slice_max(NES, n = n.top)
top.gsea$`Gene Set` <- sub("HALLMARK_","",top.gsea$Feature_set)

# vertical bar plot of top 5 sig results
nSig <- nrow(gsea[gsea$FDR_q_value <= 0.25 & gsea$p_value <= 0.05,])
nTotal <- nrow(gsea)
title <- paste0("Hallmark Gene Sets\n(",nSig,"/",nTotal, " Enriched)")
setOrder <- top.gsea[order(top.gsea$NES),]$`Gene Set`
gsea.bar <- ggplot(top.gsea, aes(x=`Gene Set`, y=NES, fill="red", alpha=0.5)) + geom_col()+
  ggplot2::scale_x_discrete(limits = setOrder) +
  theme_classic(base_size = 12) + ggtitle(title) + coord_flip()
gsea.bar
ggsave(paste0("GSEA_Hallmark_top_",n.top,"_upregulated_sets_barPlot_",Sys.Date(),".pdf"), gsea.bar, width=5, height=3)


#### 3. Diffexp dot plot: leading genes ####
# load input
goi <- stringr::str_split(top.gsea[top.gsea$`Gene Set` == "TNFA_SIGNALING_VIA_NFKB",]$Leading_edge, ", ")[[1]] # 29
goi.diffexp <- diffexp[diffexp$Gene %in% goi,]
sig.goi.diffexp <- goi.diffexp[goi.diffexp$adj.P.Val<=0.05,] # 11
nSig <- nrow(sig.goi.diffexp)
nTotal <- nrow(goi.diffexp)
title <- paste0("NFKB signaling genes\n(",nSig,"/",nTotal, " with adjusted p <= 0.05)")
sig.goi.diffexp$`-Log(FDR)` <- -log(sig.goi.diffexp$adj.P.Val, base=10)
absMaxLog2FC <- max(abs(sig.goi.diffexp$Log2FC))
maxLogFDR <- ceiling(max(sig.goi.diffexp$`-Log(FDR`))
sig.goi.diffexp$Significant <- TRUE
setOrder <- sig.goi.diffexp[order(sig.goi.diffexp$Log2FC),]$Gene
diffexp.dot <- ggplot(sig.goi.diffexp, aes(x=Gene, y="CD14+ vs. CD34+", color=Log2FC, size=`-Log(FDR)`)) + geom_point()+
  ggplot2::scale_x_discrete(limits = setOrder) +
  scale_size(limits=c(1,maxLogFDR), range = c(0.5,7)) +
  scale_color_gradient2(low="blue",high="red", mid="grey", limits=c(-absMaxLog2FC, absMaxLog2FC)) +
  geom_point(data = subset(sig.goi.diffexp, Significant), col = "black", stroke = 1.5, shape = 21) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text = element_text(size=16)) + 
  theme_classic(base_size = 12) + ggtitle(title) + coord_flip()
diffexp.dot
ggsave(paste0("NFKB_leadingEdge_sig_diffexp_dotPlot_",Sys.Date(),".pdf"), diffexp.dot, width=3, height=4)

diffexp.dot <- ggplot(sig.goi.diffexp, aes(x=Gene, y="CD14+ vs. CD34+", color=Log2FC, size=`-Log(FDR)`)) + geom_point()+
  ggplot2::scale_x_discrete(limits = rev(setOrder)) +
  scale_size(limits=c(1,maxLogFDR), range = c(0.5,7)) +
  scale_color_gradient2(low="blue",high="red", mid="grey", limits=c(-absMaxLog2FC, absMaxLog2FC)) +
  geom_point(data = subset(sig.goi.diffexp, Significant), col = "black", stroke = 1.5, shape = 21) +
  theme_classic(base_size = 12) + ggtitle(title) +
  theme(axis.title.y=element_blank(), #axis.text = element_text(size=16),
        axis.text.x = element_text(angle = 45, vjust=1, hjust=1))
diffexp.dot
ggsave(paste0("NFKB_leadingEdge_sig_diffexp_dotPlot_horizontal_",Sys.Date(),".pdf"), diffexp.dot, width=5, height=2)
ggsave(paste0("NFKB_leadingEdge_sig_diffexp_dotPlot_horizontal_taller_",Sys.Date(),".pdf"), diffexp.dot, width=5, height=3)
#### 2. DMEA bar plot ####
# load input
moa.results <- read.csv(synapser::synGet("syn69928448")$path) # was syn64606616 before considering cell type, sort type, patient factors in differential expression
top.gsea <- moa.results %>% slice_max(abs(NES), n = n.top)
top.gsea$Significance <- "FDR > 0.25"
top.gsea[top.gsea$FDR_q_value <= 0.25 & top.gsea$p_value <= 0.05,]$Significance <- "FDR <= 0.25"
top.gsea$Significance <- factor(top.gsea$Significance, levels = c("FDR <= 0.25", "FDR > 0.25"))
top.gsea$`Drug Mechanism` <- top.gsea$Drug_set

# vertical bar plot of top 5 results
nSig <- nrow(top.gsea[top.gsea$Significance == "FDR <= 0.25",])
nTotal <- nrow(moa.results)
title <- paste0("Drug Mechanisms (",nSig,"/",nTotal, " Enriched)")
setOrder <- top.gsea[order(top.gsea$NES, decreasing = TRUE),]$`Drug Mechanism`
gsea.bar <- ggplot(top.gsea, aes(x=`Drug Mechanism`, y=NES, fill=Significance)) + geom_col()+
  ggplot2::scale_x_discrete(limits = setOrder) +
  scale_fill_manual(values=c("red","grey"), breaks = c("FDR <= 0.25", "FDR > 0.25"))+
  theme_classic(base_size = 12) + ggtitle(title) + coord_flip()
ggsave(paste0("DMEA_top_",n.top,"_MOAs_barPlot_",Sys.Date(),".pdf"), gsea.bar, width=7, height=7)

# only focus on significant MOAs for resistant CD14+ samples
sig.moa.results <- moa.results[moa.results$FDR_q_value <= 0.25 & moa.results$p_value <= 0.05,]
top.gsea <- sig.moa.results %>% slice_max(NES, n = n.top)
top.gsea$`Drug Mechanism` <- top.gsea$Drug_set

# vertical bar plot of top 5 results
nSig <- nrow(sig.moa.results)
nTotal <- nrow(moa.results)
title <- paste0("Drug Mechanisms\n(",nSig,"/",nTotal, " Enriched)")
setOrder <- top.gsea[order(top.gsea$NES, decreasing = TRUE),]$`Drug Mechanism`
gsea.bar <- ggplot(top.gsea, aes(x=`Drug Mechanism`, y=-NES, fill="red", alpha=0.5)) + geom_col()+
  ggplot2::scale_x_discrete(limits = setOrder) +
  theme_classic(base_size = 12) + ggtitle(title) + coord_flip()
gsea.bar
ggsave(paste0("DMEA_top_",n.top,"_MOAs_forCD14PosSamples_barPlot_",Sys.Date(),".pdf"), gsea.bar, width=4, height=2)

#### 3. Drug correlation bar plot ####
sig.moas <- unique(sig.moa.results$Drug_set)
drug.corr <- read.csv(synapser::synGet("syn69928450")$path) # was syn64606618  before considering cell type, sort type, patient factors in differential expression
drug.info <- read.csv("~/OneDrive - PNNL/Documents/PTRC2/BeatAML_single_drug_moa.csv",
                      stringsAsFactors = FALSE, fileEncoding = "latin1")
drug.info <- drug.info[,c("Drug","moa")]
drug.corr <- drug.corr[drug.corr$Drug %in% drug.info$Drug,]

maxAbsEst <- max(na.omit(abs(drug.corr$Pearson.est)))
#maxAbsEst <- 1
drug.corr <- na.omit(drug.corr)
if (any(drug.corr$Pearson.q == 0)) {
  drug.corr[drug.corr$Pearson.q == 0,]$Pearson.q <- 0.0001
}
drug.corr$minusLogP <- -log(drug.corr[,"Pearson.p"], base = 10)
drug.corr$minusLogFDR <- -log(drug.corr[,"Pearson.q"], base = 10)
drug.corr.wInfo <- merge(drug.info, drug.corr, by="Drug", all.y = TRUE)
drug.corr.wInfo$Mechanism <- "Other"

moa.results$Significant <- FALSE
moa.results[moa.results$p_value <= 0.05 & moa.results$FDR_q_value <= 0.25,]$Significant <- TRUE
sig.moas <- unique(moa.results[moa.results$Significant,]$Drug_set) # 2
mean.moa <- plyr::ddply(moa.results,.(Drug_set),
                        NES = mean(NES, na.rm=TRUE))
moaOrder <- mean.moa[order(mean.moa$NES),]$Drug_set
#drug.corr.wInfo[drug.corr.wInfo$moa %in% sig.moas,]$Mechanism <- drug.corr.wInfo[drug.corr.wInfo$moa %in% sig.moas,]$moa
drug.corr.wInfo$Mechanism <- "Other"
drug.corr.wInfo[drug.corr.wInfo$moa %in% moa.results$Drug_set,]$Mechanism <- drug.corr.wInfo[drug.corr.wInfo$moa %in% moa.results$Drug_set,]$moa

maxAbsEst <- max(na.omit(abs(drug.corr.wInfo$Pearson.est)))
maxLogFDR <- max(na.omit(drug.corr.wInfo$minusLogFDR))
library(patchwork); library(ggplot2)
drug.corr.wInfo$Significant <- FALSE
drug.corr.wInfo[drug.corr.wInfo$Pearson.q <= 0.05,]$Significant <- TRUE
drug.corr.wInfo[drug.corr.wInfo$Drug == "Ralimetinib (LY2228820)",]$moa <- "p38 MAPK inhibitor"
drug.corr.wInfo[drug.corr.wInfo$Drug == "Nilotinib",]$moa <- "Abl kinase inhibitor"
drug.corr.wInfo[drug.corr.wInfo$Drug == "AT-101",]$moa <- "BCL inhibitor"
drug.corr.wInfo$Drug <- sub(" [(].*", "", drug.corr.wInfo$Drug) # shorten drug names for plot
drug.corr.wInfo[drug.corr.wInfo$Drug == "NF-kB Activation Inhibitor",]$Drug <- "NFkB Inhibitor"

n.top.list <- c(5, 10, 15)
n.top.list <- c(10, 15, 20)
for (n.top in n.top.list) {
  temp.dot.df <- drug.corr.wInfo[drug.corr.wInfo$Significant,]
  up.dot.df <- temp.dot.df %>% slice_max(Pearson.est,n=n.top)
  dn.dot.df <- temp.dot.df %>% slice_min(Pearson.est,n=n.top)
  temp.dot.df <- rbind(up.dot.df, dn.dot.df)
  temp.dot.df$Mechanism <- temp.dot.df$moa
  temp.dot.df[grepl(", ", temp.dot.df$Mechanism),]$Mechanism <- "Other"
  MOAsInTop50 <- unique(temp.dot.df$Mechanism)

  nGenes <- length(unique(na.omit(drug.corr.wInfo)$Drug))
  nCorrGenes <- length(unique(na.omit(drug.corr.wInfo[drug.corr.wInfo$Significant,])$Drug))
  
  geneOrder <- temp.dot.df[order(temp.dot.df$Pearson.est),]$Drug
  
  plot.annot <- paste0("Drug Sensitivity", "\n(", nCorrGenes, " / ", nGenes, " Drugs Correlated)")
  
  dot.plot <- ggplot2::ggplot(temp.dot.df,
                              ggplot2::aes(
                                x = Drug, y = Pearson.est, fill = Mechanism
                              )
  ) + scale_size(limits=c(0,maxLogFDR), range = c(0.5,4)) +
    ggplot2::geom_col() +
    ggplot2::scale_x_discrete(limits = geneOrder) +
    theme_classic(base_size = 12) + scale_fill_manual(breaks=c(MOAsInTop50[MOAsInTop50 != "Other"],"Other"), 
                                        values = c(grDevices::colorRampPalette(
                                          RColorBrewer::brewer.pal(12, "Set3"))(length(MOAsInTop50)-1),"grey")) +
    ggplot2::labs(
      #x = i,
      y = "Pearson r",
      fill = "Drug Mechanism"
    ) + theme(axis.title.x = element_blank()) +
    #theme(axis.title.y=element_text(face="bold", size=16), 
    #axis.text.y = element_blank(),
    #axis.ticks.y = element_blank(),
    #axis.title.x=element_blank()) +
    theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)
          #, 
          #axis.text=element_text(size=16)
          ) +
    theme(legend.direction = "horizontal", legend.position = "bottom")
  cat(plot.annot, "\n")
  dot.plot <- dot.plot + ggtitle(plot.annot) + theme(plot.title = element_text(hjust = 0.5, face="bold", size=16))
  dot.plot
  ggplot2::ggsave(paste0("CorrelatedDrugs_barPlot_moaFill_sliceMaxPearson",n.top,"_",Sys.Date(),".pdf"), dot.plot, width=6, height=5)
  ggplot2::ggsave(paste0("CorrelatedDrugs_barPlot_moaFill_sliceMaxPearson",n.top,"_wider_",Sys.Date(),".pdf"), dot.plot, width=11, height=5)
  
  dot.plot <- ggplot2::ggplot(temp.dot.df,
                              ggplot2::aes(
                                x = Drug, y = Pearson.est, fill = Mechanism
                              )
  ) + scale_size(limits=c(0,maxLogFDR), range = c(0.5,4)) +
    ggplot2::geom_col() +
    ggplot2::scale_x_discrete(limits = geneOrder) +
    theme_classic(base_size = 16) + scale_fill_manual(breaks=c(MOAsInTop50[MOAsInTop50 != "Other"],"Other"), 
                                                      values = c(grDevices::colorRampPalette(
                                                        RColorBrewer::brewer.pal(12, "Set3"))(length(MOAsInTop50)-1),"grey")) +
    ggplot2::labs(
      #x = i,
      y = "Pearson r",
      fill = "Drug Mechanism"
    ) + theme(axis.title.x = element_blank()) +
    #theme(axis.title.y=element_text(face="bold", size=16), 
    #axis.text.y = element_blank(),
    #axis.ticks.y = element_blank(),
    #axis.title.x=element_blank()) +
    theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)
          #, 
          #axis.text=element_text(size=16)
    ) +
    theme(legend.direction = "vertical", legend.position = "right")
  cat(plot.annot, "\n")
  dot.plot <- dot.plot + ggtitle(plot.annot) + theme(plot.title = element_text(hjust = 0.5, face="bold", size=16))
  dot.plot
  ggplot2::ggsave(paste0("CorrelatedDrugs_barPlot_moaFill_sliceMaxPearson",n.top,"_legendRight_",Sys.Date(),".pdf"), dot.plot, width=8, height=5)
  ggplot2::ggsave(paste0("CorrelatedDrugs_barPlot_moaFill_sliceMaxPearson",n.top,"_wider_legendRight_",Sys.Date(),".pdf"), dot.plot, width=11, height=5)
  
  # dot.plot <- ggplot2::ggplot(temp.dot.df,
  #                             ggplot2::aes(
  #                               x = Drug, y = Pearson.est, fill = Mechanism
  #                             )
  # ) + scale_size(limits=c(0,maxLogFDR), range = c(0.5,4)) +
  #   ggplot2::geom_col() +
  #   ggplot2::scale_x_discrete(limits = geneOrder) +
  #   theme_classic() + scale_fill_manual(breaks=c(moaOrder3,"Other"), 
  #                                       values = c(grDevices::colorRampPalette(
  #                                         RColorBrewer::brewer.pal(length(MOAsInTop50), "Set3"))(length(MOAsInTop50)))) +
  #   ggplot2::labs(
  #     #x = i,
  #     y = "Pearson r",
  #     fill = "Drug Mechanism"
  #   ) + theme(axis.title.x = element_blank()) +
  #   #theme(axis.title.y=element_text(face="bold", size=16), 
  #   #axis.text.y = element_blank(),
  #   #axis.ticks.y = element_blank(),
  #   #axis.title.x=element_blank()) +
  #   theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
  #   theme(legend.direction = "horizontal", legend.position = "bottom")
  # cat(plot.annot, "\n")
  # dot.plot <- dot.plot + ggtitle(plot.annot) + theme(plot.title = element_text(hjust = 0.5, face="bold", size=16))
  # dot.plot
  # ggplot2::ggsave(paste0("CorrelatedDrugs_barPlot_moaFill_sliceMaxAbsPearson)",n.top,".pdf"), dot.plot, width=10, height=7)
}

