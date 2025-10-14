# figure 3 for NRAS manuscript: cell line ASO (exp 20)
library(plyr);library(dplyr); library(stringr); library(tidyr); library(readxl); library(scales);
library(synapser); library(ggplot2); library(MSnSet.utils); library(patchwork);library(data.table)
library(PCSF)
synapser::synLogin()
base.path <- "~/OneDrive - PNNL/Documents/GitHub/Exp21_NRAS_ASO_treated_patients/proteomics/analysis/"
setwd(base.path)
dir.create("Figures")
setwd("Figures")
dir.create("Figure6")
setwd("Figure6")

#### RAS/MAP/ERK diffexp dot plots ####
# synIDs.mix <- list("Global" = list("NRAS ASO" = "syn64732753",
#                                    "NRAS ASO + Gilt" = "syn64732932"),
#                    "Phospho" = list("NRAS ASO" = "syn64732862",
#                                     "NRAS ASO + Gilt" = "syn64733038"))
synIDs.mix <- list("Global" = list("NRAS ASO" = "syn66726954",
                                   "NRAS ASO + Gilt" = "syn66727021"),
                   "Phospho" = list("NRAS ASO" = "syn66726980",
                                    "NRAS ASO + Gilt" = "syn66726998"))
all.synIDs <- list("Mix" = synIDs.mix)
all.de <- list()
for (t in names(all.synIDs)) {
  synIDs <- all.synIDs[[t]]
  for (omics in names(synIDs)) {
    # load diffexp
    time.de <- all.de[[t]]
    if (omics %in% names(all.de)) {
      de <- time.de[[omics]]
    } else {
      de <- data.frame()
      for (contrast in names(synIDs[[omics]])) {
        temp.de <- read.csv(synapser::synGet(synIDs[[omics]][[contrast]])$path) # uses all timepoints as factor
        temp.de$Contrast <- contrast
        de <- rbind(de, temp.de)
      }
      de$Contrast <- factor(de$Contrast, levels=names(synIDs[[omics]]))
      colnames(de)[1] <- "feature"
      time.de[[omics]] <- de 
    }
    all.de[[t]] <- time.de
    
    # filter for genes with "RAS" in their name
    ras.de <- de[grepl("RAS", de$feature, ignore.case=TRUE) |
                   startsWith(de$feature, "MAP2K") | # MEK
                   startsWith(de$feature,"MAPK") | # MAPK1 is ERK
                   startsWith(de$feature,"PARP") | startsWith(de$feature,"BCL") |
                   startsWith(de$feature,"MCL") | startsWith(de$feature, "CASC") |
                   grepl("MTOR", de$feature, ignore.case=TRUE) | 
                   startsWith(de$feature, "RPS6") | startsWith(de$feature,"AKT"),]
    
    # dot plot
    sig.ras.de <- ras.de[ras.de$adj.P.Val < 0.05,]
    mean.ras <- plyr::ddply(sig.ras.de, .(feature), summarize,
                            Log2FC = mean(Log2FC, na.rm=TRUE))
    geneOrder <- na.omit(unique(mean.ras[order(mean.ras$Log2FC, decreasing=TRUE),]$feature))
    maxAbsNES <- max(abs(na.omit(ras.de[ras.de$feature %in% geneOrder,]$Log2FC)))
    maxLogFDR <- max(-log10(na.omit(ras.de[ras.de$feature %in% geneOrder,]$adj.P.Val)))
    ggplot2::ggplot(ras.de[ras.de$feature %in% geneOrder,], 
                    aes(x=Contrast, y=feature, color=Log2FC, size=-log10(adj.P.Val))) + 
      geom_point() + theme_classic(base_size=12) + scale_y_discrete(limits=geneOrder) +
      scale_size_continuous(limits=c(0,maxLogFDR)) +
      scale_color_gradient2(low="blue",high="red", mid="grey", limits=c(-maxAbsNES, maxAbsNES)) +
      geom_point(data = subset(ras.de[ras.de$feature %in% geneOrder,], adj.P.Val < 0.05), 
                 col = "black", stroke = 1.5, shape = 21) +
      labs(size="-Log(FDR)") + theme(axis.title=element_blank(), 
                                     axis.text.x=element_text(angle=45, vjust=1, hjust=1))
    ggplot2::ggsave(paste0(t, "_", omics,"_RAS-etc_diffexp_dotPlot.pdf"), 
                    width=ifelse(omics=="Global",2.75,4), 
                    height=ifelse(length(geneOrder)/4 > 49, 49, length(geneOrder)/4))
    ggplot2::ggsave(paste0(t, "_", omics,"_RAS-etc_diffexp_dotPlot_wider.pdf"), 
                    width=ifelse(omics=="Global",3,5), 
                    height=ifelse(length(geneOrder)/4 > 49, 49, length(geneOrder)/4))
    
    # too many features, so facet by family (e.g., RAS, MEK, etc.)
    ras.de$Family <- NA
    families <- c("RAS", "MAP2K", "MAPK", "PARP", "BCL", 
                  "CASC", "MCL", "MTOR", "RPS6", "AKT")
    facetPlot <- NULL
    for (fam in sort(families)) {
      ras.de[grepl(fam, ras.de$feature),]$Family <- fam
      sig.ras.de <- na.omit(ras.de[ras.de$Family == fam & ras.de$adj.P.Val < 0.05,])
      if (nrow(sig.ras.de) > 0) {
        mean.ras <- plyr::ddply(sig.ras.de, .(feature), summarize,
                                Log2FC = mean(Log2FC, na.rm=TRUE))
        geneOrder <- unique(mean.ras[order(mean.ras$Log2FC, decreasing=TRUE),]$feature)
        famPlot <- ggplot2::ggplot(na.omit(ras.de[ras.de$Family == fam & ras.de$feature %in% geneOrder,]), 
                                   aes(x=Contrast, y=reorder(feature, Log2FC), color=Log2FC, size=-log10(adj.P.Val))) + 
          geom_point() + theme_classic(base_size=12) + scale_y_discrete(limits=geneOrder) + 
          scale_size_continuous(limits=c(0,maxLogFDR)) +
          scale_color_gradient2(low="blue",high="red", mid="grey", limits=c(-maxAbsNES, maxAbsNES)) +
          geom_point(data = subset(na.omit(ras.de[ras.de$Family == fam & ras.de$feature %in% geneOrder,]), adj.P.Val < 0.05), 
                     col = "black", stroke = 1.5, shape = 21) +
          labs(size="-Log(FDR)", title=fam) + 
          theme(axis.title=element_blank(), 
                plot.title=element_text(hjust=0.5, face="bold"),
                axis.text.x=element_text(angle=45, vjust=1, hjust=1))
        if (is.null(facetPlot)) {
          facetPlot <- famPlot
        } else {
          facetPlot <- facetPlot + famPlot
        } 
      }
    }
    if (length(facetPlot) %% 3 == 0) {
      facetPlot <- facetPlot + plot_layout(nrow=3, guides="collect") 
    } else {
      facetPlot <- facetPlot + plot_layout(nrow=2, guides="collect") 
    }
    if (omics == "Phospho") {
      ggsave(paste0(t, "_", omics, "_RAS-etc_diffexp_patchDotPlot.pdf"), width = 9, height = 12)
      ggsave(paste0(t, "_", omics, "_RAS-etc_diffexp_patchDotPlot_tall.pdf"), width = 9, height = 49) 
      ggsave(paste0(t, "_", omics, "_RAS-etc_diffexp_patchDotPlot_wider.pdf"), width = 12, height = 12)
      ggsave(paste0(t, "_", omics, "_RAS-etc_diffexp_patchDotPlot_widerTall.pdf"), width = 12, height = 49)
    } else {
      ggsave(paste0(t, "_", omics, "_RAS-etc_diffexp_patchDotPlot.pdf"), width = 5.5, height = 8)
      ggsave(paste0(t, "_", omics, "_RAS-etc_diffexp_patchDotPlot_tall.pdf"), width = 5.5, height = 11)
      ggsave(paste0(t, "_", omics, "_RAS-etc_diffexp_patchDotPlot_wider.pdf"), width = 8, height = 8)
      ggsave(paste0(t, "_", omics, "_RAS-etc_diffexp_patchDotPlot_evenWider.pdf"), width = 9, height = 8)
      ggsave(paste0(t, "_", omics, "_RAS-etc_diffexp_patchDotPlot_widerTall.pdf"), width = 8, height = 11)
    }
    
    # redo with RAS separate
    families <- c("MAP2K", "MAPK", "PARP", "BCL", 
                  "CASC", "MCL", "MTOR", "RPS6", "AKT")
    facetPlot <- NULL
    for (fam in sort(families)) {
      ras.de[grepl(fam, ras.de$feature),]$Family <- fam
      sig.ras.de <- na.omit(ras.de[ras.de$Family == fam & ras.de$adj.P.Val < 0.05,])
      if (nrow(sig.ras.de) > 0) {
        mean.ras <- plyr::ddply(sig.ras.de, .(feature), summarize,
                                Log2FC = mean(Log2FC, na.rm=TRUE))
        geneOrder <- unique(mean.ras[order(mean.ras$Log2FC, decreasing=TRUE),]$feature)
        famPlot <- ggplot2::ggplot(na.omit(ras.de[ras.de$Family == fam & ras.de$feature %in% geneOrder,]), 
                                   aes(x=Contrast, y=reorder(feature, Log2FC), color=Log2FC, size=-log10(adj.P.Val))) + 
          geom_point() + theme_classic(base_size=12) + scale_y_discrete(limits=geneOrder) + 
          scale_size_continuous(limits=c(0,maxLogFDR)) +
          scale_color_gradient2(low="blue",high="red", mid="grey", limits=c(-maxAbsNES, maxAbsNES)) +
          geom_point(data = subset(na.omit(ras.de[ras.de$Family == fam & ras.de$feature %in% geneOrder,]), adj.P.Val < 0.05), 
                     col = "black", stroke = 1.5, shape = 21) +
          labs(size="-Log(FDR)", title=fam) + 
          theme(axis.title=element_blank(), 
                plot.title=element_text(hjust=0.5, face="bold"),
                axis.text.x=element_text(angle=45, vjust=1, hjust=1))
        if (is.null(facetPlot)) {
          facetPlot <- famPlot
        } else {
          facetPlot <- facetPlot + famPlot
        } 
      }
    }
    if (length(facetPlot) %% 3 == 0) {
      facetPlot <- facetPlot + plot_layout(nrow=3, guides="collect") 
    } else {
      facetPlot <- facetPlot + plot_layout(nrow=2, guides="collect") 
    }
    if (omics == "Phospho") {
      ggsave(paste0(t, "_", omics, "_MAPK_diffexp_patchDotPlot.pdf"), width = 12, height = 8)
      ggsave(paste0(t, "_", omics, "_MAPK_diffexp_patchDotPlot_tall.pdf"), width = 12, height = 49)
      ggsave(paste0(t, "_", omics, "_MAPK_diffexp_patchDotPlot_wider.pdf"), width = 14, height = 8)
      ggsave(paste0(t, "_", omics, "_MAPK_diffexp_patchDotPlot_widerTall.pdf"), width = 14, height = 49)
    } else {
      ggsave(paste0(t, "_", omics, "_MAPK_diffexp_patchDotPlot.pdf"), width = 7, height = 4.5)
      ggsave(paste0(t, "_", omics, "_MAPK_diffexp_patchDotPlot_wider.pdf"), width = 9, height = 4.5)
    }
    
    ## just RAS
    # filter for genes with "RAS" in their name
    ras.de <- de[grepl("RAS", de$feature, ignore.case=TRUE) ,]
    if (t == "Mix") {
      ras.de$Contrast2 <- 0
      ras.de[grepl("Gilt", ras.de$Contrast),]$Contrast2 <- 6
      ras.de$Contrast2 <- factor(ras.de$Contrast2, levels=c(0,6))
    } else {
      ras.de$Contrast2 <- ras.de$Contrast
    }
    
    # dot plot
    mean.ras <- plyr::ddply(ras.de, .(feature), summarize,
                            Log2FC = mean(Log2FC, na.rm=TRUE))
    geneOrder <- unique(mean.ras[order(mean.ras$Log2FC, decreasing=TRUE),]$feature)
    maxAbsNES <- max(abs(na.omit(ras.de[ras.de$feature %in% geneOrder,]$Log2FC)))
    maxLogFDR <- max(-log10(na.omit(ras.de[ras.de$feature %in% geneOrder,]$adj.P.Val)))
    ggplot2::ggplot(ras.de[ras.de$feature %in% geneOrder,], 
                    aes(x=Contrast2, y=feature, color=Log2FC, size=-log10(adj.P.Val))) + 
      geom_point() + theme_classic(base_size=12) + scale_y_discrete(limits=geneOrder) +
      scale_size_continuous(limits=c(0,maxLogFDR)) +
      scale_color_gradient2(low="blue",high="red", mid="grey", limits=c(-maxAbsNES, maxAbsNES)) +
      geom_point(data = subset(ras.de[ras.de$feature %in% geneOrder,], adj.P.Val < 0.05), 
                 col = "black", stroke = 1.5, shape = 21) +
      labs(size="-Log(FDR)", x="NRAS ASO + Gilteritinib Tx (hr)") + theme(axis.title.y=element_blank())
      #theme(axis.title=element_blank(), axis.text.x=element_text(angle=45, vjust=1, hjust=1))
    ggplot2::ggsave(paste0(t, "_", omics,"_RAS_diffexp_dotPlot.pdf"), 
                    width=ifelse(omics=="Global",2.75,4), 
                    height=ifelse(length(geneOrder)/4.6 > 49, 49, length(geneOrder)/4.6)+0.5)
    if(omics=="Global") {
      ggplot2::ggsave(paste0(t, "_", omics,"_RAS_diffexp_dotPlot_wider.pdf"), 
                    width=3, 
                    height=ifelse(length(geneOrder)/4 > 49, 49, length(geneOrder)/4)+0.5)
      ggplot2::ggsave(paste0(t, "_", omics,"_RAS_diffexp_dotPlot_width3.5.pdf"), 
                      width=3.5, 
                      height=ifelse(length(geneOrder)/4 > 49, 49, 0.5+length(geneOrder)/4)+0.5)
    }
  } 
}
saveRDS(all.de, "diffexp.rds")
all.diffexp <- data.table::rbindlist(all.de$Mix, use.names = TRUE, idcol="Omics")
write.csv(all.diffexp, "SupplementaryTable_Differential_expression.csv", row.names=FALSE)
all.de <- readRDS("diffexp.rds")

#### GSEA/KSEA dot plots of combo over time ####
# synIDs.mix <- list("Global" = list("NRAS ASO" = "syn64732761",
#                                    "NRAS ASO + Gilt" = "syn64732940"),
#                    "Phospho" = list("NRAS ASO" = "syn64732869",
#                                     "NRAS ASO + Gilt" = "syn64733043"))
data("STRINGv12")
kegg.info <- msigdbr::msigdbr(collection="C2",subcollection="CP:KEGG_LEGACY")
cellCycle <- na.omit(unique(kegg.info[kegg.info$gs_name=="KEGG_CELL_CYCLE",]$gene_symbol)) #125
mapkSignaling <- na.omit(unique(kegg.info[kegg.info$gs_description=="MAPK signaling pathway",]$gene_symbol)) #267
sphingo <- na.omit(unique(kegg.info[kegg.info$gs_name=="KEGG_SPHINGOLIPID_METABOLISM",]$gene_symbol)) # 39
oxphos <- na.omit(unique(kegg.info[kegg.info$gs_name=="KEGG_OXIDATIVE_PHOSPHORYLATION",]$gene_symbol)) # 132

synIDs.mix <- list("Global" = list("0" = "syn66728199",
                                   "6" = "syn66728346"),
                   "Phospho" = list("0" = "syn66726982",
                                    "6" = "syn66726993"))
all.synIDs <- list("Mix" = synIDs.mix)
network <- TRUE
all.ea <- list()
for (t in names(all.synIDs)) {
  synIDs <- all.synIDs[[t]]
  if (t %in% names(all.ea)) {
    time.ea <- all.ea[[t]]
  } else {
    time.ea <- list()
  }
  for (omics in names(synIDs)) {
    # load enrichment results
    if (omics %in% names(time.ea)) {
      ea <- time.ea[[omics]] 
    } else {
      ea <- data.frame()
      for (contrast in names(synIDs[[omics]])) {
        temp.ea <- read.csv(synapser::synGet(synIDs[[omics]][[contrast]])$path) # uses all timepoints as factor
        temp.ea$Contrast <- contrast
        ea <- rbind(ea, temp.ea)
      }
      ea$Contrast <- factor(ea$Contrast, levels=names(synIDs[[omics]]))
      ea$minusLogFDR <- 4
      ea[ea$FDR_q_value != 0,]$minusLogFDR <- -log10(ea[ea$FDR_q_value != 0,]$FDR_q_value)
      ea$sig <- FALSE
      ea[ea$p_value <= 0.05 & ea$FDR_q_value <= 0.25,]$sig <- TRUE
      time.ea[[omics]] <- ea 
    }
    
    # dot plot
    ea$Feature_set2 <- sub("KEGG_", "", ea$Feature_set)
    sig.ea <- ea[ea$sig,]
    if (nrow(sig.ea) > 0) {
      mean.ea <- plyr::ddply(sig.ea, .(Feature_set2), summarize,
                             NES = mean(NES, na.rm=TRUE))
      geneOrder <- na.omit(unique(mean.ea[order(mean.ea$NES, decreasing=TRUE),]$Feature_set2))
      maxAbsNES <- max(abs(na.omit(ea[ea$Feature_set2 %in% geneOrder,]$NES)))
      maxLogFDR <- max(na.omit(ea[ea$Feature_set2 %in% geneOrder,]$minusLogFDR))
      ggplot2::ggplot(ea[ea$Feature_set2 %in% geneOrder,], 
                      aes(x=Contrast, y=Feature_set2, color=NES, size=minusLogFDR)) + 
        geom_point() + theme_classic(base_size=12) + scale_y_discrete(limits=geneOrder) +
        scale_size_continuous(limits=c(0,maxLogFDR)) +
        scale_color_gradient2(low="blue",high="red", mid="grey", limits=c(-maxAbsNES, maxAbsNES)) +
        geom_point(data = subset(ea[ea$Feature_set2 %in% geneOrder,], sig), 
                   col = "black", stroke = 1.5, shape = 21) +
        labs(size="-Log(FDR)", x="NRAS ASO + Gilteritinib Tx (hr)") + theme(axis.title.y=element_blank(), legend.position="left")
      customHeight <- ifelse(omics=="Global", ifelse(length(geneOrder)/4 > 49, 49, length(geneOrder)/4),
                             ifelse(length(geneOrder)/4 < 3, 3, length(geneOrder)/4))
      ggplot2::ggsave(ifelse(omics=="Global", paste0(t, "_", omics,"_GSEA_KEGG_dotPlot.pdf"),
                             paste0(t, "_", omics,"_KSEA_dotPlot.pdf")), 
                      width=ifelse(omics=="Global",7,3), 
                      height=customHeight) 
      
      mean.ea2 <- mean.ea %>% slice_max(abs(NES), n=25)
      geneOrder <- na.omit(unique(mean.ea2[order(mean.ea2$NES, decreasing=TRUE),]$Feature_set2))
      maxAbsNES <- max(abs(na.omit(ea[ea$Feature_set2 %in% geneOrder,]$NES)))
      maxLogFDR <- max(na.omit(ea[ea$Feature_set2 %in% geneOrder,]$minusLogFDR))
      ggplot2::ggplot(ea[ea$Feature_set2 %in% geneOrder,], 
                      aes(x=Contrast, y=Feature_set2, color=NES, size=minusLogFDR)) + 
        geom_point() + theme_classic(base_size=12) + scale_y_discrete(limits=geneOrder, position="right") +
        scale_size_continuous(limits=c(0,maxLogFDR)) +
        scale_color_gradient2(low="blue",high="red", mid="grey", limits=c(-maxAbsNES, maxAbsNES)) +
        geom_point(data = subset(ea[ea$Feature_set2 %in% geneOrder,], sig), 
                   col = "black", stroke = 1.5, shape = 21) +
        labs(size="-Log(FDR)", x="NRAS ASO + Gilteritinib Tx (hr)") + theme(axis.title.y=element_blank(), legend.position="left")
      customHeight <- ifelse(omics=="Global", ifelse(length(geneOrder)/4 > 49, 49, length(geneOrder)/4),
                             ifelse(length(geneOrder)/4 < 3, 3, length(geneOrder)/4))
      ggplot2::ggsave(ifelse(omics=="Global", paste0(t, "_", omics,"_GSEA_KEGG_top25_dotPlot.pdf"),
                             paste0(t, "_", omics,"_KSEA_top25_dotPlot.pdf")), 
                      width=ifelse(omics=="Global",7,3), 
                      height=customHeight) 
      
      if (network & omics=="Phospho") {
        ksea.vert <- data.frame(geneOrder)
        colnames(ksea.vert) <- "Gene"
        ksea.vert$Group <- "Cell Cycle" # verify if KSEA is changed
        ksea.vert[grepl("MAP", ksea.vert$Gene) | grepl("RAF", ksea.vert$Gene),]$Group <- "MAPK"
        ksea.edges <- STRINGv12[STRINGv12$from %in% ksea.vert$Gene &
                                  STRINGv12$to %in% ksea.vert$Gene,]
        ksea.edges <- ksea.edges[ksea.edges$from != ksea.edges$to,]
        ksea.edges$combo <- NA
        for (edge in 1:nrow(ksea.edges)) {
          ksea.edges$combo[edge] <- paste0(sort(c(ksea.edges$from[edge], ksea.edges$to[edge])), collapse="_&_") 
        }
        ksea.edges <- ksea.edges[!duplicated(ksea.edges$combo),]
        ksea.vert <- ksea.vert[ksea.vert$Gene %in% c(ksea.edges$from, ksea.edges$to),]
        unique(sig.ras.ea$Feature_set)[!(unique(sig.ras.ea$Feature_set) %in% ksea.vert$Gene)]
        
        # turn data frame into igraph object
        temp.network <- igraph::graph_from_data_frame(
          d = ksea.edges, directed = FALSE,
          vertices = ksea.vert
        )
        RCy3::createNetworkFromIgraph(temp.network, title="exp21_sig_top25_KSEA_STRINGv12")
        
        sig.ea2 <- sig.ea[sig.ea$Feature_set %in% c(cellCycle, mapkSignaling),]
        mean.ea2 <- plyr::ddply(sig.ea2, .(Feature_set2), summarize,
                                NES = mean(NES, na.rm=TRUE))
        geneOrder <- na.omit(unique(mean.ea2[order(mean.ea2$NES, decreasing=TRUE),]$Feature_set2))
        maxAbsNES <- max(abs(na.omit(ea[ea$Feature_set2 %in% geneOrder,]$NES)))
        maxLogFDR <- max(na.omit(ea[ea$Feature_set2 %in% geneOrder,]$minusLogFDR))
        dot.df <- ea[ea$Feature_set2 %in% geneOrder,]
        dot.df$color <- rgb(81,174,174,maxColorValue = 255)
        dot.df[dot.df$Feature_set %in% mapkSignaling,]$color <- rgb(135,48,138,maxColorValue = 255)
        ggplot2::ggplot(dot.df, 
                        aes(x=Contrast, y=Feature_set2, color=NES, size=minusLogFDR)) + 
          geom_point() + theme_classic(base_size=12) + scale_y_discrete(limits=geneOrder, position="right") +
          scale_size_continuous(limits=c(0,maxLogFDR)) +
          scale_color_gradient2(low="blue",high="red", mid="grey", limits=c(-maxAbsNES, maxAbsNES)) +
          geom_point(data = subset(dot.df, sig), 
                     col = "black", stroke = 1.5, shape = 21) +
          labs(size="-Log(FDR)", x="NRAS ASO + Gilteritinib Tx (hr)") + 
          theme(axis.title.y=element_blank(), legend.position="left", axis.text.y=element_text(color=dot.df$color))
        customHeight <- ifelse(length(geneOrder)/4 < 3, 3, length(geneOrder)/4)
        ggplot2::ggsave(paste0(t, "_", omics,"_KSEA_cellCycleMAPK_dotPlot.pdf"), 
                        width=3, 
                        height=customHeight) 
        
        keptCellCycle <- cellCycle[cellCycle %in% geneOrder] # 3
        keptMapk <- mapkSignaling[mapkSignaling %in% geneOrder] # 36
        keptBoth <- keptCellCycle[keptCellCycle %in% keptMapk] # 0
        geneOrder <- c(sort(keptCellCycle), sort(keptMapk)) 
        textColors <- c(rep(rgb(81,174,174,maxColorValue = 255), length(keptCellCycle)),
                        rep(rgb(135,48,138,maxColorValue = 255), length(keptMapk)))
        ggplot2::ggplot(dot.df, 
                        aes(x=Contrast, y=Feature_set2, color=NES, size=minusLogFDR)) + 
          geom_point() + theme_classic(base_size=12) + scale_y_discrete(limits=geneOrder, position="right") +
          scale_size_continuous(limits=c(0,maxLogFDR)) +
          scale_color_gradient2(low="blue",high="red", mid="grey", limits=c(-maxAbsNES, maxAbsNES)) +
          geom_point(data = subset(dot.df, sig), 
                     col = "black", stroke = 1.5, shape = 21) +
          labs(size="-Log(FDR)", x="NRAS ASO + Gilteritinib Tx (hr)") + 
          theme(axis.title.y=element_blank(), legend.position="left", axis.text.y=element_text(color=textColors))
        customHeight <- ifelse(length(geneOrder)/5 < 3, 3, length(geneOrder)/5)
        ggplot2::ggsave(paste0(t, "_", omics,"_KSEA_cellCycleMAPK_dotPlot_v2.pdf"), 
                        width=3, 
                        height=customHeight)
        
        geneOrder <- c(sort(keptMapk), sort(keptCellCycle))
        textColors <- c(rep(rgb(135,48,138,maxColorValue = 255), length(keptMapk)),
                        rep(rgb(81,174,174,maxColorValue = 255), length(keptCellCycle)))
        ggplot2::ggplot(dot.df, 
                        aes(x=Contrast, y=Feature_set2, color=NES, size=minusLogFDR)) + 
          geom_point() + theme_classic(base_size=12) + scale_y_discrete(limits=geneOrder, position="right") +
          scale_size_continuous(limits=c(0,maxLogFDR)) +
          scale_color_gradient2(low="blue",high="red", mid="grey", limits=c(-maxAbsNES, maxAbsNES)) +
          geom_point(data = subset(dot.df, sig), 
                     col = "black", stroke = 1.5, shape = 21) +
          labs(size="-Log(FDR)", x="NRAS ASO + Gilteritinib Tx (hr)") + 
          theme(axis.title.y=element_blank(), legend.position="left", axis.text.y=element_text(color=textColors))
        customHeight <- ifelse(length(geneOrder)/5 < 3, 3, length(geneOrder)/5)
        ggplot2::ggsave(paste0(t, "_", omics,"_KSEA_cellCycleMAPK_dotPlot_v3.pdf"), 
                        width=3, 
                        height=customHeight)
        
        ksea.vert <- data.frame(geneOrder)
        colnames(ksea.vert) <- "Gene"
        ksea.vert$Group <- "Cell Cycle"
        ksea.vert[ksea.vert$Gene %in% mapkSignaling,]$Group <- "MAPK"
        ksea.edges <- STRINGv12[STRINGv12$from %in% ksea.vert$Gene &
                                  STRINGv12$to %in% ksea.vert$Gene,]
        ksea.edges <- ksea.edges[ksea.edges$from != ksea.edges$to,]
        ksea.edges$combo <- NA
        for (edge in 1:nrow(ksea.edges)) {
          ksea.edges$combo[edge] <- paste0(sort(c(ksea.edges$from[edge], ksea.edges$to[edge])), collapse="_&_") 
        }
        ksea.edges <- ksea.edges[!duplicated(ksea.edges$combo),]
        ksea.vert <- ksea.vert[ksea.vert$Gene %in% c(ksea.edges$from, ksea.edges$to),]
        unique(sig.ras.ea$Feature_set)[!(unique(sig.ras.ea$Feature_set) %in% ksea.vert$Gene)]
        
        # turn data frame into igraph object
        temp.network <- igraph::graph_from_data_frame(
          d = ksea.edges, directed = FALSE,
          vertices = ksea.vert
        )
        RCy3::createNetworkFromIgraph(temp.network, title="exp21_sig_cellCycleMAPK_KSEA_STRINGv12")
        
        # also separate out PI3K/AKT kinases
        keptPi3kAkt <- geneOrder[startsWith(geneOrder, "AKT") | startsWith(geneOrder, "STK")|startsWith(geneOrder, "PRK") |
                                   startsWith(geneOrder, "PAK") | startsWith(geneOrder, "RPS6K") |
                                   startsWith(geneOrder, "GSK3") | geneOrder %in% c("ATM","CHUK","NLK") | startsWith(geneOrder, "IKBK")] # from Sunil
        keptMAPK2 <- keptMapk[!(keptMapk %in% keptPi3kAkt)]
        keptCycle2 <- keptCellCycle[!(keptCellCycle %in% keptPi3kAkt)]
        geneOrder <- c(sort(keptMAPK2), sort(keptCycle2), sort(keptPi3kAkt))
        textColors <- c(rep(rgb(135,48,138,maxColorValue = 255), length(keptMAPK2)),
                        rep(rgb(81,174,174,maxColorValue = 255), length(keptCycle2)),
                        rep("#3C6EB6", length(keptPi3kAkt)))
        ggplot2::ggplot(dot.df, 
                        aes(x=Contrast, y=Feature_set2, color=NES, size=minusLogFDR)) + 
          geom_point() + theme_classic(base_size=12) + scale_y_discrete(limits=geneOrder, position="right") +
          scale_size_continuous(limits=c(0,maxLogFDR)) +
          scale_color_gradient2(low="blue",high="red", mid="grey", limits=c(-maxAbsNES, maxAbsNES)) +
          geom_point(data = subset(dot.df, sig), 
                     col = "black", stroke = 1.5, shape = 21) +
          labs(size="-Log(FDR)", x="NRAS ASO + Gilteritinib Tx (hr)") + 
          theme(axis.title.y=element_blank(), legend.position="left", axis.text.y=element_text(color=textColors))
        customHeight <- ifelse(length(geneOrder)/5 < 3, 3, length(geneOrder)/5)
        ggplot2::ggsave(paste0(t, "_", omics,"_KSEA_cellCycleMAPKAKT_dotPlot_v3.pdf"), 
                        width=3, 
                        height=customHeight)
        ggplot2::ggsave(paste0(t, "_", omics,"_KSEA_cellCycleMAPKAKT_dotPlot_v3_lessWide.pdf"), 
                        width=2.65, 
                        height=customHeight)
        
        sig.ea2 <- sig.ea[sig.ea$Feature_set %in% c(cellCycle,mapkSignaling,sphingo,oxphos),]
        mean.ea2 <- plyr::ddply(sig.ea2, .(Feature_set2), summarize,
                                NES = mean(NES, na.rm=TRUE))
        geneOrder <- na.omit(unique(mean.ea2[order(mean.ea2$NES, decreasing=TRUE),]$Feature_set2))
        maxAbsNES <- max(abs(na.omit(ea[ea$Feature_set2 %in% geneOrder,]$NES)))
        maxLogFDR <- max(na.omit(ea[ea$Feature_set2 %in% geneOrder,]$minusLogFDR))
        dot.df <- ea[ea$Feature_set2 %in% geneOrder,]
        dot.df$color <- rgb(81,174,174,maxColorValue = 255)
        dot.df[dot.df$Feature_set %in% mapkSignaling,]$color <- rgb(135,48,138,maxColorValue = 255)
        if (any(dot.df$Feature_set %in% sphingo)) {
          dot.df[dot.df$Feature_set %in% sphingo,]$color <- rgb(36,87,68, maxColorValue = 255)
        }
        if (any(dot.df$Feature_set %in% oxphos)) {
          dot.df[dot.df$Feature_set %in% oxphos,]$color <- rgb(73,51,142,maxColorValue = 255)
        }
        ggplot2::ggplot(dot.df, 
                        aes(x=Contrast, y=Feature_set2, color=NES, size=minusLogFDR)) + 
          geom_point() + theme_classic(base_size=12) + scale_y_discrete(limits=geneOrder, position="right") +
          scale_size_continuous(limits=c(0,maxLogFDR)) +
          scale_color_gradient2(low="blue",high="red", mid="grey", limits=c(-maxAbsNES, maxAbsNES)) +
          geom_point(data = subset(dot.df, sig), 
                     col = "black", stroke = 1.5, shape = 21) +
          labs(size="-Log(FDR)", x="NRAS ASO + Gilteritinib Tx (hr)") + 
          theme(axis.title.y=element_blank(), legend.position="left", axis.text.y=element_text(color=dot.df$color))
        customHeight <- ifelse(length(geneOrder)/4 < 3, 3, length(geneOrder)/4)
        ggplot2::ggsave(paste0(t, "_", omics,"_KSEA_cellCycleMAPKSphingoOXPHOS_dotPlot.pdf"), 
                        width=3, 
                        height=customHeight) 
        ksea.vert <- data.frame(geneOrder)
        colnames(ksea.vert) <- "Gene"
        ksea.vert$Group <- "Cell Cycle"
        ksea.vert[ksea.vert$Gene %in% mapkSignaling,]$Group <- "MAPK"
        if (any(ksea.vert$Gene %in% sphingo)) {
          ksea.vert[ksea.vert$Gene %in% sphingo,]$Group <- "Sphingolipid"
        }
        if (any(ksea.vert$Gene %in% oxphos)) {
          ksea.vert[ksea.vert$Gene %in% oxphos,]$Group <- "OXPHOS"
        }
          ksea.edges <- STRINGv12[STRINGv12$from %in% ksea.vert$Gene &
                                  STRINGv12$to %in% ksea.vert$Gene,]
        ksea.edges <- ksea.edges[ksea.edges$from != ksea.edges$to,]
        ksea.edges$combo <- NA
        for (edge in 1:nrow(ksea.edges)) {
          ksea.edges$combo[edge] <- paste0(sort(c(ksea.edges$from[edge], ksea.edges$to[edge])), collapse="_&_") 
        }
        ksea.edges <- ksea.edges[!duplicated(ksea.edges$combo),]
        ksea.vert <- ksea.vert[ksea.vert$Gene %in% c(ksea.edges$from, ksea.edges$to),]
        unique(sig.ras.ea$Feature_set)[!(unique(sig.ras.ea$Feature_set) %in% ksea.vert$Gene)]
        
        # turn data frame into igraph object
        temp.network <- igraph::graph_from_data_frame(
          d = ksea.edges, directed = FALSE,
          vertices = ksea.vert
        )
        RCy3::createNetworkFromIgraph(temp.network, title="exp21_sig_cellCycleMAPKSphingoOXPHOS_KSEA_STRINGv12")
        
      }
    }
    
    if (omics == "Phospho") {
      # filter for genes in western blot plus a few more Sunil wanted
      ras.ea <- ea[grepl("RAS", ea$Feature_set, ignore.case=TRUE) |
                     startsWith(ea$Feature_set, "MAP2K") | # MEK
                     startsWith(ea$Feature_set,"MAPK") | # MAPK1 is ERK
                     startsWith(ea$Feature_set,"PARP") | startsWith(ea$Feature_set,"BCL") |
                     startsWith(ea$Feature_set,"MCL") | startsWith(ea$Feature_set, "CASC") |
                     grepl("MTOR", ea$Feature_set, ignore.case=TRUE) |
                     startsWith(ea$Feature_set, "RPS6") | startsWith(ea$Feature_set,"AKT") |
                     ea$Feature_set %in% c("PRKCA", "PRKCB", "PRKCZ", "PIM3", "ARAF", "BRAF", "RAF1"),]
      
      # sig dot plot
      sig.ras.ea <- ras.ea[ras.ea$sig,]
      if (any(ras.ea$sig)) {
        mean.ras <- plyr::ddply(ras.ea[ras.ea$Feature_set %in% sig.ras.ea$Feature_set,], .(Feature_set), summarize,
                                NES = mean(NES, na.rm=TRUE),
                                maxNES = max(NES, na.rm=TRUE),
                                minNES = min(NES, na.rm=TRUE),
                                maxLogFDR = max(minusLogFDR, na.rm=TRUE))
        mean.ras$score <- mean.ras$maxLogFDR*mean.ras$maxNES
        geneOrder <- na.omit(unique(mean.ras[order(mean.ras$score, decreasing=TRUE),]$Feature_set))
        maxAbsNES <- max(abs(na.omit(ras.ea[ras.ea$Feature_set %in% geneOrder,]$NES)))
        maxLogFDR <- max(na.omit(ras.ea[ras.ea$Feature_set %in% geneOrder,]$minusLogFDR))
        
        if (t == "Mix") {
          ras.ea$Contrast2 <- 0
          ras.ea[grepl("6", ras.ea$Contrast),]$Contrast2 <- 6
          ras.ea$Contrast2 <- factor(ras.ea$Contrast2, levels=c(0,2,6,24))
        } else {
          ras.ea$Contrast2 <- ras.ea$Contrast
        }
        
        #pak::pak('r-lib/ragg')
        ggplot2::ggplot(ras.ea[ras.ea$Feature_set %in% geneOrder,],
                        aes(x=Contrast2, y=Feature_set, color=NES, size=minusLogFDR)) +
          geom_point() + theme_classic(base_size=12) + scale_y_discrete(limits=geneOrder, position="right") +
          scale_size_continuous(limits=c(0,maxLogFDR)) +
          scale_color_gradient2(low="blue",high="red", mid="grey", limits=c(-maxAbsNES, maxAbsNES)) +
          geom_point(data = subset(ras.ea[ras.ea$Feature_set %in% geneOrder,], sig),
                     col = "black", stroke = 1.5, shape = 21) +
          labs(size="-Log(FDR)", x="NRAS ASO + Gilteritinib Tx (hr)") + theme(axis.title.y=element_blank(),
                                                                              #theme(axis.title=element_blank(), axis.text.x=element_text(angle=45, vjust=1, hjust=1))
                                                                              legend.position = "left")
        ggplot2::ggsave(paste0(t, "_", omics,"_sig_RAS-etc_KSEA_dotPlot.pdf"),
                        width=2.75,
                        height=ifelse(length(geneOrder)/4 > 49, 49, length(geneOrder)/4))
        ggplot2::ggsave(paste0(t, "_", omics,"_sig_RAS-etc_KSEA_dotPlot_wider.pdf"),
                        width=3,
                        height=ifelse(length(geneOrder)/4.5 > 49, 49, length(geneOrder)/4.5))
        ggplot2::ggplot(ras.ea[ras.ea$Feature_set %in% geneOrder,],
                        aes(x=Contrast, y=Feature_set, color=NES, size=minusLogFDR)) +
          geom_point() + theme_classic(base_size=12) + scale_y_discrete(limits=geneOrder) +
          scale_size_continuous(limits=c(0,maxLogFDR)) + scale_x_discrete(limits=rev(names(synIDs[[omics]]))) +
          scale_color_gradient2(low="blue",high="red", mid="grey", limits=c(-maxAbsNES, maxAbsNES)) +
          geom_point(data = subset(ras.ea[ras.ea$Feature_set %in% geneOrder,], sig),
                     col = "black", stroke = 1.5, shape = 21) +
          labs(size="-Log(FDR)", x="NRAS ASO + Gilteritinib Tx (hr)") + theme(axis.title.y=element_blank()) +
          coord_flip() + theme(legend.position="bottom", legend.direction="horizontal")
        ggplot2::ggsave(paste0(t, "_", omics,"_sig_RAS-etc_KSEA_dotPlot_horizontal.pdf"),
                        width=ifelse(length(geneOrder)/4 > 49, 49, length(geneOrder)/3),
                        height=3)
        ggplot2::ggsave(paste0(t, "_", omics,"_sig_RAS-etc_KSEA_dotPlot_horizontal_shorter.pdf"),
                        width=ifelse(length(geneOrder)/4 > 49, 49, length(geneOrder)/3),
                        height=2.5)
        
        if (network) {
          ksea.vert <- data.frame(unique(sig.ras.ea$Feature_set))
          colnames(ksea.vert) <- "Gene"
          ksea.vert$Group <- "PI3K"
          ksea.vert[grepl("MAP", ksea.vert$Gene) | grepl("RAF", ksea.vert$Gene),]$Group <- "MAPK"
          ksea.edges <- STRINGv12[STRINGv12$from %in% ksea.vert$Gene &
                                    STRINGv12$to %in% ksea.vert$Gene,]
          # remove symmetry
          ksea.edges <- ksea.edges[ksea.edges$from != ksea.edges$to,]
          ksea.edges$combo <- NA
          for (edge in 1:nrow(ksea.edges)) {
            ksea.edges$combo[edge] <- paste0(sort(c(ksea.edges$from[edge], ksea.edges$to[edge])), collapse="_&_") 
          }
          ksea.edges <- ksea.edges[!duplicated(ksea.edges$combo),]
          ksea.vert <- ksea.vert[ksea.vert$Gene %in% c(ksea.edges$from, ksea.edges$to),]
          unique(sig.ras.ea$Feature_set)[!(unique(sig.ras.ea$Feature_set) %in% ksea.vert$Gene)]

          # turn data frame into igraph object
          temp.network <- igraph::graph_from_data_frame(
            d = ksea.edges, directed = FALSE,
            vertices = ksea.vert
          )
          RCy3::createNetworkFromIgraph(temp.network, title="exp21_sig_RAS-etc_KSEA_STRINGv12")
        }
      }
      
      # dot plot
      mean.ras <- plyr::ddply(ras.ea, .(Feature_set), summarize,
                              NES = mean(NES, na.rm=TRUE))
      geneOrder <- na.omit(unique(mean.ras[order(mean.ras$NES, decreasing=TRUE),]$Feature_set))
      maxAbsNES <- max(abs(na.omit(ras.ea[ras.ea$Feature_set %in% geneOrder,]$NES)))
      maxLogFDR <- max(na.omit(ras.ea[ras.ea$Feature_set %in% geneOrder,]$minusLogFDR))
      ggplot2::ggplot(ras.ea[ras.ea$Feature_set %in% geneOrder,],
                      aes(x=Contrast, y=Feature_set, color=NES, size=minusLogFDR)) +
        geom_point() + theme_classic(base_size=12) + scale_y_discrete(limits=geneOrder, position="right") +
        scale_size_continuous(limits=c(0,maxLogFDR)) +
        scale_color_gradient2(low="blue",high="red", mid="grey", limits=c(-maxAbsNES, maxAbsNES)) +
        geom_point(data = subset(ras.ea[ras.ea$Feature_set %in% geneOrder,], sig),
                   col = "black", stroke = 1.5, shape = 21) +
        labs(size="-Log(FDR)", x="NRAS ASO + Gilteritinib Tx (hr)") + theme(axis.title.y=element_blank())
      ggplot2::ggsave(paste0(t, "_", omics,"_pMAPK_KSEA_dotPlot.pdf"),
                      width=2.75,
                      height=ifelse(length(geneOrder)/4 > 49, 49, length(geneOrder)/4))
      ggplot2::ggsave(paste0(t, "_", omics,"_pMAPK_KSEA_dotPlot_wider.pdf"),
                      width=3,
                      height=ifelse(length(geneOrder)/4 > 49, 49, length(geneOrder)/4))
      ggplot2::ggsave(paste0(t, "_", omics,"_pMAPK_KSEA_dotPlot_width3.5.pdf"),
                      width=3.5,
                      height=ifelse(length(geneOrder)/4 > 49, 49, length(geneOrder)/4))
      
      # filter for just phospho in Western Blot
      ras.ea <- ea[startsWith(ea$Feature_set, "MAP2K") | # MEK
                     startsWith(ea$Feature_set,"MAPK") | # MAPK1 is ERK
                     grepl("MTOR", ea$Feature_set, ignore.case=TRUE) |
                     startsWith(ea$Feature_set, "RPS6") | startsWith(ea$Feature_set,"AKT"),]
      
      # sig dot plot
      sig.ras.ea <- ras.ea[ras.ea$sig,]
      if (any(ras.ea$sig)) {
        mean.ras <- plyr::ddply(sig.ras.ea, .(Feature_set), summarize,
                                NES = mean(NES, na.rm=TRUE))
        geneOrder <- na.omit(unique(mean.ras[order(mean.ras$NES, decreasing=TRUE),]$Feature_set))
        maxAbsNES <- max(abs(na.omit(ras.ea[ras.ea$Feature_set %in% geneOrder,]$NES)))
        maxLogFDR <- max(na.omit(ras.ea[ras.ea$Feature_set %in% geneOrder,]$minusLogFDR))
        ggplot2::ggplot(ras.ea[ras.ea$Feature_set %in% geneOrder,],
                        aes(x=Contrast, y=Feature_set, color=NES, size=minusLogFDR)) +
          geom_point() + theme_classic(base_size=12) + scale_y_discrete(limits=geneOrder, position="right") +
          scale_size_continuous(limits=c(0,maxLogFDR)) +
          scale_color_gradient2(low="blue",high="red", mid="grey", limits=c(-maxAbsNES, maxAbsNES)) +
          geom_point(data = subset(ras.ea[ras.ea$Feature_set %in% geneOrder,], sig),
                     col = "black", stroke = 1.5, shape = 21) +
          labs(size="-Log(FDR)", x="NRAS ASO + Gilteritinib Tx (hr)") + theme(axis.title.y=element_blank(), legend.position = "left")
        ggplot2::ggsave(paste0(t, "_", omics,"_sig_pMAPK_KSEA_dotPlot.pdf"),
                        width=2.75,
                        height=ifelse(length(geneOrder)/4 > 49, 49, length(geneOrder)/4))
        ggplot2::ggsave(paste0(t, "_", omics,"_sig_pMAPK_KSEA_dotPlot_wider.pdf"),
                        width=3,
                        height=ifelse(length(geneOrder)/4 > 49, 49, length(geneOrder)/4))
      }
      
      # dot plot
      mean.ras <- plyr::ddply(ras.ea, .(Feature_set), summarize,
                              NES = mean(NES, na.rm=TRUE))
      geneOrder <- na.omit(unique(mean.ras[order(mean.ras$NES, decreasing=TRUE),]$Feature_set))
      maxAbsNES <- max(abs(na.omit(ras.ea[ras.ea$Feature_set %in% geneOrder,]$NES)))
      maxLogFDR <- max(na.omit(ras.ea[ras.ea$Feature_set %in% geneOrder,]$minusLogFDR))
      ggplot2::ggplot(ras.ea[ras.ea$Feature_set %in% geneOrder,],
                      aes(x=Contrast, y=Feature_set, color=NES, size=minusLogFDR)) +
        geom_point() + theme_classic(base_size=12) + scale_y_discrete(limits=geneOrder) +
        scale_size_continuous(limits=c(0,maxLogFDR)) +
        scale_color_gradient2(low="blue",high="red", mid="grey", limits=c(-maxAbsNES, maxAbsNES)) +
        geom_point(data = subset(ras.ea[ras.ea$Feature_set %in% geneOrder,], sig),
                   col = "black", stroke = 1.5, shape = 21) +
        labs(size="-Log(FDR)") + theme(axis.title=element_blank(),
                                       axis.text.x=element_text(angle=45, vjust=1, hjust=1))
      ggplot2::ggsave(paste0(t, "_", omics,"_pMAPK_KSEA_dotPlot.pdf"),
                      width=2.75,
                      height=ifelse(length(geneOrder)/4 > 49, 49, length(geneOrder)/4))
      ggplot2::ggsave(paste0(t, "_", omics,"_pMAPK_KSEA_dotPlot_wider.pdf"),
                      width=3,
                      height=ifelse(length(geneOrder)/4 > 49, 49, length(geneOrder)/4))
      ggplot2::ggsave(paste0(t, "_", omics,"_pMAPK_KSEA_dotPlot_width3.5.pdf"),
                      width=3.5,
                      height=ifelse(length(geneOrder)/4 > 49, 49, length(geneOrder)/4))
    }
  } 
  all.ea[[t]] <- time.ea
}
saveRDS(all.ea, "enrichment.rds")
all.ea.df <- data.table::rbindlist(all.ea$Mix, use.names = TRUE, idcol="Omics")
write.csv(all.ea.df, "SupplementaryTable_Enrichment_analysis.csv", row.names=FALSE)
all.ea <- readRDS("enrichment.rds")

#### GSEA volcano plots for mixed contrasts ####
# load GSEA KEGG
kegg <- readxl::read_excel("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Exp21_NRAS_ASO_treated_patients/proteomics/analysis/unCorrectedLigPt/Figures/Figure3/Enrichment_results_SKJ_v2_BG.xlsx")
kegg <- kegg[kegg$Time == "Mix",]
kegg$minusLogFDR <- 4
kegg[kegg$FDR_q_value != 0,]$minusLogFDR <- -log10(kegg[kegg$FDR_q_value != 0,]$FDR_q_value)
kegg$minusLogP <- 4
kegg[kegg$p_value != 0,]$minusLogP <- -log10(kegg[kegg$p_value != 0,]$p_value)
kegg[is.na(kegg$Group) | kegg$Group == "NA",]$Group <- "Other"
groups <- c("Cell Cycle", "MAPK Signaling", "Metabolism", "Other")
kegg$Group <- factor(kegg$Group, levels=groups)
group.map <- dplyr::distinct(kegg[,c("Feature_set", "Group")])

kegg <- all.ea$Mix$Global
kegg$Label <- NA
kegg[kegg$Feature_set == "KEGG_SPHINGOLIPID_METABOLISM",]$Label <- "Sphingolipid\nMetabolism"
kegg[kegg$Feature_set == "KEGG_OXIDATIVE_PHOSPHORYLATION",]$Label <- "Oxidative\nPhosphorylation"
kegg <- merge(kegg, group.map, all.x=TRUE)
myColors <- c("lightseagreen", "magenta4", "#CCCCCC") # or is it "darkmagenta"? and "lightgrey"?
col2rgb(myColors)
# lightseagreen RGB: 32, 178, 170
# magenta4 RGB: 139, 0, 139
# greyish RGB: 204, 204, 204
myColors <- c(rgb(81,174,174,maxColorValue = 255),rgb(135,48,138,maxColorValue = 255),rgb(205,203,202,maxColorValue = 255))
maxAbsNES <- max(abs(kegg$NES))
kegg$minusLogP <- 4
kegg[kegg$p_value != 0,]$minusLogP <- -log10(kegg[kegg$p_value != 0,]$p_value)
maxLogP <- max(kegg$minusLogP)+0.5 # so points don't get cut off

volc.plots <- NULL
for (contrast in unique(kegg$Contrast)) {
  dot.df <- kegg[kegg$Contrast == contrast,]
  contr.plot <- ggplot2::ggplot(dot.df, aes(x=NES, y=minusLogP, color = Group)) + geom_point() +
    theme_classic(base_size=12) + labs(y="-Log(P-value)") +
    #geom_point(data = subset(dot.df, sig), col = "black", stroke = 1.5, shape = 21) +
    scale_color_manual(breaks = groups, values=c(myColors,"lightgrey")) + 
    scale_x_continuous(limits=c(-maxAbsNES, maxAbsNES)) + 
    scale_y_continuous(limits=c(0,maxLogP)) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color="darkgrey")+
    ggrepel::geom_label_repel(aes(label=ifelse(is.na(Label),"",Label)),#hjust=0,vjust=0, 
                              max.overlaps=Inf, box.padding=0.5) + labs(title=contrast) + 
    theme(plot.title=element_text(face="bold", hjust=0.5))
  if (is.null(volc.plots)) {
    volc.plots <- contr.plot
  } else {
    volc.plots <- volc.plots + contr.plot
  }
}
volc.plots <- volc.plots + plot_layout(nrow=1, guides="collect")
ggsave("GSEA_volcanoPlots_labels.pdf", volc.plots, width = 14, height = 2.5)
ggsave("GSEA_volcanoPlots_Wider.pdf", volc.plots, width = 20, height = 4)

volc.plots <- NULL
titles <- list("NRAS ASO" = "0 hr", "NRAS ASO + Gilt" = "6 hr")
for (contrast in unique(kegg$Contrast)) {
  dot.df <- kegg[kegg$Contrast == contrast & kegg$Group != "Other" & !is.na(kegg$Group),]
  contr.plot <- ggplot2::ggplot(dot.df, aes(x=NES, y=minusLogP, color = Group)) + geom_point(size = 4) +
    theme_classic(base_size=12) + labs(y="-Log(P-value)") +
    #geom_point(data = subset(dot.df, sig), col = "black", stroke = 1.5, shape = 21, size = 4) +
    scale_color_manual(breaks = groups, values=c(myColors)) +
    scale_x_continuous(limits=c(-maxAbsNES, maxAbsNES)) + 
    scale_y_continuous(limits=c(0,maxLogP)) + 
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color="darkgrey")+
    ggrepel::geom_text_repel(aes(label=ifelse(is.na(Label),"",Label)), 
                              max.overlaps=Inf, box.padding=0.5, 
                              color=ifelse(grepl("Phospho",dot.df$Label),rgb(73,51,142, maxColorValue = 255),rgb(36,87,68, maxColorValue = 255)),#"darkslateblue","darkslategray"),
                             fontface="bold",lineheight=0.7,
                              show.legend=FALSE) + 
    labs(title=titles[[contrast]]) + 
    theme(plot.title=element_text(face="bold", hjust=0.5))
  if (is.null(volc.plots)) {
    volc.plots <- contr.plot
  } else {
    volc.plots <- volc.plots + contr.plot
  }
}
volc.plots <- volc.plots + plot_layout(nrow=1, guides="collect")
ggsave("GSEA_volcanoPlots_groupsOnly.pdf", volc.plots, width = 7, height = 2.5)

volc.plots <- NULL
titles <- list("0" = "0 hr", "6" = "6 hr")
for (contrast in unique(kegg$Contrast)) {
  dot.df <- kegg[kegg$Contrast == contrast,]
  dot.df$color <- rgb(205,203,202,maxColorValue = 255)
  dot.df[dot.df$Feature_set=="KEGG_CELL_CYCLE",]$color <- rgb(81,174,174,maxColorValue = 255)
  dot.df[dot.df$Feature_set=="KEGG_MAPK_SIGNALING_PATHWAY",]$color <- rgb(135,48,138,maxColorValue = 255)
  dot.df[dot.df$Feature_set=="KEGG_OXIDATIVE_PHOSPHORYLATION",]$color <- rgb(73,51,142, maxColorValue = 255)
  dot.df[dot.df$Feature_set=="KEGG_SPHINGOLIPID_METABOLISM",]$color <- rgb(36,87,68, maxColorValue = 255)
  dot.df[dot.df$Feature_set=="KEGG_CELL_CYCLE",]$Label <- "Cell Cycle"
  dot.df[dot.df$Feature_set=="KEGG_MAPK_SIGNALING_PATHWAY",]$Label <- "MAPK Signaling"
  dot.df[dot.df$Feature_set=="KEGG_OXIDATIVE_PHOSPHORYLATION",]$Label <- "Oxidative\nPhosphorylation"
  dot.df[dot.df$Feature_set=="KEGG_SPHINGOLIPID_METABOLISM",]$Label <- "Sphingolipid\nMetabolism"
  
  contr.plot <- ggplot2::ggplot(dot.df, aes(x=NES, y=minusLogP)) + geom_point(size = 4, color=dot.df$color) +
    theme_classic(base_size=12) + labs(y="-Log(P-value)") +
    #geom_point(data = subset(dot.df, sig), col = "black", stroke = 1.5, shape = 21, size = 4) +
    scale_x_continuous(limits=c(-maxAbsNES, maxAbsNES)) + 
    scale_y_continuous(limits=c(0,maxLogP)) + 
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color="darkgrey")+
    ggrepel::geom_text_repel(aes(label=ifelse(is.na(Label),"",Label)), 
                             max.overlaps=Inf, box.padding=0.5, 
                             fontface="bold",lineheight=0.7,
                             show.legend=FALSE, color=dot.df$color) + 
    labs(title=titles[[contrast]], color="Pathway") + 
    theme(plot.title=element_text(face="bold", hjust=0.5), legend.position="none")
  if (is.null(volc.plots)) {
    volc.plots <- contr.plot
  } else {
    volc.plots <- volc.plots + contr.plot
  }
}
volc.plots <- volc.plots + plot_layout(nrow=1, guides="collect")
ggsave("GSEA_volcanoPlots_canonical.pdf", volc.plots, width = 7, height = 2.5)

volc.plots <- NULL
for (contrast in unique(kegg$Contrast)) {
  dot.df <- kegg[kegg$Contrast == contrast,]
  dot.df$color <- rgb(205,203,202,maxColorValue = 255)
  dot.df$Label <- NA
  if (nrow(dot.df[dot.df$Feature_set=="KEGG_CELL_CYCLE" & dot.df$sig,]) > 0) {
    dot.df[dot.df$Feature_set=="KEGG_CELL_CYCLE" & dot.df$sig,]$color <- rgb(81,174,174,maxColorValue = 255)
    dot.df[dot.df$Feature_set=="KEGG_CELL_CYCLE" & dot.df$sig,]$Label <- "Cell Cycle"
  }
  if (nrow(dot.df[dot.df$Feature_set=="KEGG_MAPK_SIGNALING_PATHWAY" & dot.df$sig,]) > 0) {
    dot.df[dot.df$Feature_set=="KEGG_MAPK_SIGNALING_PATHWAY" & dot.df$sig,]$color <- rgb(135,48,138,maxColorValue = 255)
    dot.df[dot.df$Feature_set=="KEGG_MAPK_SIGNALING_PATHWAY" & dot.df$sig,]$Label <- "MAPK Signaling"
  }
  if (nrow(dot.df[dot.df$Feature_set=="KEGG_OXIDATIVE_PHOSPHORYLATION" & dot.df$sig,]) > 0) {
    dot.df[dot.df$Feature_set=="KEGG_OXIDATIVE_PHOSPHORYLATION" & dot.df$sig,]$color <- rgb(73,51,142, maxColorValue = 255)
    dot.df[dot.df$Feature_set=="KEGG_OXIDATIVE_PHOSPHORYLATION" & dot.df$sig,]$Label <- "Oxidative\nPhosphorylation"
  }
  if (nrow(dot.df[dot.df$Feature_set=="KEGG_SPHINGOLIPID_METABOLISM" & dot.df$sig,]) > 0) {
    dot.df[dot.df$Feature_set=="KEGG_SPHINGOLIPID_METABOLISM" & dot.df$sig,]$color <- rgb(36,87,68, maxColorValue = 255)
    dot.df[dot.df$Feature_set=="KEGG_SPHINGOLIPID_METABOLISM" & dot.df$sig,]$Label <- "Sphingolipid\nMetabolism"
  }
  
  contr.plot <- ggplot2::ggplot(dot.df, aes(x=NES, y=minusLogP)) + geom_point(size = 4, color=dot.df$color) +
    theme_classic(base_size=12) + labs(y="-Log(P-value)") +
    #geom_point(data = subset(dot.df, sig), col = "black", stroke = 1.5, shape = 21, size = 4) +
    scale_x_continuous(limits=c(-maxAbsNES, maxAbsNES)) + 
    scale_y_continuous(limits=c(0,maxLogP)) + 
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color="darkgrey")+
    ggrepel::geom_text_repel(aes(label=ifelse(is.na(Label),"",Label)), 
                             max.overlaps=Inf, box.padding=1, 
                             fontface="bold",lineheight=0.7,
                             show.legend=FALSE, color=dot.df$color) + 
    labs(title=titles[[contrast]], color="Pathway") +
    theme(plot.title=element_text(face="bold", hjust=0.5), legend.position="none")
  if (is.null(volc.plots)) {
    volc.plots <- contr.plot
  } else {
    volc.plots <- volc.plots + contr.plot
  }
}
volc.plots <- volc.plots + plot_layout(nrow=1, guides="collect")
ggsave("GSEA_volcanoPlots_canonicalSig.pdf", volc.plots, width = 7, height = 2.5)

#### DMEA ####
### drug correlations
synIDs.mix <- list("Global" = "syn66726959", # only want to look at "0" aka NRAS ASO vs. CTRL ASO
                   "Phospho" = "syn66726973")
all.synIDs <- list("Mix" = synIDs.mix)

drug.info <- read.csv("~/OneDrive - PNNL/Documents/PTRC2/BeatAML_single_drug_moa.csv",
                      stringsAsFactors = FALSE, fileEncoding = "latin1")
drug.info <- drug.info[,c("Drug","moa")]
drug.info[drug.info$Drug == "Ralimetinib (LY2228820)",]$moa <- "p38 MAPK inhibitor"
drug.info[drug.info$Drug == "Nilotinib",]$moa <- "Abl kinase inhibitor"
drug.info[drug.info$Drug == "AT-101",]$moa <- "BCL inhibitor"
drug.info[is.na(drug.info$moa),]$moa <- "Other"
drug.info$Drug <- sub(" [(].*", "", drug.info$Drug) # shorten drug names for plot
drug.info[drug.info$Drug == "NF-kB Activation Inhibitor",]$Drug <- "NFkB Inhibitor"
drug.info$Target <- "Other"
drug.info$Target <- gsub(" inhibitor","", drug.info$moa)
drug.info[grepl("FLT3",drug.info$Target),]$Target <- "FLT3"
drug.info[grepl(",",drug.info$Target),]$Target <- "Other"
drug.info[drug.info$Target == "Aurora kinase",]$Target <- "AURK"
drug.info[drug.info$Target == "PDGFR tyrosine kinase receptor",]$Target <- "PDGFR"
drug.info[grepl("NFkB",drug.info$Target),]$Target <- "NFkB"
drug.info[grepl("TGF",drug.info$Target),]$Target <- "TGFB"
drug.info[grepl("cFMS",drug.info$Target),]$Target <- "Other"
drug.info[grepl("ALK",drug.info$Target),]$Target <- "ALK"
drug.info[grepl("BTK",drug.info$Target),]$Target <- "BTK"
drug.info[drug.info$Target %in% c("survivin","insulin growth factor receptor","tyrosine kinase","glycogen synthase kinase","oxidative stress inducer",
                                  "antineoplastic agent","gamma secretase","anticancer agent","REV-ERB agonist","DNA synthesis","exportin antagonist","lipoxygenase","carnitine palmitoyltransferase","mitochondrial respiratory chain"),]$Target <- "Other"
drug.info[drug.info$Target == "DNA methyltransferase",]$Target <- "DNMT"
drug.info[drug.info$Target == "bromodomain",]$Target <- "BRD"
drug.info[drug.info$Target == "ribonucleotide reductase",]$Target <- "RNR"
drug.info$Target <- factor(drug.info$Target, levels=c(unique(drug.info[drug.info$Target != "Other",]$Target), "Other"))

all.drug.corr <- list()
for (t in names(all.synIDs)) {
  synIDs <- all.synIDs[[t]]
  time.drug.corr <- all.drug.corr[[t]]
  if (is.null(time.drug.corr)) {
    time.drug.corr <- list()
  }
  for (omics in names(synIDs)) {
    # load diffexp
    if (omics %in% names(all.drug.corr)) {
      drug.corr <- time.drug.corr[[omics]]
    } else {
      drug.corr <- read.csv(synapser::synGet(synIDs[[omics]])$path) # uses all timepoints as factor
      drug.corr$Contrast <- "NRAS ASO"
      time.drug.corr[[omics]] <- drug.corr
    }
  }
  all.drug.corr[[t]] <- time.drug.corr
  
  # filter for significantly correlated drugs
  time.drug.corr <- data.table::rbindlist(time.drug.corr, use.names = TRUE, idcol = "Omics")
  time.drug.corr$Drug <- sub(" [(].*", "", time.drug.corr$Drug) # shorten drug names for plot
  dot.df0 <- na.omit(merge(time.drug.corr, drug.info, by="Drug", all.x=TRUE))
  
  sig.drug.corr <- dot.df0[dot.df0$Pearson.q < 0.05,]
  mean.drug <- plyr::ddply(sig.drug.corr, .(Drug), summarize,
                           Pearson.est = mean(Pearson.est, na.rm=TRUE),
                           N_sig=length(unique(Omics)))
  mean.drug <- mean.drug[mean.drug$N_sig==2,] # 35
  geneOrder <- na.omit(unique(mean.drug[order(mean.drug$Pearson.est, decreasing=TRUE),]$Drug)) # 35
  dot.df <- dot.df0[dot.df0$Drug %in% geneOrder,]
  dot.df$Target <- as.character(dot.df$Target)
  moas <- unique(dot.df$Target)
  dot.df$Target <- factor(dot.df$Target, levels=c(moas[moas!="Other"],"Other"))
  
  # require drug to be sig corr in both and
  # only fill in targets with more than 1 drug in list
  moa.df <- plyr::ddply(dot.df, .(Target), summarize,
                        nDrugs = length(unique(Drug)))
  dupMOAs <- as.character(moa.df[moa.df$nDrugs>1,]$Target)
  dupMOAs <- dupMOAs[dupMOAs != "Other"] # 5
  dot.df$Target2 <- "Other"
  dot.df[dot.df$Target %in% dupMOAs,]$Target2 <- dot.df[dot.df$Target %in% dupMOAs,]$Target
  moaOrder <- unique(dot.df[order(dot.df$Pearson.est),]$Target2)
  moaOrder <- c(moaOrder[moaOrder!="Other"],"Other")
  dot.df$Omics2 <- "Phosphoproteomics"
  dot.df[dot.df$Omics=="Global",]$Omics2 <- "Proteomics"
  dot.df$Omics2 <- factor(dot.df$Omics2, levels=c("Proteomics","Phosphoproteomics"))
  ggplot2::ggplot(dot.df, aes(x=Pearson.est, y=Drug, fill=Target2)) + 
    geom_bar(stat="identity") + facet_wrap(.~Omics2) + theme_classic(base_size=12) + 
    scale_fill_manual(values=c(scales::hue_pal()(length(dupMOAs)),rgb(205,203,202,maxColorValue = 255)), breaks=moaOrder) +
    scale_y_discrete(limits=geneOrder, position="right") + scale_x_continuous(limits=c(-1,1)) +
    labs(x="Pearson Correlation", fill="Target") + theme(axis.title.y=element_blank(), legend.position="bottom", legend.direction="horizontal")
  ggplot2::ggsave(paste0(t,"_sigCorrDrugs_dotPlot_targets.pdf"), 
                  width=4, 
                  height=ifelse(length(geneOrder)/5 > 49, 49, length(geneOrder)/5))
  ggplot2::ggsave(paste0(t,"_sigCorrDrugs_dotPlot_targets_wider.pdf"), 
                  width=5, 
                  height=ifelse(length(geneOrder)/5 > 49, 49, length(geneOrder)/5))
  
  dot.df$Direction <- "More toxic after NRAS ASO"
  dot.df[dot.df$Pearson.est>0,]$Direction <- "Less toxic after NRAS ASO"
  ggplot2::ggplot(dot.df[dot.df$Omics=="Global",], aes(y=Pearson.est, x=Drug, fill=Direction)) + 
    geom_bar(stat="identity", colour="black") + theme_classic(base_size=12) + 
    scale_x_discrete(limits=geneOrder) + #scale_y_continuous(limits=c(-1,1)) +
    scale_fill_manual(values=c(rgb(237,181,174, maxColorValue=255), rgb(189,203,242, maxColorValue=255)), 
                      breaks=c("Less toxic after NRAS ASO","More toxic after NRAS ASO"))+
    labs(y="Pearson r", fill="") + theme(axis.title.x=element_blank(), axis.line.x=element_blank(),
                                         axis.ticks.x = element_blank(),
                                                    legend.position=c(0.25,0.2),legend.direction="vertical",
                                         legend.background=element_blank(),
                                         axis.text.x = element_text(angle=45, hjust=1, vjust=1))
  ggplot2::ggsave(paste0(t,"_sigCorrDrugs_barPlot_global.pdf"), 
                  width=7, 
                  height=3)
  
  dot.df.all <- merge(dot.df0, dot.df, all.x=TRUE, by=c("Drug", "Omics", 
                                                        "Pearson.est", "Pearson.p",
                                                        "Pearson.q","Spearman.est", 
                                                        "Spearman.p","Spearman.q",
                                                        "Slope","Intercept", "R.squared",
                                                        "Rank.slope", "Rank.intercept",
                                                        "Rank.R.squared", "Contrast", 
                                                        "moa", "Target", "N"))
  targetColors <- c(scales::hue_pal()(length(dupMOAs)),rgb(205,203,202,maxColorValue = 255))
  dot.df.all$color <- rgb(205,203,202,maxColorValue = 255)
  for (i in 1:length(moaOrder)) {
    dot.df.all[dot.df.all$Target2 == moaOrder[i],]$color <- targetColors[i]
  }
  maxAbsNES <- max(abs(dot.df.all$Pearson.est))
  maxLogP <- max(-log10(dot.df.all$Pearson.p))
  dot.df.all[is.na(dot.df.all$Target2),]$Target2 <- "Other"
  ggplot2::ggplot(dot.df.all, aes(x=Pearson.est, y=-log10(Pearson.p), color=Target2)) + 
    geom_point() + facet_wrap(.~Omics) + theme_classic(base_size=12) + 
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color="darkgrey")+
    scale_x_continuous(limits=c(-1, 1)) + 
    scale_y_continuous(limits=c(0,maxLogP)) + 
    geom_point(data = subset(dot.df.all, Pearson.q<=0.05), col = "black", stroke = 1.5, shape = 21) +
    ggrepel::geom_text_repel(aes(label=ifelse(Target2 != "Other", ifelse(Pearson.q<=0.05,Drug, ""), "")), 
                             max.overlaps=Inf, 
                             box.padding=0.3, 
                             fontface="bold",lineheight=0.7,
                             show.legend=FALSE, color=dot.df.all$color) + 
    scale_color_manual(values=c(scales::hue_pal()(length(dupMOAs)),rgb(205,203,202,maxColorValue = 255)), breaks=moaOrder) +
    labs(x="Pearson Correlation", y="-Log(P-value)", color="Target") + theme(legend.position="bottom", legend.direction="horizontal")
  ggplot2::ggsave(paste0(t,"_corrDrugs_volcanoPlot_targets.pdf"), 
                  width=8, 
                  height=4)
  ggplot2::ggsave(paste0(t,"_corrDrugs_volcanoPlot_targets_lessWide.pdf"), 
                  width=5, 
                  height=4)
  
  ggplot2::ggplot(dot.df.all, aes(x=Pearson.est, y=-log10(Pearson.p), color=Target2)) + 
    geom_point() + facet_wrap(.~Omics) + theme_classic(base_size=12) + 
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color="darkgrey")+
    scale_x_continuous(limits=c(-1, 1)) + 
    scale_y_continuous(limits=c(0,maxLogP)) + 
    geom_point(data = subset(dot.df.all, Pearson.q<=0.05), col = "black", stroke = 1.5, shape = 21) +
    ggrepel::geom_text_repel(aes(label=ifelse(abs(Pearson.est) > 0.5, ifelse(Pearson.q<=0.05,Drug, ""), "")), 
                             max.overlaps=Inf, 
                             box.padding=0.3, 
                             fontface="bold",lineheight=0.7,
                             show.legend=FALSE, color=dot.df.all$color) + 
    scale_color_manual(values=c(scales::hue_pal()(length(dupMOAs)),rgb(205,203,202,maxColorValue = 255)), breaks=moaOrder) +
    labs(x="Pearson Correlation", y="-Log(P-value)", color="Target") + theme(legend.position="bottom", legend.direction="horizontal")
  ggplot2::ggsave(paste0(t,"_corrDrugs_volcanoPlot_targets_labelPearson0.5.pdf"), 
                  width=8, 
                  height=4)
  ggplot2::ggsave(paste0(t,"_corrDrugs_volcanoPlot_targets_lessWide_labelPearson0.5.pdf"), 
                  width=5, 
                  height=4)
  
  drugsToLabel <- unique(dot.df[dot.df$Pearson.q <= 0.05,]$Drug) # 34
  ggplot2::ggplot(dot.df.all, aes(x=Pearson.est, y=-log10(Pearson.p), color=Target2)) + 
    geom_point() + facet_wrap(.~Omics) + theme_classic(base_size=12) + 
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color="darkgrey")+
    scale_x_continuous(limits=c(-1, 1)) + 
    scale_y_continuous(limits=c(0,maxLogP)) + 
    geom_point(data = subset(dot.df.all, Pearson.q<=0.05), col = "black", stroke = 1.5, shape = 21) +
    ggrepel::geom_text_repel(aes(label=ifelse(Drug %in% drugsToLabel, Drug, "")), 
                             max.overlaps=Inf, 
                             box.padding=0.3, 
                             fontface="bold",lineheight=0.7,
                             show.legend=FALSE, color=dot.df.all$color) + 
    scale_color_manual(values=c(scales::hue_pal()(length(dupMOAs)),rgb(205,203,202,maxColorValue = 255)), breaks=moaOrder) +
    labs(x="Pearson Correlation", y="-Log(P-value)", color="Target") + theme(legend.position="bottom", legend.direction="horizontal")
  ggplot2::ggsave(paste0(t,"_corrDrugs_volcanoPlot_targets_sig2Omics.pdf"), 
                  width=8, 
                  height=4)
  ggplot2::ggsave(paste0(t,"_corrDrugs_volcanoPlot_targets_lessWide_sig2Omics.pdf"), 
                  width=5, 
                  height=4)
  
  drugsToLabel <- unique(c(dot.df[dot.df$Pearson.q <= 0.05 & dot.df$Target2 != "Other",]$Drug,
                    dot.df[dot.df$Pearson.est == max(dot.df$Pearson.est),]$Drug,
                    dot.df[dot.df$Pearson.est == min(dot.df$Pearson.est),]$Drug)) # 18
  dot.df.all$Significance <- "q > 0.05"
  dot.df.all[dot.df.all$Pearson.q <= 0.05,]$Significance <- "q < 0.05" # none == 0.05
  dot.df.all$Significance <- factor(dot.df.all$Significance, levels=c("q < 0.05", "q > 0.05"))
  dot.df.all[dot.df.all$Omics=="Global",]$Omics2 <- "Proteomics"
  dot.df.all[dot.df.all$Omics=="Phospho",]$Omics2 <- "Phosphoproteomics"
  ggplot2::ggplot(dot.df.all, aes(x=Pearson.est, y=-log10(Pearson.p))) + 
    geom_point(aes(fill=Target2, shape=Target2, linewidth=Significance, colour=Significance), size=3) + facet_wrap(.~Omics2) + theme_classic(base_size=12) + 
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color="darkgrey")+
    scale_x_continuous(limits=c(-1, 1)) + scale_linewidth_manual(values=c(1,0), breaks=c("q < 0.05", "q > 0.05"))+
    scale_color_manual(values=c("black", rgb(205,203,202,maxColorValue = 255)), breaks=c("q < 0.05", "q > 0.05"))+
    scale_y_continuous(limits=c(0,maxLogP)) + scale_shape_manual(values=c(seq(21,25),21), breaks=moaOrder) +
    #geom_point(data = subset(dot.df.all, Pearson.q<=0.05), col = "black", stroke = 1.5, shape = 21) +
    ggrepel::geom_text_repel(aes(label=ifelse(Drug %in% drugsToLabel, Drug, "")), 
                             max.overlaps=Inf, 
                             box.padding=0.5, 
                             fontface="bold",lineheight=0.7, size=3,
                             show.legend=FALSE, color=dot.df.all$color) + 
    scale_fill_manual(values=c(scales::hue_pal()(length(dupMOAs)),rgb(205,203,202,maxColorValue = 255)), breaks=moaOrder) +
    labs(x="Pearson Correlation", y="-Log(P-value)", fill="Target", shape="Target") + theme(legend.position="bottom", legend.direction="horizontal")
  ggplot2::ggsave(paste0(t,"_corrDrugs_volcanoPlot_targets_topAndSigSharedTargets_size3andShapes.pdf"), 
                  width=8, 
                  height=4)
  ggplot2::ggsave(paste0(t,"_corrDrugs_volcanoPlot_targets_lessWide_topAndSigSharedTargets_size3andShapes.pdf"), 
                  width=5, 
                  height=4)
  
  ggplot2::ggplot(dot.df.all[dot.df.all$Omics=="Global",], aes(x=Pearson.est, y=-log10(Pearson.p))) + 
    geom_point(aes(fill=Target2, shape=Target2, linewidth=Significance, colour=Significance), size=3) + theme_classic(base_size=12) + 
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color="darkgrey")+
    scale_x_continuous(limits=c(-1, 1)) + scale_linewidth_manual(values=c(1,0), breaks=c("q < 0.05", "q > 0.05"))+
    scale_color_manual(values=c("black", rgb(205,203,202,maxColorValue = 255)), breaks=c("q < 0.05", "q > 0.05"))+
    scale_y_continuous(limits=c(0,maxLogP)) + scale_shape_manual(values=c(seq(21,25),21), breaks=moaOrder) +
    #geom_point(data = subset(dot.df.all, Pearson.q<=0.05), col = "black", stroke = 1.5, shape = 21) +
    ggrepel::geom_text_repel(aes(label=ifelse(Drug %in% drugsToLabel, Drug, "")), 
                             max.overlaps=Inf, 
                             box.padding=0.5,
                             fontface="bold",lineheight=0.7, size=3,
                             show.legend=FALSE, color=dot.df.all[dot.df.all$Omics=="Global",]$color) + 
    scale_fill_manual(values=c(scales::hue_pal()(length(dupMOAs)),rgb(205,203,202,maxColorValue = 255)), breaks=moaOrder) +
    labs(x="Pearson Correlation", y="-Log(P-value)", fill="Target", shape="Target") + theme(legend.position="bottom"#, legend.direction="horizontal"
                                                                                            )
  ggplot2::ggsave(paste0(t,"_corrDrugs_volcanoPlot_targets_topAndSigSharedTargets_size3andShapes_global.pdf"), 
                  width=6, 
                  height=4)
  ggplot2::ggsave(paste0(t,"_corrDrugs_volcanoPlot_targets_lessWide_topAndSigSharedTargets_size3andShapes_global.pdf"), 
                  width=5, 
                  height=4)
  
  ggplot2::ggplot(dot.df.all, aes(x=Pearson.est, y=-log10(Pearson.p), color=Target2)) + 
    geom_point() + facet_grid(Omics~.) + theme_classic(base_size=12) + 
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color="darkgrey")+
    scale_x_continuous(limits=c(-1, 1)) + 
    scale_y_continuous(limits=c(0,maxLogP)) + 
    geom_point(data = subset(dot.df.all, Pearson.q<=0.05), col = "black", stroke = 1.5, shape = 21) +
    ggrepel::geom_text_repel(aes(label=ifelse(Drug %in% drugsToLabel, Drug, "")), 
                             max.overlaps=Inf, 
                             box.padding=0.5, 
                             fontface="bold",lineheight=0.7,
                             show.legend=FALSE, color=dot.df.all$color) + 
    scale_color_manual(values=c(scales::hue_pal()(length(dupMOAs)),rgb(205,203,202,maxColorValue = 255)), breaks=moaOrder) +
    labs(x="Pearson Correlation", y="-Log(P-value)", color="Target")# + theme(legend.position="bottom", legend.direction="horizontal")
  ggplot2::ggsave(paste0(t,"_corrDrugs_volcanoPlot_targets_topAndSigSharedTargets_vertical.pdf"), 
                  width=4, 
                  height=6)
  
  ggplot2::ggplot(dot.df.all[dot.df.all$Omics=="Global",], 
                  aes(x=Pearson.est, y=-log10(Pearson.p), color=Target2)) + 
    geom_point() + theme_classic(base_size=12) + 
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color="darkgrey")+
    scale_x_continuous(limits=c(-1, 1)) + 
    scale_y_continuous(limits=c(0,maxLogP)) + 
    geom_point(data = subset(dot.df.all[dot.df.all$Omics=="Global",], Pearson.q<=0.05), col = "black", stroke = 1.5, shape = 21) +
    ggrepel::geom_text_repel(aes(label=ifelse(Drug %in% drugsToLabel, Drug, "")), 
                             max.overlaps=Inf, 
                             #box.padding=0.5, 
                             fontface="bold",lineheight=0.7,
                             show.legend=FALSE, color=dot.df.all[dot.df.all$Omics=="Global",]$color) + 
    scale_color_manual(values=c(scales::hue_pal()(length(dupMOAs)),rgb(205,203,202,maxColorValue = 255)), breaks=moaOrder) +
    labs(x="Pearson Correlation", y="-Log(P-value)", color="Target")# + theme(legend.position="bottom", legend.direction="horizontal")
  ggplot2::ggsave(paste0(t,"_corrDrugs_volcanoPlot_targets_topAndSigSharedTargets_global.pdf"), 
                  width=6, 
                  height=4)
  
  ggplot2::ggplot(dot.df.all[dot.df.all$Omics=="Global",], 
                  aes(x=Pearson.est, y=-log10(Pearson.p), color=Target2)) + 
    geom_point() + theme_classic(base_size=12) + 
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color="darkgrey")+
    scale_x_continuous(limits=c(-1, 1)) + 
    scale_y_continuous(limits=c(0,maxLogP)) + 
    geom_point(data = subset(dot.df.all[dot.df.all$Omics=="Global",], Pearson.q<=0.05), col = "black", stroke = 1.5, shape = 21) +
    ggrepel::geom_text_repel(aes(label=ifelse(Drug %in% drugsToLabel, Drug, "")), 
                             max.overlaps=Inf, 
                             #box.padding=0.5, 
                             fontface="bold",lineheight=0.7,
                             show.legend=FALSE, color=dot.df.all[dot.df.all$Omics=="Global",]$color) + 
    scale_color_manual(values=c(scales::hue_pal()(length(dupMOAs)),rgb(205,203,202,maxColorValue = 255)), breaks=moaOrder) +
    labs(x="Pearson Correlation", y="-Log(P-value)", color="Target") + theme(legend.position="bottom", legend.direction="horizontal")
  ggplot2::ggsave(paste0(t,"_corrDrugs_volcanoPlot_targets_topAndSigSharedTargets_global_horizLegend.pdf"), 
                  width=6, 
                  height=4)
}
saveRDS(all.drug.corr, "drugCorr.rds")
all.drug.df <- data.table::rbindlist(all.drug.corr$Mix, use.names = TRUE, idcol="Omics")
write.csv(all.drug.df, "SupplementaryTable_Drug_sensitivity_correlations.csv", row.names=FALSE)
all.drug.corr <- readRDS("drugCorr.rds")
