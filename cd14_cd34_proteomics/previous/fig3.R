# how does mono sig compare to ven res signature?
library(plyr);library(dplyr);library(ggplot2);library(ggvenn)
# DIA bead
setwd("~/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/")
mono <- read.csv("analysis/combined24-27/using_cellType-sortType-patient_factors/DIA_2batches_noOutliers_noMSC/Sort Type_Bead/CD14_Pos_vs_Neg/global/Differential_expression/Differential_expression_results.csv")
ven <- read.csv("analysis/combined24-27/using_cellType-sortType-patient_factors/DIA_2batches_noOutliers_noMSC/Sort Type_Bead/Ven_Sensitive_vs_Resistant/global/Differential_expression/Differential_expression_results.csv")
ven$Log2FC <- -ven$Log2FC # make res vs sens instead of sens vs res
mono$Contrast <- "CD14+ vs. CD34+"
ven$Contrast <- "Ven Res vs. Sens"
venn.list <- list("CD14+ vs. CD34+" = mono[mono$adj.P.Val<0.05,c("Contrast","Gene","Log2FC")]$Gene,
                  "Ven Res vs. Sens" = ven[ven$adj.P.Val<0.05,c("Contrast","Gene","Log2FC")]$Gene)
sigs <- merge(mono[mono$adj.P.Val<0.05,c("Contrast","Gene","Log2FC")], 
              ven[ven$adj.P.Val<0.05,c("Contrast","Gene","Log2FC")], 
              by="Gene", suffixes=c(": CD14+ vs. CD34+", ": Ven Res vs. Sens"))
ggvenn::ggvenn(venn.list, show_percentage=FALSE, set_name_size=5, text_size=5)
ggplot2::ggsave("mono_ven_signatures_vennDiagram.pdf",width=7, height=7)

gsea.cor <- cor.test(sigs$`Log2FC: CD14+ vs. CD34+`, sigs$`Log2FC: Ven Res vs. Sens`)
# Pearson's product-moment correlation
# 
# data:  sigs$`NES: CD14+ vs. CD34+` and sigs$`NES: Ven Res vs. Sens`
# t = 20.091, df = 18, p-value = 8.904e-14
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.9451077 0.9916057
# sample estimates:
#       cor 
# 0.9784229 

Pearson.est <- gsea.cor$estimate
Pearson.p <- gsea.cor$p.value
stats_pearson <- substitute(
  r == est * "," ~ ~"p" ~ "=" ~ p,
  list(
    est = format(as.numeric(Pearson.est), digits = 3),
    p = format(Pearson.p, digits = 3)
  )
)
maxVal <- max(abs(c(sigs$`Log2FC: Ven Res vs. Sens`, sigs$`Log2FC: CD14+ vs. CD34+`)))
ggplot(sigs, aes(x=`Log2FC: CD14+ vs. CD34+`, y=`Log2FC: Ven Res vs. Sens`)) + geom_point() + theme_minimal() + 
  scale_x_continuous(limits=c(-maxVal, maxVal)) + scale_y_continuous(limits=c(-maxVal, maxVal)) + 
  labs(#x="Lasry RNA-seq: LRCC25", y = "Sorted Proteomics: NCF2", 
    title=paste0("Differential Expression (n = ",nrow(sigs),")")) +
  geom_smooth(method="lm", se=FALSE, linetype="dashed") + 
  ggrepel::geom_label_repel(data=sigs, aes(label=Gene), size=2.5) + 
  theme(plot.title=element_text(hjust=0.5)) +
  ggplot2::geom_text(
    x = -Inf, y = 0.5, vjust = "inward", hjust = "inward",
    colour = "blue", parse = TRUE, 
    label = as.character(as.expression(stats_pearson)), size = 4.5
  ) 
ggsave("Mono_vs_VenRes_DEGs.pdf", width=4, height=4)


av <- read.csv("analysis/combined24-27using_cellType-sortType-patient_factors/DIA_2batches_noOutliers_noMSC/Sort Type_Bead/Aza.Ven_Sensitive_vs_Resistant/global/Differential_expression/Differential_expression_results.csv")
# not there?

# what about GSEA?
mono <- read.csv("analysis/combined24-27/using_cellType-sortType-patient_factors/DIA_2batches_noOutliers_noMSC/Sort Type_Bead/CD14_Pos_vs_Neg/global/GSEA/GSEA_Hallmark/GSEA_results.csv")
ven <- read.csv("analysis/combined24-27/using_cellType-sortType-patient_factors/DIA_2batches_noOutliers_noMSC/Sort Type_Bead/Ven_Sensitive_vs_Resistant/global/GSEA/GSEA_Hallmark/GSEA_results.csv")
ven$NES <- -ven$NES # make res vs sens instead of sens vs res
mono$Contrast <- "CD14+ vs. CD34+"
ven$Contrast <- "Ven Res vs. Sens"
venn.list <- list("CD14+ vs. CD34+" = mono[mono$FDR_q_value<0.25 & mono$p_value<0.05,c("Contrast","Feature_set","NES")]$Feature_set,
                  "Ven Res vs. Sens" = ven[ven$FDR_q_value<0.25 & ven$p_value<0.05,c("Contrast","Feature_set","NES")]$Feature_set)
sigs <- merge(mono[mono$FDR_q_value<0.25 & mono$p_value<0.05,c("Contrast","Feature_set","NES")], 
              ven[ven$FDR_q_value<0.25 & ven$p_value<0.05,c("Contrast","Feature_set","NES")], 
              by="Feature_set", suffixes=c(": CD14+ vs. CD34+", ": Ven Res vs. Sens")) # 20 overlapping
ggvenn::ggvenn(venn.list, show_percentage = FALSE, set_name_size = 5, text_size=12)
ggsave("GSEA_venn_mono_venRes.pdf",width=5,height=5)
ggvenn::ggvenn(venn.list, show_percentage = FALSE, set_name_size = 5, text_size=5)
ggsave("GSEA_venn_mono_venRes_size5.pdf",width=5,height=5)
gsea.cor <- cor.test(sigs$`NES: CD14+ vs. CD34+`, sigs$`NES: Ven Res vs. Sens`)
# Pearson's product-moment correlation
# 
# data:  sigs$`NES: CD14+ vs. CD34+` and sigs$`NES: Ven Res vs. Sens`
# t = 20.091, df = 18, p-value = 8.904e-14
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.9451077 0.9916057
# sample estimates:
#       cor 
# 0.9784229 

Pearson.est <- gsea.cor$estimate
Pearson.p <- gsea.cor$p.value
stats_pearson <- substitute(
  r == est * "," ~ ~"p" ~ "=" ~ p,
  list(
    est = format(as.numeric(Pearson.est), digits = 3),
    p = format(Pearson.p, digits = 3)
  )
)
maxVal <- max(abs(c(sigs$`NES: Ven Res vs. Sens`, sigs$`NES: CD14+ vs. CD34+`)))
sigs$Feature_set <- sub("HALLMARK_","", sigs$Feature_set)
sigs$inMain <- FALSE
sigs[sigs$Feature_set %in% c("COMPLEMENT","TNFA_SIGNALING_VIA_NFKB","COAGULATION",
                             "INTERFERON_GAMMA_RESPONSE","INFLAMMATORY_RESPONSE",
                             "ALLOGRAFT_REJECTION","IL6_JAK_STAT3_SIGNALING","APICAL_JUNCTION",
                             "OXIDATIVE_PHOSPHORYLATION",
                             "MYC_TARGETS_V2"),]$inMain <- TRUE
ggplot(sigs, aes(x=`NES: CD14+ vs. CD34+`, y=`NES: Ven Res vs. Sens`)) + geom_point() + theme_minimal() + 
  scale_x_continuous(limits=c(-maxVal, maxVal)) + scale_y_continuous(limits=c(-maxVal, maxVal)) + 
  labs(#x="Lasry RNA-seq: LRCC25", y = "Sorted Proteomics: NCF2", 
    title=paste0("Gene Set Enrichment (n = ",nrow(sigs),")")) +
  geom_smooth(method="lm", se=FALSE, linetype="dashed") + 
  ggrepel::geom_label_repel(data=subset(sigs, inMain), aes(label=Feature_set), size=2.5) + 
  theme(plot.title=element_text(hjust=0.5)) +
  ggplot2::geom_text(
    x = -Inf, y = 0.5, vjust = "inward", hjust = "inward",
    colour = "blue", parse = TRUE, 
    label = as.character(as.expression(stats_pearson)), size = 4.5
  ) 
ggsave("Mono_vs_VenRes_GSEA_NES.pdf", width=4, height=4)

#### let's also look at the top DEGs for Ven Res vs. Sens and how this agrees with previous literature ####
##### top Ven res vs. sens degs
n.top <- 5
sig.diffexp <- ven[ven$adj.P.Val < 0.05,]
top.pos.diffexp <- sig.diffexp %>% slice_max(Log2FC, n = n.top)
top.neg.diffexp <- sig.diffexp %>% slice_min(Log2FC, n = n.top)
top.diffexp <- rbind(top.pos.diffexp, top.neg.diffexp)
top.diffexp$Significance <- "Upregulated"
top.diffexp[top.diffexp$Log2FC < 0,]$Significance <- "Downregulated"
top.diffexp$Significance <- factor(top.diffexp$Significance, levels=c("Upregulated","Downregulated"))

# vertical bar plot of top 5 sig results
nSig <- nrow(sig.diffexp)
nTotal <- nrow(ven)
title <- paste0("Differentially expressed genes (",nSig,"/",nTotal, " with adjusted p < 0.05)")
#setOrder <- top.diffexp[order(top.diffexp$Log2FC*-log10(top.diffexp$adj.P.Val)),]$Gene
setOrder <- top.diffexp[order(top.diffexp$Log2FC),]$Gene

# bold CD14 and CD34
geneFace <- rep("plain", nrow(top.diffexp))
names(geneFace) <- setOrder
geneFace[c("CD14","CD34")] <- "bold"

title <- paste0("Differentially expressed genes\n(",nSig,"/",nTotal, " with adjusted p < 0.05)")
top.diffexp$`-Log(FDR)` <- -log(top.diffexp$adj.P.Val, base=10)
absMaxLog2FC <- max(abs(top.diffexp$Log2FC))
maxLogFDR <- ceiling(max(top.diffexp$`-Log(FDR`))
top.diffexp$Significant <- TRUE
diffexp.dot <- ggplot(top.diffexp, aes(x=Gene, y="Ven Res vs. Sens", color=Log2FC, size=`-Log(FDR)`)) + geom_point()+
  ggplot2::scale_x_discrete(limits = setOrder) +
  theme(axis.text.x=element_text(face=geneFace)) +
  scale_size(limits=c(1,maxLogFDR), range = c(0.5,7)) +
  scale_color_gradient2(low="blue",high="red", mid="grey", limits=c(-absMaxLog2FC, absMaxLog2FC)) +
  geom_point(data = subset(top.diffexp, Significant), col = "black", stroke = 1.5, shape = 21) +
  theme_classic(base_size = 12) + ggtitle(title) + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), 
        axis.text = element_text(size=16)) + coord_flip()
diffexp.dot
ggsave(paste0("Ven_diffexp_top_",n.top,"_sig_absLog2FC_dotPlot_",Sys.Date(),".pdf"), diffexp.dot, width=2.5, height=4)


##### venn diagram
ven.res.beatAML <- c("ITGAM","FCGR3A", "FCGR1A","HLA-DRB1") # FCGR3B and NCAM1 are not detected in our proteomics
ven.sens.beatAML <- "KIT"
ven.beatAML <- c(ven.res.beatAML, ven.sens.beatAML)
venn.list2 <- append(venn.list, list("Beat AML: Ven Res vs. Sens" = ven.beatAML))
ggvenn::ggvenn(venn.list2, show_percentage=FALSE, set_name_size=5, text_size=5)
ggplot2::ggsave("mono_ven_BeatAML_signatures_vennDiagram.pdf",width=7, height=7)

# do up and down regulated separately
mono <- read.csv("analysis/combined24-27/using_cellType-sortType-patient_factors/DIA_2batches_noOutliers_noMSC/Sort Type_Bead/CD14_Pos_vs_Neg/global/Differential_expression/Differential_expression_results.csv")
ven <- read.csv("analysis/combined24-27/using_cellType-sortType-patient_factors/DIA_2batches_noOutliers_noMSC/Sort Type_Bead/Ven_Sensitive_vs_Resistant/global/Differential_expression/Differential_expression_results.csv")
ven$Log2FC <- -ven$Log2FC # make res vs sens instead of sens vs res
mono$Contrast <- "CD14+ vs. CD34+"
ven$Contrast <- "Ven Res vs. Sens"

venn.list.res <- list("CD14+" = mono[mono$adj.P.Val<0.05 & mono$Log2FC>0,c("Contrast","Gene","Log2FC")]$Gene,
                  "Ven Res" = ven[ven$adj.P.Val<0.05 & ven$Log2FC>0,c("Contrast","Gene","Log2FC")]$Gene,
                  "Beat AML: Ven Res" = ven.res.beatAML)
ggvenn::ggvenn(venn.list.res, show_percentage=FALSE, set_name_size=5, text_size=5)
ggplot2::ggsave("mono_venRes_BeatAML_signatures_vennDiagram.pdf",width=7, height=7)

venn.list.sens <- list("CD34+" = mono[mono$adj.P.Val<0.05 & mono$Log2FC<0,c("Contrast","Gene","Log2FC")]$Gene,
                      "Ven Sens" = ven[ven$adj.P.Val<0.05 & ven$Log2FC<0,c("Contrast","Gene","Log2FC")]$Gene,
                      "Beat AML: Ven Sens" = ven.sens.beatAML)
ggvenn::ggvenn(venn.list.sens, show_percentage=FALSE, set_name_size=5, text_size=5)
ggplot2::ggsave("mono_venSens_BeatAML_signatures_vennDiagram.pdf",width=7, height=7)

# how many of the overlapping DEGs are in the same direction?
# all of the overlapping Beat AML markers are in the same direction as CD14+ vs. CD34+ signature
n.match.res <- venn.list.res[["CD14+"]][venn.list.res[["CD14+"]] %in% venn.list.res[["Ven Res"]]] # 18
n.match.sens <- venn.list.sens[["CD34+"]][venn.list.sens[["CD34+"]] %in% venn.list.sens[["Ven Sens"]]] # 142
n.overlap <- venn.list[["CD14+ vs. CD34+"]][venn.list[["CD14+ vs. CD34+"]] %in% venn.list[["Ven Res vs. Sens"]]] # 180
n.match <- length(n.match.res) + length(n.match.sens) # 160
perc.match <- 100*n.match/length(n.overlap) # 88.88...%
n.mismatch <- length(n.overlap) - n.match # 20
perc.mismatch <- 100*n.mismatch/length(n.overlap) # 11.11...%
