# adapted by Belinda B. Garana, first 32 lines sourced from Peter van Galen, 230411
# Check expression of genes of interest in data from single-cell AML paper (van Galen, Hovestadt et al, Cell 2019)

# Load required libraries
library(tidyverse)
library(Seurat)
library(ggforce)
library(cowplot)
external.path <- "~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/data/externalSignatures"

#### 0. Get protein-coding genes ####
# # source: https://support.bioconductor.org/p/62441/#62451
# library(Biostrings) ## dna to rna
# library(biomaRt)
# ensembl = useMart("ensembl", dataset=
#                     "hsapiens_gene_ensembl")

# source: https://bioconductor.org/packages/release/bioc/vignettes/ensembldb/inst/doc/proteins.html
library(ensembldb)
BiocManager::install("EnsDb.Hsapiens.v86")
library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86
## Evaluate whether we have protein annotation available
hasProteinData(edb)
listTables(edb)

#### 1. van Galen et al ####
# Start with a clean slate
rm(list=ls())

# Frequently used function
cutf <- function(x, f=1, d="/") sapply(strsplit(x, d), function(i) paste(i[f], collapse=d))

# Load AML paper data. This first needs to be downloaded as described in README.md
setwd("~/OneDrive - PNNL/Documents/GitHub/reanalyze-vanGalen-aml2019")
aml <- readRDS("Seurat_AML.rds")

# Change MUTZ3 labels because the data for the different labels are similar
aml$orig.ident <- gsub("MUTZ3.*", "MUTZ3", aml$orig.ident)
aml$orig.ident <- factor(aml$orig.ident, levels = unique(aml$orig.ident))

# Generate a metadata tibble and add gene expression (choose one gene)
metadata <- as_tibble(aml@meta.data, rownames = "cell")

# Filter for AML cells at Dx. Also remove "healthy" cells from AML patients since their gene expression may be aberrant
metadata.filter <- metadata %>% filter(grepl("AML.*D0", orig.ident)) %>%
  mutate(Donor = ifelse(grepl("BM", orig.ident), yes = "Healthy", no = "AML")) %>%
  mutate(Donor = factor(Donor, levels = c("Healthy", "AML")))

#metadata.filter <- metadata.filter[metadata.filter$PredictionRefined == "malignant",]
# rest is inspired by: https://hbctraining.github.io/scRNA-seq_online/lessons/pseudobulk_DESeq2_scrnaseq.html
# prep data ----------------------------------------------------------------------------------
counts <- aml@assays$RNA$counts
#counts2 <- class(aml@assays$RNA$counts)

DefaultAssay(aml) # RNA
length(Cells(aml)) # 44823
length(levels(metadata.filter$CellType)) # 21
length(levels(metadata.filter$orig.ident)) # 43

# aggregate counts to the sample level for each cell type
# output: one aggregated counts matrix per cell type with genes along rows and samples along columns

groups <- dplyr::distinct(metadata.filter[,c("CellType","orig.ident")]) # 68 rows
agg.counts <- AggregateExpression(aml, group.by = c("CellType", "orig.ident")) 

data.table::tstrsplit(colnames(agg.counts$RNA), "_") %>% str()

cluster.counts <- list()
clusters <- levels(metadata.filter$CellType)
cluster.meta <- list()
for (i in clusters) {
  # filter for cluster
  colIDs <- colnames(agg.counts$RNA)[data.table::tstrsplit(colnames(agg.counts$RNA), "_")[[1]] == i]
  
  # also filter for AML D0 samples
  colIDs <- colIDs[grepl("AML.*D0", colIDs)]
  
  # store data for each cluster
  cluster.counts[[i]] <- agg.counts$RNA[,colIDs]
  cluster.meta[[i]] <- metadata.filter[metadata.filter$CellType == i,]
}
all(names(cluster.counts) == names(cluster.meta))
#all.counts <- data.table::rbindlist(cluster.counts, use.names = TRUE, idcol = "CellType")

# reformat for diffexp
factor.df <- as.data.frame(dplyr::distinct(metadata.filter[,c("CellType","orig.ident")]))
rownames(factor.df) <- paste0(factor.df$CellType,"_",factor.df$orig.ident)
#factor.df$orig.ident <- NULL

# mono-like vs. prog-like
cell.selection <- c("Mono-like", "Prog-like")
factor.df1 <- factor.df[factor.df$CellType %in% cell.selection,]
#factor.df1$CellType <- factor(factor.df1$CellType, levels = cell.selection)
# Error in limma::makeContrasts(contrasts = cts, levels = design) : 
#   The levels must by syntactically valid names in R, see help(make.names).  Non-valid names: Mono-like,Prog-like
factor.df1[factor.df1$CellType == "Mono-like",]$CellType <- "Mono"
factor.df1[factor.df1$CellType == "Prog-like",]$CellType <- "Prog"
factor.df1$CellType <- factor(factor.df1$CellType, levels = c("Mono", "Prog"))
factor.df1$orig.ident <- NULL

input.df <- as.data.frame(agg.counts$RNA)
keepCols <- colnames(input.df)
keepCols <- keepCols[grepl("Mono-like", keepCols) | grepl("Prog-like", keepCols)]
keepCols <- keepCols[grepl("AML.*D0", keepCols)]
keepCols <- keepCols[keepCols %in% rownames(factor.df1)]
input.df <- input.df[,keepCols]
input.df$Gene <- rownames(input.df)
input.list <- list("van_Galen_AML_D0" = input.df)
diffexp <- panSEA::mDEG(input.list, factor.df1)
write.csv(diffexp$all.results$van_Galen_AML_D0, file.path(external.path,"formatted/notFilteredForMalignant/Differential_expression_van_Galen_AML_D0_Mono-like_vs_Prog-like.csv"), row.names = FALSE)

diffexp.result <- read.csv(file.path(external.path,"formatted/notFilteredForMalignant/Differential_expression_van_Galen_AML_D0_Mono-like_vs_Prog-like.csv")) # 27899
write.csv(na.omit(diffexp.result), file.path(external.path,"formatted/notFilteredForMalignant/Differential_expression_van_Galen_AML_D0_Mono-like_vs_Prog-like_noNA.csv"), row.names = FALSE)
diffexp.result <- na.omit(diffexp.result) # 17075
# gb <- getBM(attributes=c("hgnc_symbol","gene_biotype"),filters = c("hgnc_symbol","biotype"), values=list(diffexp.result$Gene,"protein_coding"), mart=ensembl)
txs <- transcripts(edb, filter=GeneNameFilter(diffexp.result$Gene), columns = "tx_biotype")
protein.coding.genes <- txs[txs$tx_biotype == "protein_coding",]$gene_name
diffexp.result <- diffexp.result[diffexp.result$Gene %in% protein.coding.genes,] # 13828
diffexp.result$adj.P.Val <- 1
diffexp.result$adj.P.Val <- qvalue::qvalue(p = diffexp.result$P.Value, pi0 = 1)$qvalues
write.csv(diffexp.result, file.path(external.path,"formatted/notFilteredForMalignant/Differential_expression_van_Galen_AML_D0_Mono-like_vs_Prog-like_protein-coding.csv"), row.names = FALSE)
nrow(diffexp.result[diffexp.result$adj.P.Val <= 0.05,]) # 13

# mono-like vs. other
factor.df2 <- factor.df
#factor.df1$CellType <- factor(factor.df1$CellType, levels = cell.selection)
# Error in limma::makeContrasts(contrasts = cts, levels = design) : 
#   The levels must by syntactically valid names in R, see help(make.names).  Non-valid names: Mono-like,Prog-like
factor.df2$CellType <- sub("-like.*","",factor.df2$CellType)
factor.df2[factor.df2$CellType != "Mono",]$CellType <- "Other"
factor.df2$CellType <- factor(factor.df2$CellType, levels = c("Mono", "Other"))
factor.df2$orig.ident <- NULL

input.df <- as.data.frame(agg.counts$RNA)
keepCols <- colnames(input.df)
keepCols <- keepCols[grepl("AML.*D0", keepCols)]
keepCols <- keepCols[keepCols %in% rownames(factor.df2)]
input.df <- input.df[,keepCols]
input.df$Gene <- rownames(input.df)
input.list <- list("van_Galen_AML_D0" = input.df)
diffexp_other <- panSEA::mDEG(input.list, factor.df2)
write.csv(diffexp_other$all.results$van_Galen_AML_D0, file.path(external.path,"formatted/notFilteredForMalignant/Differential_expression_van_Galen_AML_D0_Mono-like_vs_Other.csv"), row.names = FALSE)
diffexp.result <- read.csv(file.path(external.path,"formatted/notFilteredForMalignant/Differential_expression_van_Galen_AML_D0_Mono-like_vs_Other.csv")) # 27899
#gb <- getBM(attributes=c("hgnc_symbol","gene_biotype"),filters = c("hgnc_symbol","biotype"), values=list(diffexp.result$Gene,"protein_coding"), mart=ensembl)
txs <- transcripts(edb, filter=GeneNameFilter(diffexp.result$Gene), columns = "tx_biotype")
protein.coding.genes <- txs[txs$tx_biotype == "protein_coding",]$gene_name
diffexp.result <- diffexp.result[diffexp.result$Gene %in% protein.coding.genes,] # 18972
diffexp.result$adj.P.Val <- 1
diffexp.result$adj.P.Val <- qvalue::qvalue(p = diffexp.result$P.Value, pi0 = 1)$qvalues
write.csv(diffexp.result, file.path(external.path,"formatted/notFilteredForMalignant/Differential_expression_van_Galen_AML_D0_Mono-like_vs_Other_protein-coding.csv"), row.names = FALSE)
nrow(diffexp.result[diffexp.result$adj.P.Val <= 0.05,]) # 4754

#### 2. Lasry et al ####
#lasry <- Matrix::readMM("~/Downloads/lasry/RNA_soupX_data.mtx.gz")
lmeta <- read.csv("~/Documents/misc/lasry/metadata_clustering_w_header_upd.csv")
lmeta <- lmeta[2:nrow(lmeta),]
lmeta <- lmeta[lmeta$ap_aml_age == "adult_AML" & lmeta$malignant=="malignant",]
#lcells <- read.csv("~/Downloads/lasry/cells_RNA_soupX_data.csv", header = FALSE)
#lgenes <- read.csv("~/Downloads/lasry/features_RNA_soupX_data.csv", header = FALSE)

# create seurat object
lasry <- Seurat::ReadMtx(mtx = "~/Documents/misc/lasry/RNA_soupX_data.mtx.gz", 
                         features = "~/Documents/misc/lasry/features_RNA_soupX_data.csv",
                         cells = "~/Documents/misc/lasry/cells_RNA_soupX_data.csv",
                         feature.column = 1)
rownames(lmeta) <- lmeta$NAME
lasry.seur <- CreateSeuratObject(counts = lasry, meta.data = lmeta)

# aggregate counts for cell type_patient level
lasry.df <- AggregateExpression(lasry.seur, group.by = c("Broad_cell_identity", "orig.ident"))

# extract data for relevant samples
factor.df <- as.data.frame(dplyr::distinct(lmeta[,c("Broad_cell_identity","orig.ident")]))
rownames(factor.df) <- paste0(factor.df$Broad_cell_identity,"_",factor.df$orig.ident)
#factor.df$orig.ident <- NULL

# mono-like vs. prog-like
#cell.selection <- c("CD14+ monocyte", "HSC")
cell.selection <- c("CD14+ monocyte", "MPP")
factor.df1 <- factor.df[factor.df$Broad_cell_identity %in% cell.selection,]
factor.df1[factor.df1$Broad_cell_identity == "CD14+ monocyte",]$Broad_cell_identity <- "Mono"
#factor.df1$Broad_cell_identity <- factor(factor.df1$Broad_cell_identity, levels = c("Mono", "HSC"))
factor.df1$Broad_cell_identity <- factor(factor.df1$Broad_cell_identity, levels = c("Mono", "MPP"))

input.df <- as.data.frame(lasry.df$RNA)
keepCols <- colnames(input.df)
keepCols <- keepCols[keepCols %in% rownames(factor.df1)]
factor.df1 <- factor.df1[colnames(input.df),]
factor.df1$orig.ident <- NULL
input.df <- input.df[,keepCols]
input.df$Gene <- rownames(input.df)
input.list <- list("Lasry" = input.df)
factor.df1 <- na.omit(factor.df1)
diffexp <- panSEA::mDEG(input.list, factor.df1)
write.csv(na.omit(diffexp$all.results$Lasry), file.path(external.path,"formatted/Differential_expression_Lasry_AML_CD14PosMonocyte_vs_MPP.csv"), row.names = FALSE)

#diffexp.result <- read.csv(file.path(external.path,"formatted/Differential_expression_Lasry_AML_CD14PosMonocyte_vs_HS.csv"))
diffexp.result <- na.omit(diffexp$all.results$Lasry) # 22744
txs <- transcripts(edb, filter=GeneNameFilter(diffexp.result$Gene), columns = "tx_biotype")
protein.coding.genes <- txs[txs$tx_biotype == "protein_coding",]$gene_name
diffexp.result <- diffexp.result[diffexp.result$Gene %in% protein.coding.genes,] # 15078
write.csv(diffexp.result, file.path(external.path,"formatted/Differential_expression_Lasry_AML_CD14PosMonocyte_vs_MPP_protein-coding.csv"), row.names = FALSE)
nrow(diffexp.result[diffexp.result$adj.P.Val <= 0.05,]) # 910

# # mono-like vs. other
# cell.selection <- c("CD14+ monocyte", "HSC", "MPP")
# factor.df2 <- factor.df[factor.df$Broad_cell_identity %in% cell.selection,]
# factor.df2[factor.df2$Broad_cell_identity == "CD14+ monocyte",]$Broad_cell_identity <- "Mono"
# factor.df2[factor.df2$Broad_cell_identity != "CD14+ monocyte",]$Broad_cell_identity <- "Other"
# factor.df2$Broad_cell_identity <- factor(factor.df2$Broad_cell_identity, levels = c("Mono", "Other"))
# 
# input.df <- as.data.frame(lasry.df$RNA)
# keepCols <- colnames(input.df)
# keepCols <- keepCols[keepCols %in% rownames(factor.df2)]
# factor.df2 <- factor.df2[rownames(factor.df2)[rownames(factor.df2) %in% colnames(input.df)],]
# factor.df2$orig.ident <- NULL
# input.df <- input.df[,keepCols]
# input.df$Gene <- rownames(input.df)
# input.list <- list("Lasry" = input.df)
# diffexp <- panSEA::mDEG(input.list, factor.df2)
# #write.csv(na.omit(diffexp$all.results$Lasry), file.path(external.path,"formatted/Differential_expression_Lasry_AML_CD14PosMonocyte_vs_HSCandMPP.csv"), row.names = FALSE)
# # Coefficients not estimable: Mono 
# # Error in limma::contrasts.fit(fit, cont.matrix) : 
# #   trying to take contrast of non-estimable coefficient
# # were there only mono-like instead of monocyte in van Galen because I filtered for malignant?

#### 2.1 Lasry et al - not filtered for malignant ####
#lasry <- Matrix::readMM("~/Downloads/lasry/RNA_soupX_data.mtx.gz")
lmeta <- read.csv("~/Documents/misc/lasry/metadata_clustering_w_header_upd.csv")
lmeta <- lmeta[2:nrow(lmeta),]
lmeta <- lmeta[lmeta$ap_aml_age == "adult_AML",]
#lcells <- read.csv("~/Downloads/lasry/cells_RNA_soupX_data.csv", header = FALSE)
#lgenes <- read.csv("~/Downloads/lasry/features_RNA_soupX_data.csv", header = FALSE)

# create seurat object
lasry <- Seurat::ReadMtx(mtx = "~/Documents/misc/lasry/RNA_soupX_data.mtx.gz", 
                         features = "~/Documents/misc/lasry/features_RNA_soupX_data.csv",
                         cells = "~/Documents/misc/lasry/cells_RNA_soupX_data.csv",
                         feature.column = 1)
rownames(lmeta) <- lmeta$NAME
lasry.seur <- CreateSeuratObject(counts = lasry, meta.data = lmeta)

# aggregate counts for cell type_patient level
lasry.df <- AggregateExpression(lasry.seur, group.by = c("Broad_cell_identity", "orig.ident"))

# extract data for relevant samples
factor.df <- as.data.frame(dplyr::distinct(lmeta[,c("Broad_cell_identity","orig.ident")]))
rownames(factor.df) <- paste0(factor.df$Broad_cell_identity,"_",factor.df$orig.ident)
#factor.df$orig.ident <- NULL
unique(factor.df$Broad_cell_identity)
# [1] "HSC"               "MPP"               "Ery"               "CD14+ monocyte"    "HLA-II+ monocyte"  "CD4+ T"            "Granulocyte"       "CD8+ T"           
# [9] "GMP"               "CD11c+"            "NK"                "MEP"               "B"                 "DC precursor"      "CD16+ monocyte"    "cDC2"             
# [17] "MAIT"              "LymP"              "cDC1"              "gd T"              "Plasmablast"       "Plasma cell"       "pDC"               "Megakaryocyte"    
# [25] "Pre-B"             "Perivascular cell" "Pro-B" 

# mono-like vs. prog-like
#cell.selection <- c("CD14+ monocyte", "HSC")
cell.selection <- c("CD14+ monocyte", "MPP")
factor.df1 <- factor.df[factor.df$Broad_cell_identity %in% cell.selection,]
factor.df1[factor.df1$Broad_cell_identity == "CD14+ monocyte",]$Broad_cell_identity <- "Mono"
#factor.df1$Broad_cell_identity <- factor(factor.df1$Broad_cell_identity, levels = c("Mono", "HSC"))
factor.df1$Broad_cell_identity <- factor(factor.df1$Broad_cell_identity, levels = c("Mono", "MPP"))

input.df <- as.data.frame(lasry.df$RNA)
keepCols <- colnames(input.df)
keepCols <- keepCols[keepCols %in% rownames(factor.df1)]
factor.df1 <- factor.df1[colnames(input.df),]
factor.df1$orig.ident <- NULL
factor.df1 <- na.omit(factor.df1)
input.df <- input.df[,keepCols]
input.df$Gene <- rownames(input.df)
input.list <- list("Lasry" = input.df)
diffexp <- panSEA::mDEG(input.list, factor.df1)
write.csv(na.omit(diffexp$all.results$Lasry), file.path(external.path,"formatted/notFilteredForMalignant/Differential_expression_Lasry_AML_CD14PosMonocyte_vs_MPP.csv"), row.names = FALSE)

#diffexp.result <- read.csv(file.path(external.path,"formatted/Differential_expression_Lasry_AML_CD14PosMonocyte_vs_HS.csv")) # 27899
diffexp.result <- na.omit(diffexp$all.results$Lasry) # 22896
txs <- transcripts(edb, filter=GeneNameFilter(diffexp.result$Gene), columns = "tx_biotype")
protein.coding.genes <- txs[txs$tx_biotype == "protein_coding",]$gene_name
diffexp.result <- diffexp.result[diffexp.result$Gene %in% protein.coding.genes,] # 15136
write.csv(diffexp.result, file.path(external.path,"formatted/notFilteredForMalignant/Differential_expression_Lasry_AML_CD14PosMonocyte_vs_MPP_protein-coding.csv"), row.names = FALSE)
nrow(diffexp.result[diffexp.result$adj.P.Val <= 0.05,]) # 922

#### make sure Triana genes are protein-coding ####
diffexp.result <- read.csv(file.path(external.path,"formatted/Triana_RNA_AML_100PercentCells_Classical-Monocytes_vs_HSCs-and-MPPs_differentialExpression.csv")) # 355, no NAs
txs <- transcripts(edb, filter=GeneNameFilter(diffexp.result$Gene), columns = "tx_biotype")
protein.coding.genes <- txs[txs$tx_biotype == "protein_coding",]$gene_name
diffexp.result <- diffexp.result[diffexp.result$Gene %in% protein.coding.genes,] # 339
write.csv(diffexp.result, file.path(external.path,"formatted/Triana_RNA_AML_100PercentCells_Classical-Monocytes_vs_HSCs-and-MPPs_differentialExpression_protein-coding.csv"), row.names = FALSE)
