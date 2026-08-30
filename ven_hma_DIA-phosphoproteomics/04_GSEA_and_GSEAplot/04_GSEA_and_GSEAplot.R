library(dplyr)
library(ggplot2)
library(openxlsx)
library(msigdbr)
library(fgsea)
library(scales)
library(synapser)

## synapse login with .Renviron
synLogin()

## define folder ID and the target files
folder_id <-"syn68733653"

file_id <- synFindEntityId("PTRC_EXP28_Phospho_Stats_Results.xlsx", parent = folder_id)
file_entity <- synGet(file_id)
df <- read.xlsx(file_entity$path)

## parse out 2 comparisons from Limma output file
NRRef <- df %>% filter(contrast == "Norelapse-Refractory")   
RRef  <- df %>% filter(contrast == "Relapse-Refractory")  

# suject NRRef, RRef to the following code for GSEA
# prepare input for GSEA
names(NRRef)[names(NRRef) == "logFC"] <- "log2FC"
names(NRRef)[names(NRRef) == "P.Value"] <- "pvalue"
NRRef$gene <- sub("-.*", "", NRRef$feature)

gene_level_NRRef <- NRRef[,c(1,4,9)] %>%
  filter(!is.na(log2FC), !is.na(pvalue), !is.na(gene)) %>%
  group_by(gene) %>%
  slice_max(order_by = abs(log2FC), n = 1, with_ties = FALSE) %>%
  ungroup()

# Build ranking statistic for GSEA: sign(log2FC) * -log10(p)
gene_level_NRRef <- gene_level_NRRef %>%
  mutate(rank_stat = sign(log2FC) * -log10(pvalue))

# Create named numeric vector, sorted decreasing
ranks_NRRef <- gene_level_NRRef$rank_stat
names(ranks_NRRef) <- gene_level_NRRef$gene
ranks_NRRef <- sort(ranks_NRRef, decreasing = TRUE)

# Hallmark pathways (H collection), human
m_df <- msigdbr(species = "Homo sapiens", collection = "H")

# Convert to a list: names = pathway, each element = vector of genes
pathways <- split(m_df$gene_symbol, m_df$gs_name)

set.seed(123)

fgsea_res_NRRef <- fgsea(
  pathways = pathways,
  stats    = ranks_NRRef,
  minSize  = 10,
  maxSize  = 500,
  nperm    = 10000
)

# calculate GeneRatio for plotting purpose
fgsea_res_NRRef <- fgsea_res_NRRef %>%
  mutate(
    GeneRatio = lengths(leadingEdge) / size,
    Count     = lengths(leadingEdge)       
  )

# filter out pathways with p-value <0.05 
pfgsea <- fgsea_res_NRRef %>% filter(pval <0.05) %>% arrange(desc(GeneRatio))
# clean up pathway labels
pfgsea$pathway <-sub("^HALLMARK_", "", pfgsea$pathway)
pfgsea$pathway <-gsub("_", " ", pfgsea$pathway)
pfgsea$pathway <- factor(pfgsea$pathway, levels=rev(pfgsea$pathway))

ggplot(pfgsea, aes(x = GeneRatio, y = pathway)) +
  geom_point(aes(size = -log10(pval), color = NES)) +
  scale_y_discrete(labels = label_wrap(20))+
  scale_color_gradient2(low="blue", mid="white", high="red", midpoint=0, name="NES") +
  scale_size_continuous(name = "-log10(pval)", range=c(6,10)) +
  labs(
    x = "GeneRatio",
    y = NULL,
    title = ""
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 22, face="bold", color="black"),
    plot.title  = element_text(hjust = 0.5),
    axis.title.x=element_text(size=30, face="bold"),
    axis.title.y=element_text(size=27, face="bold", color="black"),
    axis.text.x=element_text(size=27, face="bold", color="black"),
    legend.title=element_text(size=22, face="bold"),
    legend.text=element_text(size=20),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  )
