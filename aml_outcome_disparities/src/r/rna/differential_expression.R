library(apeglm)
library(DESeq2)
library(EnhancedVolcano)
library(here)
library(synapser)

here::i_am("src/r/rna/differential_expression.R")

# Import data
data <- read.csv(
  here("src", "r", "rna", "uncorrected.csv"),
  row.names = 1
)
meta <- read.csv(
  here("src", "r", "rna", "meta.csv"),
  row.names = 1
)

# Trim to Black and White patients
# More race comparisons would be interesting, but only Black and White patients
# have enough measurements for any useful interpretation
meta <- meta[
  meta$Race %in% c("Black", "White"),
]
data <- data[
  row.names(meta),
]

# Setup DESeq2 object and run
dds <- DESeqDataSetFromMatrix(
  t(data),
  colData = meta,
  design = ~ Source + Race
)

# Establish race as a factor, limit to Black & White
dds$Race <- factor(dds$Race, levels = c("Black", "White"))

# Run DESeq, encode Black v. White contrast
dds <- DESeq(dds)
res <- results(dds, contrast=c("Race","Black","White"))

# Run LFC shrinkage, trim to LFC > 0
res <- lfcShrink(dds, coef="Race_White_vs_Black")
res <- res[abs(res$log2FoldChange) > 1e-4,]

EnhancedVolcano(
  res,
  lab = rownames(res),
  x = 'log2FoldChange',
  y = 'padj',
  pCutoff = 1e-02
)
