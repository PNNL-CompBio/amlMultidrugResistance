library(apeglm)
library(DESeq2)
library(EnhancedVolcano)
library(here)
library(synapser)

#' DESeq2 differential expression analysis
#'
#' Uses DESeq2 to compare transcriptomic expression between black and white
#' patients.
#' @param data DataFrame with unnormalized counts.
#' @param meta DataFrame with patient meta-data, including race.
#' @param p_cutoff Float corresponding to p-value cutoff for significance.
#' Defaults to 1e-02.
#' @param fc_cutoff Float corresponding to LFC cutoff for significance. Defaults
#' to 1.
#' @param save_plot String corresponding to where file should be saved.
#' @return [list] with the following elements:
#' \itemize{
#'  \item dds: DESeq2 results object
#'  \item volcano: Volcano plot comparing black & white patients
#' }
run_deseq2 <- function(
  data,
  meta,
  p_cutoff = 1e-02,
  fc_cutoff = 1,
  save_plot = NULL
) {
  # Trim to Black and White patients
  # More race comparisons would be interesting, but only Black and White 
  # patients have enough measurements for any useful interpretation
  meta <- meta[
    meta$Race %in% c("Black", "White"),
  ]
  shared_patients <- intersect(
    row.names(meta),
    row.names(data)
  )
  data <- data[shared_patients,]
  meta <- meta[shared_patients,]
  
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
  results <- results(dds, contrast=c("Race","Black","White"))
  
  # Run LFC shrinkage, trim to LFC > 0
  results <- lfcShrink(dds, coef="Race_White_vs_Black")
  
  # Plot volcano
  volcano <- EnhancedVolcano(
    results,
    lab = rownames(results),
    x = 'log2FoldChange',
    y = 'padj',
    pCutoff = p_cutoff,
    FCcutoff = fc_cutoff,
    title = NULL,
    subtitle = NULL,
    legendLabels = NULL,
    caption = NULL,
    ylim = c(0, max(-log10(results[["padj"]]), na.rm = TRUE) + 1)
  )
  
  if (!is.null(save_plot)) {
    ggsave(save_plot)
  }
  
  return(list(results = results, volcano = volcano))
}