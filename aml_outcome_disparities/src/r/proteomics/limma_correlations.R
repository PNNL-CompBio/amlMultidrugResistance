library(Biobase)
library(EnhancedVolcano)
library(limma)

#' LIMMA-based proteomic comparison
#'
#' Uses LIMMA to compare proteomic expression between black and white patients.
#' @param data DataFrame with expression measurements.
#' @param meta DataFrame with patient meta-data, including race.
#' @param p_cutoff Float corresponding to p-value cutoff for significance.
#' Defaults to 1e-02.
#' @param fc_cutoff Float corresponding to LFC cutoff for significance. Defaults
#' to 1.
#' @param save_plot String corresponding to where file should be saved.
#' @return [list] with the following elements:
#' \itemize{
#'  \item results: LIMMA results object
#'  \item volcano: Volcano plot comparing black & white patients
#' }
run_limma <- function(
  data, 
  meta,
  p_cutoff = 1e-02,
  fc_cutoff = 1,
  save_plot = NULL
) {
  # Build design matrix
  design <- model.matrix(
    ~0 + Race + Study + Age,
    data = meta
  )
  
  # Adjust PacificIslander category for syntax
  colnames(design)[
    which(colnames(design) == "RacePacific Islander")
  ] = "RacePacificIslander"
  
  # Trim to patients with meta-data, format to expression matrix
  data <- data[row.names(design),]
  data <- t(as.matrix(data))
  exp_set <- ExpressionSet(data)
  cont.matrix <- makeContrasts(RaceBlack - RaceWhite, levels = design)
  
  # Fit LIMMA, get contrasts by race
  data_fit <- lmFit(exp_set, design)
  race_contrasts <- contrasts.fit(data_fit, cont.matrix)
  race_contrasts <- eBayes(race_contrasts)
  results <- topTable(race_contrasts, number = Inf, adjust.method = "BH")
  
  # Plots volcano
  volcano = EnhancedVolcano(
    results,
    lab = rownames(results),
    x = "logFC",
    y = "adj.P.Val",
    pCutoff = p_cutoff,
    FCcutoff = fc_cutoff,
    title = NULL,
    subtitle = NULL,
    legendLabels = NULL,
    caption = NULL,
    ylim = c(0, max(-log10(results[["adj.P.Val"]]), na.rm = TRUE) + 1)
  )
  
  if (!is.null(save_plot)) {
    ggsave(save_plot)
  }
  
  return(list(results = results, volcano = volcano))
}
