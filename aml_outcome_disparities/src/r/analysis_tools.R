library(EnhancedVolcano)
library(limma)
library(mixOmics)


#' Runs PLSR using mixOmics's NA-friendly version.
#'
#' This assumes data is z-scored prior to analysis, so scale measurements (both
#' data and labels) prior to function call.
#' 
#' @param data DataFrame. Data to regress.
#' @param labels DataFrame. Labels to regress data against.
#' @param n_comp Integer. Number of components.
#' @return [list] with the following elements:
#' \itemize{
#'  \item strip_plot: Strip plot of scores for each Race, mutation status
#'  \item result: LIMMA result object
#' }
run_plsr <- function(data, labels, n_comp = 2, new_data = NULL) {
  spls.result <- pls(data, labels, ncomp = n_comp, scale = FALSE)
  x_scores <- data.frame(spls.result$variates$X)
  
  if (!is.null(new_data)) {
    external <- predict(spls.result, new_data)
    colnames(external$variates) <- colnames(x_scores) 
    x_scores <- rbind(x_scores, external$variates)
  }
  
  return(
    list(
      x_scores = x_scores,
      y_scores = data.frame(spls.result$variates$Y),
      x_loadings = data.frame(spls.result$loadings$X),
      y_loadings = data.frame(spls.result$loadings$Y),
      var_explained = data.frame(spls.result$prop_expl_var)
    )
  )
}

#' Runs PLSR using mixOmics's NA-friendly version, returns deflated matrices.
#'
#' @param data DataFrame. Data to regress.
#' @param labels DataFrame. Labels to regress data against.
#' @param n_comp Integer. Number of components.
#' @return [list] with the following elements:
#' \itemize{
#'  \item strip_plot: Strip plot of scores for each Race, mutation status
#'  \item result: LIMMA result object
#' }
deflated_plsr <- function(data, labels, n_comp = 2, new_data = NULL) {
  if (!is.null(new_data)) {
    scaled <- rbind(data, new_data)
    scaled <- scale(scaled)
    data <- scaled[rownames(data),]
    new_data <- scaled[rownames(new_data),]
  } else {
    data <- scale(data)
  }
  labels <- scale(labels)
  
  spls.result <- pls(data, labels, ncomp = n_comp, scale = FALSE)
  x_scores <- data.frame(spls.result$variates$X)
  
  if (!is.null(new_data)) {
    external <- predict(spls.result, new_data)
    colnames(external$variates) <- colnames(x_scores) 
    x_scores <- rbind(x_scores, external$variates)
  }
  
  return(
    list(
      x_scores = x_scores,
      y_scores = data.frame(spls.result$variates$Y),
      x_loadings = data.frame(spls.result$loadings$X),
      y_loadings = data.frame(spls.result$loadings$Y),
      var_explained = data.frame(spls.result$prop_expl_var)
    )
  )
}

#' Runs LIMMA for correlations.
#'
#' @param data DataFrame. Data to regress.
#' @param values DataFrame. Values to regress against.
#' @return DataFrame with correlations between values and measurements in data.
limma_correlation <- function(data, values) {
  model_formula <- sprintf("~ %s", colnames(values)[1])
  design <- model.matrix(as.formula(model_formula), data = values)
  
  # Fit the linear model
  fit <- lmFit(t(data), design)
  fit <- eBayes(fit)
  results <- topTable(
    fit, 
    number = Inf, 
    adjust.method = "BH", 
    coef = colnames(values)[1]
  )
  
  return(results)
}

#' LIMMA-based proteomic comparison
#'
#' Uses LIMMA to compare expression across a mutation
#' @param data DataFrame with expression measurements.
#' @param meta DataFrame with patient meta-data, including race.
#' @param gene Character string corresponding to mutation gene to compare
#' @param p_cutoff Float corresponding to p-value cutoff for significance.
#' Defaults to 1e-02.
#' @return [list] with the following elements:
#' \itemize{
#'  \item results: LIMMA results object
#'  \item volcano: Volcano plot comparing black & white patients
#' }
limma_mutation <- function(
  data,
  meta,
  gene, 
  p_cutoff = 5e-02, 
  fc_cutoff = 1,
  save_plot = NULL
) {
  meta$Gene <- as.factor(meta[,gene])
  
  # Build design matrix
  design <- model.matrix(
    ~0 + Gene + Study,
    data = meta,
    contrasts.arg = list(Gene = contrasts(meta$Gene, contrasts = FALSE))
  )
  
  # Trim to patients with meta-data, format to expression matrix
  data <- data[row.names(design),]
  data <- t(as.matrix(data))
  
  # Fit LIMMA, get contrasts by race
  data_fit <- lmFit(data, design)
  
  if (!is.null(save_plot)) {
    cont.matrix <- makeContrasts(
      GeneMutant - GeneWT, 
      levels = design
    )
    mutation_contrasts <- contrasts.fit(data_fit, cont.matrix)
    mutation_contrasts <- eBayes(mutation_contrasts)
    results <- topTable(mutation_contrasts, number = Inf, adjust.method = "BH")
    
    volcano = EnhancedVolcano(
      results,
      lab = rownames(results),
      x = 'logFC',
      y = 'adj.P.Val',
      pCutoff = p_cutoff,
      FCcutoff = fc_cutoff,
      title = save_plot,
      subtitle = sprintf("%s Mutation", gene)
    )
    x_lim <- max(abs(results$logFC)) + 1
    volcano <- volcano + coord_cartesian(
      xlim = c(-x_lim, x_lim), 
      ylim = c(0, -log10(min(results$adj.P.Val)) + 0.5)
    )
    ggsave(sprintf("%s_%s.pdf", save_plot, gene))
  }
  
  return(list(results = results, volcano = volcano))
}
