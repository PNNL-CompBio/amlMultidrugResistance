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
