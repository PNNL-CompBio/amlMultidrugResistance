PCA3D <- function(eset, phenotype = NULL, phenotype2 = NULL, z_score = TRUE, 
                       show_ellipse = TRUE, components = 1:2, standardize = TRUE, 
                  num_features = 6L, show_NA = TRUE, legend_title = phenotype, 
                      ...) 
{
  if (!is.null(phenotype)) {
    if (!phenotype %in% colnames(pData(eset))) {
      stop(sprintf("'%s' is not the name of a column in pData(eset).", 
                   phenotype))
    }
    colorBy <- pData(eset)[, phenotype]
    if (!is.null(phenotype2)) {
      if (!phenotype %in% colnames(pData(eset))) {
      stop(sprintf("'%s' is not the name of a column in pData(eset).", 
                   phenotype))
      }
    }
    alphaBy <- pData(eset)[, phenotype2]
    if (!show_NA) {
      idx <- !is.na(colorBy)
      eset <- eset[, idx]
      colorBy <- colorBy[idx]
      alphaBy <- alphaBy[idx]
    }
  }
  else {
    stop("phenotype and phenotype2 are required")
  }
  if (length(components) != 2) {
    stop(sprintf("components must be a vector of length 2, not %d.", 
                 length(components)))
  }
  if (!all(components %in% 1:ncol(eset))) {
    stop(sprintf("The values of components must be between 1 and %d.", 
                 ncol(eset)))
  }
  complete_rows <- complete.cases(exprs(eset))
  if (sum(complete_rows) < 2) {
    stop("There are fewer than 2 rows with non-missing data.")
  }
  message(sprintf("Subsetting to %d complete rows for PCA.", 
                  sum(complete_rows)))
  eset <- eset[complete_rows, ]
  if (z_score) {
    z <- t(scale(exprs(eset), center = TRUE, scale = TRUE))
  }
  else {
    z <- t(exprs(eset))
  }
  pca_res <- prcomp(z)
  u <- pca_res$x
  v <- pca_res$rotation
  if (standardize) {
    n <- nrow(u)
    lam <- pca_res$sdev * sqrt(n)
    u <- t(t(u)/lam)
    v <- t(t(v) * lam)
  }
  u_range <- apply(u[, components], 2, function(x) abs(range(x)))
  v_range <- apply(v[, components], 2, function(x) abs(range(x)))
  ratio <- max(v_range/u_range)
  v <- v/ratio
  df.u <- as.data.frame(u[, components])
  df.v <- as.data.frame(v[, components])
  d <- pca_res$sdev
  var_expl <- round(100 * d^2/sum(d^2), digits = 2)[components]
  axis_labs <- sprintf("PC%d (%g%%)", components, var_expl)
  if (!is.null(colorBy)) {
    df.u$colorBy <- colorBy
  }
  p <- plot_ly(df.u, x= ~df.u[, 1], y = ~df.u[, 2], z = ~alphaBy, color = ~colorBy)
  p <- p %>% add_markers()
  p <- p %>% layout(scene=list(xaxis=list(title=axis_labs[1]), yaxis=list(title=axis_labs[2]),
                               zaxis=list(title=phenotype2)))
  return(p)
}
