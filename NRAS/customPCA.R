customPCA <- function(eset, phenotype = NULL, phenotype2 = NULL, label = NULL, z_score = TRUE, 
                       show_ellipse = TRUE, components = 1:2, biplot = FALSE, biplot_labels = NULL, 
                       standardize = TRUE, num_features = 6L, show_NA = TRUE, legend_title = phenotype, 
                       arrow_args = list(), label_args = list(), ...) 
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
    show_ellipse <- FALSE
    colorBy <- NULL
    alphaBy <- NULL
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
  p <- ggplot(data = df.u, mapping = aes(x = df.u[, 1], y = df.u[, 
                                                                 2], color = colorBy)) + geom_hline(yintercept = 0, lty = "longdash", 
                                                                                                    color = "darkgrey") + geom_vline(xintercept = 0, lty = "longdash", 
                                                                                                                                     color = "darkgrey") + labs(x = axis_labs[1], y = axis_labs[2]) + 
    theme_bw() + theme(aspect.ratio = 1)
  if (show_ellipse & !is.numeric(colorBy)) {
    p <- p + ggforce::geom_mark_ellipse(aes(fill = colorBy, alpha = alphaBy,
                                        color = colorBy)) + 
      scale_alpha_discrete(range=c(0.1,0.6))# + labs(alpha=alphaBy)
  }
  if (is.null(label)) {
    p <- p + geom_point(...)
  }
  else {
    if (!label %in% colnames(pData(eset))) {
      stop(sprintf("'%s' is not the name of a column in pData(eset).", 
                   label))
    }
    labels <- pData(eset)[, label]
    p <- p + geom_text(mapping = aes(label = labels), ...)
  }
  p <- p + guides(color = guide_legend(title = legend_title), 
                  fill = guide_legend(title = legend_title),
                  alpha = guide_legend(title=phenotype2))
  if (is.numeric(colorBy)) {
    p <- p + guides(color = guide_colorbar(title = legend_title))
  }
  if (biplot) {
    top_features <- lapply(1:2, function(i) {
      order(abs(df.v)[, i], decreasing = TRUE)[1:num_features]
    })
    top_features <- unique(unlist(top_features))
    df.v <- df.v[top_features, ]
    colnames(df.v) <- c("xend", "yend")
    df.v$x <- df.v$y <- 0
    if (is.null(biplot_labels)) {
      df.v$labels <- rownames(df.v)
    }
    else {
      df.v$labels <- fData(eset)[top_features, biplot_labels]
    }
    scale_args <- list(expand = expansion(mult = rep(0.1, 
                                                     2)), sec.axis = sec_axis(~. * ratio))
    arrow_args <- list(mapping = aes(x = x, y = y, xend = xend, 
                                     yend = yend), arrow = arrow(length = unit(0.5, "line")), 
                       data = df.v, color = "red3") %>% modifyList(val = arrow_args, 
                                                                   keep.null = TRUE)
    label_args <- list(mapping = aes(x = xend, y = yend, 
                                     label = labels), data = df.v, color = arrow_args[["color"]], 
                       max.overlaps = Inf, min.segment.length = 0, fill = alpha("white", 
                                                                                0.5)) %>% modifyList(val = label_args, keep.null = TRUE)
    p <- p + do.call(scale_x_continuous, scale_args) + do.call(scale_y_continuous, 
                                                               scale_args) + do.call(geom_segment, arrow_args) + 
      do.call(geom_label_repel, label_args) + theme(axis.text.y.right = element_text(color = arrow_args[["color"]]), 
                                                    axis.text.x.top = element_text(color = arrow_args[["color"]]), 
                                                    axis.ticks.y.right = element_line(color = arrow_args[["color"]]), 
                                                    axis.ticks.x.top = element_line(color = arrow_args[["color"]]))
  }
  return(p)
}
