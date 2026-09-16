library(MSnSet.utils)
library(dplyr)
library(KSEAapp)
library(ggplot2)
library(scales)
library(ggpubr)


## Modification of MSnSetUtils function.
plot_pca <- function(eset, phenotype = NULL, shape = NULL, label = NULL, z_score = TRUE,
                     princomp_center = TRUE, show_ellipse = TRUE, components = 1:2, biplot = FALSE,
                     biplot_labels = NULL, standardize = TRUE, save_dfs = NULL,
                     num_features = 6L, show_NA = TRUE, label_size = 3, output_type = NA,
                     legend_title = phenotype,
                     arrow_args = list(), label_args = list(), ...) {
   
   # Handling coloring by phenotype. Do this first, in case
   # rows are removed when show_NA = FALSE
   if (!is.null(phenotype)) {
      colorBy <- pData(eset)[, phenotype]
      # If not showing missing values, remove those samples
      if (!show_NA) {
         idx <- !is.na(colorBy)
         eset <- eset[, idx]
         colorBy <- colorBy[idx]
      }
   } else {
      show_ellipse <- FALSE
      colorBy <- NULL
   }
   if (!is.null(shape)){
      shapeBy <- pData(eset)[, shape]
      if (!show_NA) {
         idx <- !is.na(shapeBy)
         eset <- eset[, idx]
         shapeBy <- shapeBy[idx]
      }
   } else {
      shapeBy <- NULL
   }
   
   # Check that components are valid
   if (length(components) != 2) {
      stop(sprintf("components must be a vector of length 2, not %d.",
                   length(components)))
   }
   if (!all(components %in% 1:ncol(eset))) {
      stop(sprintf("The values of components must be between 1 and %d.",
                   ncol(eset)))
   }
   
   complete_rows <- complete.cases(exprs(eset))
   
   # Check that there are enough complete rows for PCA
   if (sum(complete_rows) < 2) {
      stop("There are fewer than 2 rows with non-missing data.")
   }
   
   message(sprintf("Subsetting to %d complete rows for PCA.",
                   sum(complete_rows)))
   
   # Subset to complete rows
   eset <- eset[complete_rows, ]
   
   # If z_score, convert to Z-Scores by sample (row when transposed)
   if (z_score) {
      z <- t(scale(exprs(eset), center = TRUE, scale = TRUE))
   } else {
      z <- t(exprs(eset))
   }
   
   ## PCA
   # By default, center = TRUE, scale. = FALSE
   pca_res <- prcomp(z, center = princomp_center)
   
   u <- pca_res$x # Scores
   v <- pca_res$rotation # Eigenvectors
   
   if (standardize) {
      n <- nrow(u)
      lam <- pca_res$sdev * sqrt(n)
      
      # Scale u down and v up. Product is still the same
      u <- t(t(u) / lam)
      v <- t(t(v) * lam)
   }
   
   # Determine ratio between scale of v and u
   u_range <- apply(u[, components], 2, function(x) abs(range(x)))
   v_range <- apply(v[, components], 2, function(x) abs(range(x)))
   
   ratio <- max(v_range / u_range) # ratio for scaling v and secondary axes
   v <- v / ratio # scale v
   
   if (!is.null(save_dfs)){
      u <- u %>% as.data.frame()
      # If colorBy is not NULL, add that column to df
      if (!is.null(colorBy)) {
         u$colorBy <- colorBy
      }
      if (!is.null(shapeBy)) {
         u$shapeBy <- shapeBy
      }
      
      assign(save_dfs, list("sample_decomposition" = u, "feature_decomposition" = v %>% as.data.frame() %>%
                               mutate(feature = rownames(.))), envir = globalenv())
   }
   
   # Data frames for plotting
   df.u <- as.data.frame(u[, components])
   df.v <- as.data.frame(v[, components])
   
   # Percent of variance explained by each PC
   d <- pca_res$sdev # Standard deviations
   var_expl <- round(100 * d ^ 2 / sum(d ^ 2), digits = 2)[components]
   axis_labs <- sprintf("PC%d (%g%%)", #"%sPC%d (%g%%)",
                        # ifelse(obs.scale == 0, "Standardized ", ""),
                        components,
                        var_expl)
   
   # If colorBy is not NULL, add that column to df
   if (!is.null(colorBy)) {
      df.u$colorBy <- colorBy
   }
   if (!is.null(shapeBy)) {
      df.u$shapeBy <- shapeBy
   }
   
   ## Visualization
   # Base plot
   p <- ggplot(data = df.u, mapping = aes(x = df.u[, 1], y = df.u[, 2], color = colorBy, shape = shapeBy)) +
      geom_hline(yintercept = 0, lty = "longdash", color = "darkgrey") +
      geom_vline(xintercept = 0, lty = "longdash", color = "darkgrey") +
      labs(x = axis_labs[1], y = axis_labs[2]) +
      theme_bw() +
      theme(aspect.ratio = 1)
   
   # 50% confidence ellipse layer first so they are
   # beneath the layer of points or labels.
   if (show_ellipse & !is.numeric(colorBy)) {
      p <- p +
         stat_ellipse(mapping = aes(fill = colorBy, color = NULL),
                      geom = "polygon", type = "norm",
                      level = 0.5, alpha = 0.1, show.legend = TRUE)
   }
   
   # If label is NULL, add points. Otherwise, add labels
   if (is.null(label)) {
      p <- p +
         geom_point(...)
   } else {
      labels <- pData(eset)[, label]
      p <- p + geom_point(...) + 
         ggrepel::geom_label_repel(mapping = aes(label = labels), 
                                   size = label_size, ...)
   }
   
   # Set titles for color and fill legend
   p <- p +
      guides(color = guide_legend(title = legend_title),
             fill = guide_legend(title = legend_title))
   
   # If colorBy is numeric, use a colorbar
   if (is.numeric(colorBy)) {
      p <- p +
         guides(color = guide_colorbar(title = legend_title))
   }
   
   ## Biplot
   if (biplot) {
      # Get the indices of the top influential features
      # from each principal component. num_features determines how
      # many to select from each component.
      top_features <- lapply(1:2, function(i) {
         order(abs(df.v)[, i], decreasing = TRUE)[1:num_features]
      })
      top_features <- unique(unlist(top_features))
      
      # Subset loadings to top features and rename columns
      df.v <- df.v[top_features, ]
      colnames(df.v) <- c("xend", "yend")
      df.v$x <- df.v$y <- 0
      
      # If biplot_labels is not provided, default to row names
      if (is.null(biplot_labels)) {
         df.v$labels <- rownames(df.v)
      } else {
         df.v$labels <- fData(eset)[top_features, biplot_labels]
      }
      
      scale_args <- list(expand = expansion(mult = rep(0.1, 2)),
                         sec.axis = sec_axis(~ . * ratio))
      
      # Arguments for geom_segment
      arrow_args <- list(mapping = aes(x = x, y = y, xend = xend, yend = yend),
                         arrow = arrow(length = unit(0.5, "line")),
                         data = df.v, color = "red3") %>%
         # Allow user-supplied args to overwrite defaults
         modifyList(val = arrow_args, keep.null = TRUE)
      
      # Arguments for geom_label_repel
      label_args <- list(mapping = aes(x = xend, y = yend, label = labels),
                         data = df.v,
                         color = arrow_args[["color"]],
                         max.overlaps = Inf,
                         min.segment.length = 0,
                         fill = alpha("white", 0.5)) %>%
         # Allow user-supplied args to overwrite defaults
         modifyList(val = label_args, keep.null = TRUE)
      
      # Add segments with arrows and text labels
      p <- p +
         # Add extra padding around plot area and secondary axes for v units
         do.call(scale_x_continuous, scale_args) +
         do.call(scale_y_continuous, scale_args) +
         do.call(geom_segment, arrow_args) +
         do.call(geom_label_repel, label_args) +
         theme(axis.text.y.right = element_text(color = arrow_args[["color"]]),
               axis.text.x.top = element_text(color = arrow_args[["color"]]),
               axis.ticks.y.right = element_line(color = arrow_args[["color"]]),
               axis.ticks.x.top = element_line(color = arrow_args[["color"]]))
   }
   
   if (output_type == "full"){
      return(list("plot" = p, "data.u" = as.data.frame(u[, ]), "data.v" = as.data.frame(v[, ]),
                  "axis_labs" = axis_labs))
   } else{
      return(p)
   }
}

diffexp_helper <- function(m, contrast_var, contrasts){
  pData(m)$bgd_ <- pData(m)[[contrast_var]]
  pData(m)$Sample <- sampleNames(m)
  
  all_results <- data.frame()
  
  for (contrast in contrasts){
    contrast_groups = strsplit(contrast, "-")[[1]]
    contrast = paste0("bgd_", contrast_groups[[1]], "-bgd_", contrast_groups[[2]])
    limma_res <- limma_contrasts(m, model.str = "~0 + bgd_", 
                                 coef.str = "bgd_", contrasts = contrast) %>% as.data.frame()
    rownames(limma_res) <- limma_res$feature
    counter = 1
    
    m_contrast <- m[, m$bgd_ %in% contrast_groups]
    p_values_t_test <- vector(mode="character", length = nrow(limma_res))
    p_values_welch_test <- vector(mode="character", length = nrow(limma_res))
    p_values_wilcox_test <- vector(mode="character", length = nrow(limma_res))
    for (feature in limma_res$feature){
      data_df <- data.frame(value = exprs(m_contrast)[feature, ],
                            Sample = colnames(exprs(m_contrast))) %>%
        filter(!is.na(value)) %>%
        merge(pData(m_contrast) %>% select(Sample, bgd_), by = "Sample")
      
      p_values_t_test[[counter]] <- tryCatch({t.test(value ~ bgd_, data = data_df, 
                                                     alternative = "two.sided", var.equal = TRUE)[[3]]}, 
                                             error = function(e) {NA}) 
      p_values_welch_test[[counter]] <- tryCatch({t.test(value ~ bgd_, data = data_df, 
                                                         alternative = "two.sided", var.equal = FALSE)[[3]]}, 
                                                 error = function(e) {NA}) 
      p_values_wilcox_test[[counter]] <- tryCatch({wilcox.test(value ~ bgd_, data = data_df, 
                                                         alternative = "two.sided")[[3]]}, 
                                                 error = function(e) {NA}) 
      counter = counter + 1
    }
    
    limma_res <- limma_res %>%
      mutate(t_test_pval = as.numeric(p_values_t_test),
             t_test_adj = p.adjust(t_test_pval, method = "BH"),
             welch_pval = as.numeric(p_values_welch_test),
             welch_adj = p.adjust(welch_pval, method = "BH"),
             wilcox_pval = as.numeric(p_values_wilcox_test),
             wilcox_adj = p.adjust(wilcox_pval, method = "BH"))
    all_results <- rbind(all_results, limma_res)
  }
  
  return(all_results)
}


plot_features <- function (m, features, feature_name_col = NULL, color_by = NULL, mode = "box", 
          order_by = color_by){
  p_data <- pData(m) %>% select(-matches("sample name"))
  p_data[["sample name"]] <- rownames(p_data)
  p_data <- p_data[order(p_data[[order_by]]), ]
  p_data[['sample name']] <- factor(p_data[['sample name']], levels = p_data[['sample name']])
  if (is.null(feature_name_col)) 
    feature_names <- featureNames(m)
  else feature_names <- fData(m)[[feature_name_col]]
  
  idx <- c()
  # all(grepl("^.*-[STY]{1}[0-9]+.*$", feature_names)) finds if all features in msnset are phosphosites
  # !any(grepl("^.*-[STY]{1}[0-9]+.*$", features)) finds if NONE of the chosen features provided are phosphosites
  # when both true, we plot the phosphosites from feature_names belonging to the chosen features (proteins).
  if (all(grepl("^.*-[STY]{1}[0-9]+.*$", feature_names)) & !any(grepl("^.*-[STY]{1}[0-9]+.*$", features))){
    print("detected phospho data, collecting phosphosites of the supplied features")
    for (feature in features){
      idx <- c(idx, which(startsWith(feature_names, feature)))
    }
  } else{
    idx <- which(featureNames(m) %in% features)
  }
  
  x <- exprs(m[idx, ]) %>% as.data.frame() %>% tibble::rownames_to_column("feature") %>%
    tidyr::pivot_longer(-feature, names_to = "sample name", values_to = "abundance") %>% 
    inner_join(p_data, by = "sample name") %>%
    mutate(`sample name` = factor(`sample name`, levels = p_data[['sample name']]))
  if (length(features) > 1){
    text_element <- element_blank()
  } else {
    text_element <- element_text(angle = 45, 
                                 hjust = 1)
  }
  if (mode == "point"){
    p <- x %>% ggplot() + aes(x = `sample name`, y = abundance) + 
      geom_point(size = 3) + theme_bw() + theme(axis.text.x = text_element) + 
      facet_wrap(~ feature, scales = "free") 
    if (!is.null(color_by)) 
      p <- p + aes_string(color = color_by)
    p
  } else if (mode == "box") {
    p <- x %>% ggplot() + aes(x = color_by, y = abundance) + 
      geom_boxplot(notch = TRUE) + theme_bw() + theme(axis.text.x = text_element) + 
      facet_wrap(~ feature, scales = "free") 
    if (!is.null(color_by)) 
      p <- p + aes_string(color = color_by)
    p
  }

}


volcano_function <- function(diffexp_res, chosen_terms, t2g, label_genes = c("Tkfc")){
  sig_genes <- diffexp_res %>% filter(t_test_adj < 0.05) %>% pull(feature) %>% unique()
  t2g_annotation <- t2g %>%
    filter(gene_symbol %in% sig_genes,
           gs_name %in% chosen_terms) %>%
    group_by(gs_name) %>%
    mutate(total = n()) %>%
    ungroup() %>%
    group_by(gene_symbol) %>%
    top_n(n = 1, wt = -total) %>%
    slice(1) %>%
    select(feature = gene_symbol, pathway = gs_name)
  
  ## Remove later
  t2g_annotation$pathway <- gsub("_", " ", t2g_annotation$pathway)
  t2g_annotation$pathway <- gsub("HALLMARK ", "", t2g_annotation$pathway)
  ###
  
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  
  diffexp_res_annotated <- left_join(diffexp_res, t2g_annotation, by = "feature") %>%
    mutate(pathway = case_when(!is.na(pathway) ~ pathway,
                               TRUE ~ 'NA')) %>%
    mutate(pathway = factor(pathway, levels = c(t2g_annotation$pathway %>% unique(), 'NA')),
           size = case_when(pathway != 'NA' ~ 4,
                            TRUE ~ 1.3)) %>%
    mutate(feat_lab = case_when(pathway != 'NA' ~ feature,
                                TRUE ~ '')) %>%
    arrange(t_test_adj)
  diffexp_res_annotated$rank <- 1:nrow(diffexp_res_annotated)
  
  diffexp_res_annotated$alpha[diffexp_res_annotated$feat_lab == ''] <- 0.85
  diffexp_res_annotated$alpha[diffexp_res_annotated$feat_lab != ''] <- 1
  diffexp_res_annotated$feat_lab[diffexp_res_annotated$feat_lab == ''] <- NA
  diffexp_res_annotated$feat_lab[!(diffexp_res_annotated$feature %in% label_genes)] <- NA
  
  
  p <- ggplot(diffexp_res_annotated, aes(x = logFC, y = -log10(t_test_adj), 
                                         color = pathway, alpha = alpha, size = size)) + geom_point() #+
    # scale_color_manual(values = c(gg_color_hue(length(unique(t2g_annotation$pathway))), "grey30"), 
    #                    name = "Pathway", limits = unique(t2g_annotation$pathway)) #+ 
    ## Add again later
    # geom_label_repel(aes(label = feat_lab), size = 10, force = 5) +
    # scale_size_identity() + scale_alpha(guide = 'none') + 
    # theme(# axis.text = element_text(size = 20), legend.text = element_text(size = 20),
    #       text = element_text(size = 17)) 
  
  return(list(p, diffexp_res_annotated))
}

