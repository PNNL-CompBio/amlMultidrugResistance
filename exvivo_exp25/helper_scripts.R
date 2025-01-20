library(MSnSet.utils)
library(dplyr)
library(KSEAapp)
library(ggplot2)
library(scales)
library(ggpubr)

load_ksdb <- function(kind = "PSP_networkin_2016", output_as = "original", networkin_cutoff = Inf, organism = "human"){
  if (kind == "PSP_networkin_2016"){
    KSDB <- read.csv('C:/Users/poss982/Documents/Github/camilo_helper/PSP&NetworKIN_Kinase_Substrate_Dataset_July2016.csv',
                     stringsAsFactors = FALSE)
  } else if (kind == "PSP May 23 2023"){
    KSDB <- read.table("C:/Users/poss982/Documents/GitHub/camilo_helper/PSP_Kinase_Substrate.txt", 
                       sep = "\t", header = T) %>%
      mutate(Source = "PhosphoSitePlus",
             networkin_score = "Inf") %>%
      filter(KIN_ORGANISM == organism,
             SUB_ORGANISM == organism)
    KSDB <- KSDB[, c("KINASE", "KIN_ACC_ID", "GENE", "KIN_ORGANISM", "SUBSTRATE", "SUB_GENE_ID", "SUB_ACC_ID", "SUB_GENE",
                     "SUB_ORGANISM", "SUB_MOD_RSD", "SITE_GRP_ID", "SITE_...7_AA", "networkin_score", "Source")]
  }
  if (output_as == "ora"){
    KSDB <- KSDB %>%
      filter(networkin_score >= networkin_cutoff) %>%
      mutate(gene_symbol = paste(SUB_GENE, SUB_MOD_RSD, sep = "-")) %>%
      select(gs_name = GENE,
             gene_symbol)
  }
  return(KSDB)
}


GSEA_helper_single <- function(m, contrasts, t2g, t2g_name, coef.str){
  if (!all(grepl(paste0("^", coef.str), contrasts))){
    print("Adding coef.str to contrast names")

    group1 <- sub("-.*$", "", contrasts)
    group2 <- sub("^.*-", "", contrasts)
    
    l_contrasts <- paste0(coef.str, group1, "-", coef.str, group2)
  } else {
    l_contrasts = contrasts
  }
  limma_res <- limma_contrasts(m, model.str = paste("~0 +", coef.str), 
                               coef.str = coef.str, contrasts = l_contrasts)
  
  plot_title <- gsub(coef.str, "", contrasts) %>% gsub("_", " ", .) %>% gsub("-", " vs ", .)
  group1 <- sub(" vs .*$", "", plot_title)
  group2 <- sub("^.* vs ", "", plot_title)
  plot_path <- gsub(" ", "_", plot_title) %>% paste0("GSEA_", ., paste0("_", t2g_name, ".png"))
  tbl_path <- sub(".png", ".txt", plot_path)
  
  ## For some pathways, gsea is unable to assess a pvalue. These are removed by the function
  ## Internally. This is why some pathways (like "GOBP NCRNA METABOLISM") are excluded from
  ## the result.
  
  if (file.exists(tbl_path)){
    print("using saved results")
    gsea_res <- read.table(tbl_path, sep = "\t")
  } else {
    fold_change <- limma_res$logFC
    names(fold_change) <- limma_res$feature
    fold_change <- sort(fold_change, decreasing = TRUE)
    set.seed(69)
    
    gsea_res <- clusterProfiler::GSEA(fold_change, eps = 1e-16, minGSSize = 10, 
                     pvalueCutoff = 1, TERM2GENE = t2g)@result %>%
      dplyr::select(Description, setSize, NES, pvalue, p.adjust, core_enrichment) %>%
      dplyr::mutate(contrast = contrasts)
    write.table(x = gsea_res, file = tbl_path, 
                sep = "\t", quote = F)
  }
}


GSEA_helper <- function(m, contrasts, t2g, t2g_name, coef.str){
  for (contrast in contrasts){
    GSEA_helper_single(m, contrast, t2g, t2g_name, coef.str)
  }

  file_pattern <- paste0("^GSEA_.*_vs_.*_", t2g_name, ".txt$")
  combined <- data.frame()
  for (file_path in list.files(pattern = file_pattern)){
    xx <- read.table(file_path, sep = "\t")
    combined <- rbind(combined, xx)
  }
  write.table(x = combined, file = paste0("GSEA_", t2g_name, "_combined.txt"), 
              sep = "\t", quote = F)
}




KSEA_helper_single <- function(m, contrasts, coef.str, prefix = "", psp_db = "PSP_networkin_2016", ...){
  plot_title <- gsub(coef.str, "", contrasts) %>% gsub("_", " ", .) %>% gsub("-", " vs ", .)
  plot_title <- paste(prefix, plot_title)
  group1 <- sub(" vs .*$", "", plot_title)
  group2 <- sub("^.* vs ", "", plot_title)
  plot_path <- gsub(" ", "_", plot_title) %>% paste0("KSEA_", ., paste0("_", "NetworKIN_Inf", ".png"))
  tbl_path <- sub(".png", ".txt", plot_path)
  
  KSDB <- load_ksdb(psp_db, ...)
  
  if (!all(grepl(paste0("^", coef.str), contrasts))){
    print("Adding coef.str to contrast names")
    
    group1 <- sub("-.*$", "", contrasts)
    group2 <- sub("^.*-", "", contrasts)
    
    l_contrasts <- paste0(coef.str, group1, "-", coef.str, group2)
  } else {
    l_contrasts = contrasts
  }
  
  limma_res <- limma_contrasts(m, model.str = paste("~0 +", coef.str), 
                               coef.str = coef.str, contrasts = l_contrasts)
  
  fold_change <- limma_res$logFC
  fold_change <- 2**fold_change
  PX <- data.frame(Protein = "NULL", Gene = limma_res$feature, Peptide = "NULL", 
                   Residue.Both = limma_res$feature, p = "NULL", FC = fold_change) %>%
    dplyr::mutate(Residue.Both = sub("^.*-", "", Residue.Both)) %>%
    dplyr::mutate(Residue.Both = gsub("[a-z]", ";", Residue.Both)) %>%
    dplyr::mutate(Residue.Both = gsub(";$", "", Residue.Both),
                  Gene = sub("^(.*)-.*$", "\\1", Gene))
  
  ksea_res <- KSEA.Scores(KSDB, PX, NetworKIN = TRUE, NetworKIN.cutoff = Inf) %>%
    dplyr::select(Kinase.Gene, m, FDR, z.score) %>%
    dplyr::rename(pathway = Kinase.Gene, enrichment = z.score,
                  adj_p_val = FDR, set_size = m) %>%
    dplyr::mutate(contrast = contrasts) %>%
    filter(set_size >= 3)
  write.table(x = ksea_res, file = tbl_path, 
              sep = "\t", quote = F)
  
}


KSEA_helper <- function(m, contrasts, coef.str, prefix = "", ...){
  for (contrast in contrasts){
    print(contrast)
    KSEA_helper_single(m, contrast, coef.str, prefix, ...)
  }
  
  file_pattern <- paste0("^KSEA_", prefix, ".*_vs_.*_", "NetworKIN_Inf", ".txt$")
  combined <- data.frame()
  for (file_path in list.files(pattern = file_pattern)){
    xx <- read.table(file_path, sep = "\t")
    combined <- rbind(combined, xx)
  }
  write.table(x = combined, file = paste0("KSEA_", prefix, "_NetworKIN_Inf", "_combined.txt"), 
              sep = "\t", quote = F)
}


site_splitter <- function(m){
  f_data <- fData(m) %>%
    mutate(og_site = rownames(.)) %>%
    mutate(site_end = sub("^.*-([A-Za-z0-9]+)$", "\\1", og_site),
           site_base = sub("-[A-Za-z0-9]+$", "", og_site))
  
  f_data_expand <- data.frame()
  
  for (og_site_i in f_data$og_site){
    endings <- str_split(f_data[og_site_i, "site_end"], "[a-z]")[[1]] %>% head(-1)
    f_data_expand <- f_data_expand %>% 
      rbind(data.frame(og_site = og_site_i, split_site = paste(f_data[og_site_i, "site_base"], endings, sep = "-")))
  }
  
  m_long <- exprs(m) %>%
    as.data.frame() %>% mutate(og_site = rownames(.)) %>%
    tidyr::pivot_longer(-og_site, names_to = "Sample", values_to = "value") %>%
    filter(!is.na(value)) %>%
    left_join(f_data_expand, by = "og_site") %>%
    select(split_site, Sample, value) %>%
    mutate(value = as.numeric(value))
  
  new_mat <- tidyr::pivot_wider(m_long, values_from = "value", names_from = "Sample", values_fn = mean) %>%
    as.data.frame()
  rownames(new_mat) <- new_mat$split_site
  new_mat <- new_mat[, -1]
  new_mat <- new_mat[, colnames(exprs(m))]
  
  f_data_split <- f_data_expand %>%
    group_by(split_site) %>% summarize(og_sites = list(unlist(og_site))) %>% as.data.frame()
  rownames(f_data_split) <- f_data_split$split_site
  
  m_split <- MSnSet(exprs = as.matrix(new_mat), pData = pData(m), fData = f_data_split[rownames(new_mat), ])
  return(m_split)
}




## Modification of MSnSetUtils function.
plot_pca <- function(eset, phenotype = NULL, shape = NULL, label = NULL, z_score = TRUE,
                     princomp_center = TRUE, show_ellipse = TRUE, components = 1:2, biplot = FALSE,
                     biplot_labels = NULL, standardize = TRUE, save_dfs = NULL,
                     num_features = 6L, show_NA = TRUE, label_size = 3,
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
  
  return(p)
}


arrow_plotter <- function(arrow_df, top_n = 10, text_size = 4){
  
  arrow_df$arrow_len <- sqrt(arrow_df$x_c**2 + arrow_df$y_c**2)
  arrow_df <- arrow_df %>% arrange(-arrow_len) %>% head(top_n)
  range <- max(abs(arrow_df$x_c), abs(arrow_df$y_c))
  
  p <- ggplot(arrow_df, aes(xend = x_c, yend = y_c)) + 
    geom_segment(x = 0, y = 0, arrow = arrow(angle = 25, length = unit(0.25, "cm"))) + 
    ggrepel::geom_text_repel(arrow_df, mapping = aes(x = x_c, y = y_c, label = label_name), 
                    nudge_x = 0.005, nudge_y = 0.005, box.padding = 0.5, 
                    segment.color = "black", segment.alpha = 1, size = text_size) +
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position="none",
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank()) + xlim(-range, range) + ylim(-range, range)
  
  return(p)
}


ORA_helper <- function(diffexp, universe, t2g, adj_pval_cutoff = 0.05, logFC_cutoff = 0, override_phospho_check = F){
  results <- data.frame()
  
  if (all(grepl("^.*-[A-Z]+[0-9]+.*$", diffexp$feature)) & !override_phospho_check){
    print("Detected phospho data, splitting sites.\n")
    
    unique_features <- unique(diffexp$feature)
    m_dummy <- MSnSet(exprs = matrix(unique_features, nrow = length(unique_features), ncol = 2, 
                                     dimnames = list(unique_features, c("dummy_1", "dummy_2"))), 
                      pData = data.frame("dummy" = c("dummy_1", "dummy_2"), 
                                         row.names = c("dummy_1", "dummy_2")))
    print("splitting sites\n")
    m_dummy <- site_splitter(m_dummy)
    f_data <- tidyr::unnest(fData(m_dummy)) %>% as.data.frame() %>% select(feature = og_sites, split_site)
    diffexp <- merge(diffexp, f_data, by = "feature") %>%
      group_by(split_site) %>%
      top_n(1, wt = -adj_pval) %>%
      select(-feature) %>%
      select(feature = split_site, everything())
    universe = diffexp$feature
  }
  
  for (contrast_i in unique(diffexp$contrast)){
    print(contrast_i)
    feature_df <- diffexp %>% filter(contrast == contrast_i,
                                    adj_pval < adj_pval_cutoff,
                                    abs(logFC) >= logFC_cutoff) %>%
      mutate(sign = sign(logFC)) %>% select(feature, sign)
      
    xx_positive <- try(clusterProfiler::enricher(feature_df %>% filter(sign == 1) %>% pull(feature), 
                                             universe = universe, TERM2GENE = t2g,
                                             maxGSSize = 500, pvalueCutoff = 0.05)@result %>%
      mutate(sign = 'positive', contrast = contrast_i))
    
    xx_negative <- try(clusterProfiler::enricher(feature_df %>% filter(sign == -1) %>% pull(feature), 
                                             universe = universe, TERM2GENE = t2g,
                                             maxGSSize = 500, pvalueCutoff = 0.05)@result %>%
      mutate(sign = 'negative', contrast = contrast_i))
    
    xx_both <- try(clusterProfiler::enricher(feature_df %>% pull(feature), 
                                             universe = universe, TERM2GENE = t2g,
                                             maxGSSize = 500, pvalueCutoff = 0.05)@result %>%
                     mutate(sign = 'both', contrast = contrast_i))
    
    if (class(xx_positive) == "try-error"){
      xx_positive <- data.frame()
    }
    if (class(xx_negative) == "try-error"){
      xx_negative <- data.frame()
    }
    if (class(xx_both) == "try-error"){
      xx_both <- data.frame()
    }
    results <- rbind(results, xx_positive) %>%
      rbind(xx_negative) %>% rbind(xx_both)
  }
  
  return(results)
}

# ann_colors <- list(c("darkorange2", "firebrick2", "mediumpurple", "dodgerblue3", "darkolivegreen4"),
#                    c("gray99", "gray42", "gray69"))

plot_heatmap <- function(m, features, annotation_cols, ann_colors = NA, heatmap_title = "", ...){
  meta_df <- pData(m)[, annotation_cols] %>% as.data.frame()
  colnames(meta_df) <- annotation_cols
  rownames(meta_df) <- rownames(pData(m))
  mat <- exprs(m)
  
  mat <- mat[intersect(features, rownames(mat)), ]
  mat <- sweep(mat, 1, apply(mat, 1, mean, na.rm = T), FUN = '-')
  
  if (any(!is.na(ann_colors))){
     for (index in 1:ncol(meta_df)){
        if (class(meta_df[[index]]) != "numeric"){
           # ann_col = colnames(meta_df)[[index]]
           # names(ann_colors[[index]]) <- unique(meta_df[[index]])
        }
     }
     names(ann_colors) <- colnames(meta_df)
  }
  
  
  color_pal <- rev(RColorBrewer::brewer.pal(n = 9, name = "RdYlBu"))
  color_pal <- c(colorRampPalette(color_pal[1:4])(15), colorRampPalette(color_pal[5:9])(15))
  y_max = (abs(min(mat, na.rm = T)) + abs(max(mat, na.rm = T)))/1.7
  breaks_buddy <- c(seq(-y_max, -0.001, length.out = 15),
                    0,
                    seq(0.001, y_max, length.out = 15))
  
  p <- pheatmap::pheatmap(mat, annotation_col = meta_df, breaks = breaks_buddy,
                          color = color_pal, annotation_colors = ann_colors, 
                          main = heatmap_title, ...)
  return(p)
}

plot_pathway_heatmap <- function(m, pathway, t2g, annotation_cols, ann_colors = NA, ...){
  features <- t2g %>%
    filter(gs_name == pathway) %>%
    pull(gene_symbol) %>% unique()
  
  plot_heatmap(m, features, annotation_cols, ann_colors, heatmap_title = pathway, ...)
}  

plot_diffexp_genes <- function(diffexp, features, sig_cutoff = 0.05, 
                               sig_symbol = "*", heatmap_title = "Diffexp heatmap", ...){
  mat <- diffexp %>%
    filter(feature %in% features) %>%
    select(feature, logFC, contrast) %>%
    tidyr::pivot_wider(names_from = "contrast", values_from = logFC) %>% as.data.frame()
  rownames(mat) <- mat$feature
  mat <- as.matrix(mat %>% select(-feature))
  
  mat_pval <- diffexp %>%
    filter(feature %in% features) %>%
    select(feature, adj_pval, contrast) %>%
    tidyr::pivot_wider(names_from = "contrast", values_from = adj_pval) %>% as.data.frame()
  rownames(mat_pval) <- mat_pval$feature
  mat_pval <- as.matrix(mat_pval %>% select(-feature))
  mat_pval[mat_pval < sig_cutoff] <- sig_symbol
  mat_pval[mat_pval >= sig_cutoff] <- ""
  mat_pval[is.na(mat_pval)] <- ""
  
  color_pal <- rev(RColorBrewer::brewer.pal(n = 9, name = "RdYlBu"))
  color_pal <- c(colorRampPalette(color_pal[1:4])(15), colorRampPalette(color_pal[5:9])(15))
  y_max = (abs(min(mat, na.rm = T)) + abs(max(mat, na.rm = T)))/1.7
  breaks_buddy <- c(seq(-y_max, -0.001, length.out = 15),
                    0,
                    seq(0.001, y_max, length.out = 15))
  
  p <- pheatmap::pheatmap(mat, breaks = breaks_buddy, display_numbers = mat_pval,
                          angle_col = 45, color = color_pal, main = heatmap_title, ...)
  return(p)
}

plot_diffexp_pathway <- function(diffexp, pathway, t2g, sig_cutoff = 0.05, 
                                 sig_symbol = "*", heatmap_title = "Diffexp heatmap", ...){
  features <- t2g %>%
    filter(gs_name == pathway) %>%
    pull(gene_symbol) %>% unique()
  pathway_title <- paste("Differential expression of", pathway)
  
  plot_diffexp_genes(diffexp, features, sig_cutoff, sig_symbol,
                     heatmap_title = pathway_title, ...)
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
      counter = counter + 1
    }
    
    limma_res <- limma_res %>%
      mutate(t_test_pval = as.numeric(p_values_t_test),
             t_test_adj = p.adjust(t_test_pval, method = "BH"),
             welch_pval = as.numeric(p_values_welch_test),
             welch_adj = p.adjust(welch_pval, method = "BH"))
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
      geom_boxplot() + theme_bw() + theme(axis.text.x = text_element) + 
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




library(dplyr)
library(ggplot2)
library(gridExtra)

## enrichment_df needs a pathway column, adj_p_val column, and enrichment column (NES from GSEA). Enrichment is expected to be a continuous variable.
## top is the number of pathways (up regulated and down regulated) to take.

plot_enrichment_bars <- function(enrichment_result, enrichment_label, top = 10,
                                 adj_pval_cutoff = 0.05, enrich_width = 9, sig_width = 4,
                                 enrichment_title = "Enrichment"){
   
   #' Used to make reversed logarithmic scales
   #' @import scales
   reverselog_trans <- function(base = exp(1)) {
      library(scales)
      trans <- function(x) -log(x, base)
      inv <- function(x) base^(-x)
      scales::trans_new(paste0("reverselog-", format(base)), trans, inv,
                        breaks = scales::log_breaks(base = base),
                        domain = c(1e-100, Inf))
   }
   
   
   plot_df <- enrichment_result %>%
      dplyr::select(pathway, adj_p_val, enrichment) %>%
      filter(adj_p_val < adj_pval_cutoff) %>%
      mutate(sign = case_when(enrichment > 0 ~ "Up",
                              TRUE ~ "Down")) %>%
      group_by(sign) %>%
      arrange(adj_p_val) %>%
      slice(1:top) %>%
      ungroup()
   
   p_enrichment <- ggplot(plot_df, aes(x = enrichment, y = reorder(pathway, enrichment))) +
      geom_bar(stat = 'identity', aes(fill = sign)) +
      scale_fill_manual(values = c("Down" = "dodgerblue3", "Up" = "firebrick2", "Not significant" = "black")) +
      scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40)) +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
            legend.position = "none",
            axis.title.x = element_text(size=16),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 12),
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank()) +
      labs(x = enrichment_label) +
      ggtitle(enrichment_title)
   
   
   p_sig <- ggplot(plot_df, aes(x = adj_p_val, y = reorder(pathway, enrichment))) +
      geom_bar(stat='identity') +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
            legend.position="none",
            axis.title.x = element_text(size=16),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = 14),
            axis.ticks.y = element_blank(),
            axis.line.y = element_line(color = "black"),
            axis.text.y = element_blank()) +
      scale_x_continuous(trans = reverselog_trans(10)) +
      labs(x = "Adjusted p-value") +
      ggtitle("Significance")
   
   arrange_matrix <- as.matrix(c(rep(1, enrich_width), rep(2, sig_width))) %>% t()
   p_both <- grid.arrange(p_enrichment, p_sig, layout_matrix = arrange_matrix)
}
