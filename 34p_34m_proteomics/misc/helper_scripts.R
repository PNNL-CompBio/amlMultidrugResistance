library(MSnSet.utils)
library(dplyr)
library(KSEAapp)
library(ggplot2)
library(scales)
library(ggpubr)

load_ksdb <- function(kind = "PSP_networkin_2016", output_as = "original", networkin_cutoff = Inf, organism = "human"){
   if (kind == "PSP_networkin_2016"){
      # KSDB <- read.csv('C:/Users/poss982/Documents/Github/camilo_helper/PSP&NetworKIN_Kinase_Substrate_Dataset_July2016.csv',
      #                  stringsAsFactors = FALSE)
      KSDB <- read.csv(syn$get("syn72663630")$path, stringsAsFactors = FALSE)
   } else if (kind == "PSP May 23 2023"){
      # KSDB <- read.table("C:/Users/poss982/Documents/GitHub/camilo_helper/PSP_Kinase_Substrate.txt",
      #                    sep = "\t", header = T) %>%syn72663630
      KSDB <- read.table(syn$get("syn72663629")$path, sep = "\t", header = T) %>%
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



