library(arrow)
library(limma)
library(dplyr)

# working directory should be: aml_outcome_disparities/analysis/mutations/
# 1. Load Parquets
data <- read_parquet("./_cache/combined_omics_data.parquet")
contrasts_df <- read_parquet("./_cache/limma/contrast.parquet")
design_df <- read_parquet("./_cache/limma/design.parquet")

# 2. Safely format the Contrast Matrix
# Handle whether index is stored in a column or in arrow metadata
if ("Group" %in% colnames(contrasts_df)) {
  rownames_contrasts <- contrasts_df$Group
  contrasts_df <- contrasts_df %>% select(-Group)
} else if ("__index_level_0__" %in% colnames(contrasts_df)) {
  rownames_contrasts <- contrasts_df$`__index_level_0__`
  contrasts_df <- contrasts_df %>% select(-`__index_level_0__`)
} else {
  rownames_contrasts <- rownames(contrasts_df)
}

contrasts <- as.matrix(contrasts_df)
contrasts <- apply(contrasts, 2, as.numeric)
rownames(contrasts) <- rownames_contrasts


# 3. Safely format the Design Matrix
if ("Sample" %in% colnames(design_df)) {
  rownames_design <- design_df$Sample
  design_df <- design_df %>% select(-Sample)
} else if ("__index_level_0__" %in% colnames(design_df)) {
  rownames_design <- design_df$`__index_level_0__`
  design_df <- design_df %>% select(-`__index_level_0__`)
} else {
  rownames_design <- rownames(design_df)
}

design_full <- as.matrix(design_df)
# FORCE strictly numeric (float/double) to eliminate strings/booleans
mode(design_full) <- "numeric"
rownames(design_full) <- rownames_design

results <- list()

for (om in unique(data$OmicsType)) {
  # 1. Filter to this omics layer
  sub_data <- data %>% filter(OmicsType == om)
  
  # 2. Extract expression matrix (Features x All Samples)
  sample_cols <- intersect(colnames(sub_data), rownames(design_full))
  expr_mat <- as.matrix(sub_data[, sample_cols])
  rownames(expr_mat) <- sub_data$Feature
  
  # 3. Find samples that ACTUALLY have data for this omics layer
  # (drop samples that are entirely NA for RNA)
  valid_samples <- sample_cols[colSums(!is.na(expr_mat)) > 0]
  
  # 4. Subset both expression and design to the active samples
  expr_sub <- expr_mat[, valid_samples]
  design_sub <- design_full[valid_samples, , drop = FALSE]
  
  # 5. Sanity check: drop design columns (groups) with 0 samples in this omics layer
  active_groups <- colSums(design_sub) > 0
  if (!all(active_groups)) {
    warning(paste(
      "Omics", om, "is missing samples for groups:", 
      paste(names(active_groups)[!active_groups], collapse = ", ")
    ))
  }
  
  # 6. Run limma
  fit <- lmFit(expr_sub, design_sub)
  fit_contrast <- contrasts.fit(fit, contrasts)
  fit_bayes <- eBayes(fit_contrast)
  
  # 7. Collect results
  for (c_name in colnames(contrasts)) {
    tt <- topTable(fit_bayes, coef = c_name, number = Inf, sort.by = "none")
    tt$Feature <- rownames(tt)
    tt$OmicsType <- om
    tt$Contrast <- c_name
    results[[paste(om, c_name, sep = "_")]] <- tt
  }
}

write_parquet(bind_rows(results), "./_cache/limma/results.parquet")