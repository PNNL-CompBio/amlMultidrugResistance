#!/usr/local/bin/Rscript --vanilla
# MomDiablo.R
# Dylan Ross (dylan.ross@pnnl.gov)
# Description:
#   Run DIABLO analysis from mixOmics

# nolint start: object_name_linter 

library(argparser, quietly = TRUE)
library(magrittr, quietly = TRUE)
suppressMessages(library(mixOmics, quietly = TRUE))


# sets up an argument parser with description and expected arguments
.SetupParser <- function() {
  p <- arg_parser(
    "Perform DIABLO analysis on multiomic data using mixOmics",
    name = "MomDiablo",
    hide.opts = TRUE
  )
  p <- add_argument(
    p,
    "cmd",
    "specify whether to tune parameters for DIABLO analysis ('tune') or run the final DIABLO analysis ('final')",
    nargs = 1
  )
  p <- add_argument(
    p,
    "--block-csvs",
    "input data files in CSV format (one for each block, all should have the same columns)",
    nargs = Inf,  # accept indefinite number of files
    short = "-bc"
  )
  p <- add_argument(
    p,
    "--block-names",
    "block names (one for each block)",
    nargs = Inf,  # accept indefinite number of files
    short = "-bn"
  )
  p <- add_argument(
    p,
    "--design-csv",
    "input file with design matrix, rows/columns should be in same order as specified blocks",
    nargs = 1
  )
  p <- add_argument(
    p,
    "--target-csv",
    "a CSV with target classification label for all samples (rows correspond to columns in block data)",
    nargs = 1,
    short = "-t"
  )
  p <- add_argument(
    p,
    "--target-csv-column",
    "Select the column from the CSV with the target classification label",
    nargs = 1,
    short = "-tc"
  )
  p <- add_argument(
    p,
    "--plot-dir",
    default = "./_figures/",
    "directory to save all plots under",
    nargs = 1,
    short = "-plt"
  )
  p <- add_argument(
    p,
    "--distance-metric",
    "distance metric (must be one of 'max.dist', 'centroids.dist', 'mahalanobis.dist', used in 'final' analysis)",
    nargs = 1,
    short = "-m"
  )
  p <- add_argument(
    p,
    "--n-components",
    "number of components for DIABLO (used in 'final' analysis)",
    nargs = 1,
    short = "-n",
    type = "integer"
  )
  p <- add_argument(
    p,
    "--keepx-csv",
    "a CSV with the count of features to keep for each component from each block (used in 'final' analysis)",
    nargs = 1,
    short = "-k"
  )
  return(p)
}


# load all of the block data from CSV files,
# transpose to shape: samples x features,
# combine into a list with specified block names
.LoadBlocks <- function(block_files, block_names) {
  block_data <- setNames(block_files, block_names)
  block_data %<>%
    lapply(read.csv, header = TRUE, row.names = 1) %>%
    lapply(t)
  return(block_data)
}


# load the target variable from a specified column in a specified CSV
# the rows are samples, like the columns from the block data files
.LoadTarget <- function(target_file, target_column) {
  target <- factor(
    read.csv(
      target_file,
      header = TRUE
    )[[target_column]]
  )
  return(target)
}


# load a design matrix from file to define relationships between blocks
.LoadDesign <- function(design_file) {
  design <- as.matrix(
    read.csv(
      design_file,
      header = TRUE,
      row.names = 1
    )
  )
  return(design)
}


# run analysis parameter tuning
# determine the number of components to use, distance metric,
# and how many features to keep from each dataset
.RunTuning <- function(argv) {
  cat("TUNING DIABLO MODEL", "\n")
  # load data
  X <- .LoadBlocks(argv$block_csvs, argv$block_names)
  cat("loaded blocks", "\n")
  print(names(X))
  Y <- .LoadTarget(argv$target_csv, argv$target_csv_column)
  cat("loaded targets", "\n")
  design <- .LoadDesign(argv$design_csv)
  cat("loaded design", "\n")
  # set up the DIABLO analysis
  suppressMessages(
    diablo <- block.plsda(X, Y,
                          ncomp = 8,
                          design = design,
                          near.zero.var = TRUE)
  )
  # characterize the performance
  suppressMessages(
    perf.diablo <- perf(diablo,
                        validation = 'Mfold',
                        folds = 5,
                        nrepeat = 10,)
  )
  # plot the results
  png(paste0(argv$plot_dir, "diablotune_perf.png"),
      width = 4.5, height = 6., units = "in", res = 300)
  plot(perf.diablo)
  dev.off()
  # run the model tuning
  test.keepX <- rep(list(c(8, 16, 32)), length(argv$block_names))
  names(test.keepX) <- argv$block_names
  # perform tuning with all three distance metrics
  for (dist in c("max.dist", "centroids.dist", "mahalanobis.dist")) {
    cat("--------------------", "\n")
    cat("dist: ", dist, "\n")
    # get and report the optimal number of components by weighted vote using this metric
    ncomp <- perf.diablo$choice.ncomp$WeightedVote["Overall.BER", dist]
    cat("ncomp: ", ncomp, "\n")
    # run the tuning to see how many features to keep from each block for each component
    suppressMessages(
      tune.diablo <- tune.block.splsda(X, Y,
                                       ncomp = ncomp,
                                       test.keepX = test.keepX,
                                       design = design,
                                       validation = 'Mfold',
                                       folds = 5,
                                       nrepeat = 5,
                                       dist = dist,
                                       near.zero.var = TRUE,
                                       BPPARAM = BiocParallel::SnowParam(workers = 10))
    )
    # report the number of selected features for each component
    cat("list.keepX:", "\n")
    print(tune.diablo$choice.keepX)
  }
  cat("--------------------", "\n")
}


# load the number of features to keep for each component from each block
.LoadKeepX <- function(keepx_file) {
  df <- read.csv(
    keepx_file,
    header = FALSE,
    row.names = 1
  )
  keep_x <- lapply(split(df, row.names(df)), unlist)
  return(keep_x)
}


# runs the final DIABLO analysis and produces all the plots
.RunFinal <- function(argv) {
  cat("FINAL DIABLO ANALYSIS", "\n")
  cat("--------------------", "\n")
  # load data
  X <- .LoadBlocks(argv$block_csvs, argv$block_names)
  Y <- .LoadTarget(argv$target_csv, argv$target_csv_column)
  design <- .LoadDesign(argv$design_csv)
  keepX <- .LoadKeepX(argv$keepx_csv)
  # fit the final DIABLO model
  suppressMessages(
    final.diablo <- block.splsda(X, Y,
                                 ncomp = argv$n_components,
                                 keepX = keepX,
                                 design = design,
                                 near.zero.var = TRUE)
  )
  # plot the DIABLO results
  # arrow plot
  # correlation for all components
  for (ncomp in 1:argv$n_components) {
    png(paste0(argv$plot_dir, "diablofinal_comp", ncomp, "_correlation.png"),
        width = 12, height = 12, units = "in", res = 300)
    plotDiablo(final.diablo, ncomp = ncomp)
    dev.off()
  }
  # projections from all blocks
  png(paste0(argv$plot_dir, "diablofinal_block_projections.png"),
      width = 12, height = 12, units = "in", res = 300)
  plotIndiv(final.diablo, ind.names = FALSE, legend = TRUE)
  dev.off()
  # correlation circle plot
  png(paste0(argv$plot_dir, "diablofinal_corr_circle.png"),
      width = 8, height = 8, units = "in", res = 300)
  plotVar(final.diablo, var.names = FALSE, style = 'graphics', legend = TRUE)
  dev.off()
  # loadings for each component
  for (ncomp in 1:argv$n_components) {
    for (block in argv$block_names) {
      png(paste0(argv$plot_dir, "diablofinal_comp", ncomp, "_block-", block, "_loadings.png"),
          width = 7, height = 4, units = "in", res = 300)
      plotLoadings(final.diablo,
                  comp = ncomp,
                  contrib = "max",
                  method = "median",
                  names.var = NULL,
                  block = block,
                  size.title = 1)
      dev.off()
    }
  }
  block_colors <- c("#8AE0AE", "#ADE558", "#7FCFDE", "#B954E0", "#DC6780", "#AF8ED6")
  # truncate block colors to the same length as the number of blocks
  n_blocks <- length(argv$block_names)
  block_colors <- block_colors[1:n_blocks]
  # arrow plot
  png(paste0(argv$plot_dir, "diablofinal_arrows.png"),
      width = 8, height = 8, units = "in", res = 300)
  # NOTE: This might fail silently, I have no clue why.
  plotArrow(final.diablo, ind.names = FALSE, legend = TRUE)
  # circos plot, all components then individual components
  png(paste0(argv$plot_dir, "diablofinal_circos.png"),
      width = 8, height = 8, units = "in", res = 300)
  circosPlot(final.diablo, cutoff = 0.8, line = TRUE,
             #color.blocks = rep("#606060", length(argv$block_names)),
             color.blocks = block_colors,
             color.cor = c("darkorchid", "lightgreen"), size.labels = 1.5)
  dev.off()
  for (ncomp in 1:argv$n_components) {
    png(paste0(argv$plot_dir, "diablofinal_comp", ncomp, "_circos.png"),
        width = 8, height = 8, units = "in", res = 300)
    suppressWarnings(
      circosPlot(final.diablo, comp = c(ncomp), cutoff = 0.8, line = TRUE,
                 color.blocks = block_colors,
                 color.cor = c("darkorchid", "lightgreen"), size.labels = 1.5)
    )
    dev.off()
  }
  # cluster heatmap
  png(paste0(argv$plot_dir, "diablofinal_clust_heatmap.png"),
      width = 12, height = 12, units = "in", res = 300)
  invisible(capture.output(
    cimDiablo(final.diablo, 
              color.blocks = block_colors,
              comp = 1, margin=c(8,20), legend.position = "right")
  ))
  dev.off()
  # compute final performance metrics
  perf.final.diablo <- perf(final.diablo,  validation = 'Mfold', folds = 5,
                            nrepeat = 10, dist = "all")

  # plot results
  png(paste0(argv$plot_dir, "diablofinal_perf.png"),
      width = 5, height = 5, units = "in", res = 300)
  plot(perf.final.diablo)
  dev.off()
  # compute and report confusion matrix
  predict.diablo <- predict(final.diablo, newdata = X)
  predicted <- predict.diablo$WeightedVote[[argv$distance_metric]][, 1]
  confusion_mat <- get.confusion_matrix(truth = Y,
                                        predicted = predicted)
  print(confusion_mat)
  cat("BER: ", get.BER(confusion_mat), "\n")
}


.Main <- function() {
  # create a parser and parse the args
  argv <- parse_args(.SetupParser())
  # proceed with either tuning or final DIABLO analysis
  switch(argv$cmd,
         "tune" = .RunTuning(argv),
         "final" = .RunFinal(argv),
         stop("Unknown `cmd` (should be 'tune' or 'final')", call. = FALSE)
  )

}

# nolint end: object_name_linter 

# run the main function
.Main()
