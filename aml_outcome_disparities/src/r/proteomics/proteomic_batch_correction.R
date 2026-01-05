library(Biostrings)
library(dplyr)
library(here)
library(MSnSet.utils)
library(synapser)

here::i_am("src/r/proteomics/proteomic_batch_correction.R")
# Assumes an authorization token for Synapse is stored in the root directory
# with name 'auth_token.txt'
token <- readLines(here("auth_token.txt"), warn = FALSE)
synapser::synLogin(authToken = token)

#' Proteomic batch correction
#'
#' Corrects phospho & global proteomic measurements via ComBat.
#' @param phospho Boolean. TRUE processes phospho measurements; FALSE processes
#' global measurements.
#' @param output_path Path to save processed data. Defaults to 
#' "src/python/data".
#' @return NULL (invisibly).
combat_proteomics <- function(phospho = TRUE, output_path = NULL) {
  # Setup storage path
  if (is.null(output_path)) {
    output_path <- here(
      "src",
      "python",
      "data"
    )
  }
  
  # Import data from Synapse
  if (phospho) {
    ba <- t(read.csv(
      synapser::synGet("syn25714936")$path,
      sep = "\t"
    ))
    pilot <- t(read.csv(
      synapser::synGet("syn69075545")$path,
      sep = "\t"
    ))
  }
  else {
    ba <- t(read.csv(
      synapser::synGet("syn25714254")$path,
      sep = "\t"
    ))
    pilot <- t(read.csv(
      synapser::synGet("syn69075554")$path,
      sep = "\t"
    ))
  }
  
  # Trim to measurements in both datasets
  shared <- intersect(
    colnames(ba),
    colnames(pilot)
  )
  ba <- ba[,shared]
  pilot <- pilot[,shared]
  
  # Concatenate, define batches
  data <- t(rbind(
    ba,
    pilot
  ))
  batch <- integer(ncol(data))
  batch[-1:-nrow(ba)] <- 1
  
  # Correct using PNNL's missing-friendly ComBat
  corrected <- ComBat.NA(
    data,
    batch,
    par.prior = TRUE
  )
  corrected <- corrected$`corrected data`
  
  if (phospho) {
    file_name <- "phospho"
  }
  else {
    file_name <- "global"
  }
  
  # Separate Pilot from BeatAML samples
  ba_corrected <- corrected[,row.names(ba)]
  pilot_corrected <- corrected[,row.names(pilot)]
  
  # Write data
  write.table(
    ba_corrected,
    here(
      output_path,
      sprintf("ba_%s_corrected.csv", file_name)
    ),
    quote=F,
    sep=","
  )
  write.table(
    pilot_corrected,
    here(
      output_path,
      sprintf("pilot_%s_corrected.csv", file_name)
    ),
    quote=F,
    sep=","
  )
}

combat_proteomics(phospho = TRUE)
combat_proteomics(phospho = FALSE)
