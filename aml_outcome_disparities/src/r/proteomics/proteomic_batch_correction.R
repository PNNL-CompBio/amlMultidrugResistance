library(Biostrings)
library(dplyr)
library(here)
library(MSnSet.utils)

here::i_am("src/r/proteomics/proteomic_batch_correction.R")

#' Proteomic batch correction
#'
#' Corrects phospho & global proteomic measurements via ComBat.
#' @param ba DataFrame. Measurements for BeatAML cohort.
#' @param pilot DataFrame. Measurements for Pilot cohort.
#' @param meta DataFrame. Meta-data for patients. Row names should be same as
#' those in ba and pilot DataFrames. Must contain columns 'Race' and 'Study'.
#' @return [list] with the following elements:
#' \itemize{
#'  \item ba: Batch-corrected data for BeatAML cohort
#'  \item pilot: Batch-corrected data for Pilot cohort.
#' }
combat_proteomics <- function(ba, pilot, meta) {
  # Group non-Black/non-White patients, given limited sample size
  meta[!(meta$Race %in% c("Black", "White")), "Race"] <- "Other"
  
  # Trim to measurements in both datasets
  shared <- intersect(
    colnames(ba),
    colnames(pilot)
  )
  ba <- ba[,shared]
  pilot <- pilot[,shared]
  
  # Concatenate
  data <- t(rbind(
    ba,
    pilot
  ))
  
  # Trim meta-data to dataset patients, reorder patients
  meta <- meta[intersect(row.names(meta), colnames(data)),]
  data <- data[,row.names(meta)]
  
  # Define race as a covariate to preserve
  model <- model.matrix(
    ~as.factor(Race),
    meta
  )
  
  # Correct using PNNL's missing-friendly ComBat
  corrected <- ComBat.NA(
    data,
    as.factor(meta$Study),
    par.prior = TRUE,
    mod = model
  )
  corrected <- corrected$`corrected data`
  
  # Separate Pilot from BeatAML samples
  ba_corrected <- data.frame(
    corrected[,row.names(meta[meta$Study == "BeatAML",])]
  )
  pilot_corrected <- data.frame(
    corrected[,row.names(meta[meta$Study == "pilotStudy",])]
  )
  
  return(list(ba = ba_corrected, pilot = pilot_corrected))
}
