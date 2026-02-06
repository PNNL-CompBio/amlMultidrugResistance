library(EnhancedVolcano)
library(here)
library(KSEAapp)

here::i_am("src/r/proteomics/run_ksea.R")

#' KSEA for phosphoproteomic kinase-enrichment analysis
#'
#' Uses KSEA to evaluate phosphoproteomic measurements for concerted trends
#' in kinase function.
#' @param phospho DataFrame with phosphoproteomic measurements.
#' @param meta DataFrame with patient meta-data, including race.
#' @return [list] with the following elements:
#' \itemize{
#'  \item table: DataFrame of KSEA scores, LFC, and p-values
#'  \item links: DataFrame of found Kinase-Substrate links
#' }
run_ksea <- function(phospho, meta) {
  # Load KSEA Kinase-Substrate link data, if needed
  if (!file.exists(here("src/r/proteomics/ks_data.csv"))) {
    download.file(
      "https://raw.githubusercontent.com/casecpb/KSEA/refs/heads/master/
      PSP%26NetworKIN_Kinase_Substrate_Dataset_July2016.csv",
      here("src/r/proteomics/ks_data.csv")
    )
  }
  ks_data <- read.csv(here("src/r/proteomics/ks_data.csv"))
  
  # Import LIMMA function
  source(here("src/r/proteomics/limma_correlations.R"))
  
  # Trim meta to patients with phosphoproteomic measurements
  meta <- meta[
    rownames(phospho),
  ]
  
  # Run LIMMA, get LFC
  limma_res <- run_limma(phospho, meta)$results
  
  # Adjust names with additional periods
  period_counts <- gregexpr("\\.", rownames(limma_res))
  period_counts <- lengths(regmatches(rownames(limma_res), period_counts))
  rownames(limma_res)[period_counts == 2] <- sub(
    "\\.",
    "-",
    rownames(limma_res)[period_counts == 2]
  )
  
  # Isolate gene and residue names
  phospho_names <- strsplit(rownames(limma_res), "\\.")
  
  # Swap residue format to match KSEA format
  residue <- unlist(phospho_names)[c(FALSE, TRUE)]
  residue <- gsub("[sty]", ";", residue)
  residue <- substr(residue, start=1, stop=nchar(residue) - 1)
  
  # Build KSEA table
  ksea_table <- data.frame(
    Protein = "NULL",
    Gene = unlist(phospho_names)[c(TRUE, FALSE)],
    Peptide = "NULL",
    Residue.Both = residue,
    p = limma_res$adj.P.Val,
    FC = 2 ** limma_res$logFC
  )

  # Run KSEA, get enrichment scores and kinase-substrate links
  result <- KSEA.Scores(
    ks_data, 
    ksea_table, 
    NetworKIN = FALSE
  )
  links <- KSEA.KS_table(
    ks_data, 
    ksea_table, 
    NetworKIN = FALSE
  )
  
  # Plots volcano
  volcano <- EnhancedVolcano(
    result,
    lab = result$Kinase.Gene,
    x = 'z.score',
    y = 'FDR',
    pCutoff = 5e-02,
    xlab = "KSEA Z-Score",
    ylab = "FDR q-value",
    legendLabels = NULL,
    selectLab = c("PRKACA")
  )
  
  return(list(table = result, links = links, volcano = volcano))
}
