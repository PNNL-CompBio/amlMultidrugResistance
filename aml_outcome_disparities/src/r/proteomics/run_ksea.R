library(dplyr)
library(EnhancedVolcano)
library(forcats)
library(ggplot2)
library(here)
library(KSEAapp)
library(synapser)

here::i_am("src/r/proteomics/run_ksea.R")
token <- readLines(here("auth_token.txt"), warn = FALSE)
synapser::synLogin(authToken = token)


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
  # Load KSEA Kinase-Substrate link data from Synapse
  ks_data <- read.csv(
    synapser::synGet("syn73849653")$path
  )  
  
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
  
  # Plot bar chart of kinase enrichment scores
  # Extract table of kinase z-scores, trim to those with an FDR < 0.05 and
  # at least 10 measured substrates
  df <- result
  df <- df[df$FDR < 0.05,]
  df <- df[df$m > 10,]
  df <- df[order(df$z.score),]
  
  # Plot bar chart, sort by z-scores
  bar_chart <- df %>%
    mutate(Kinase = fct_reorder(Kinase.Gene, z.score)) %>%
    ggplot(
      aes(
        x = Kinase, 
        y = z.score
      )
    ) + geom_bar(
      stat = "identity"
    ) + coord_flip()
  bar_chart <- bar_chart + ylab("Kinase Z-Score")
  bar_chart <- bar_chart + theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16)
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
  
  return(
    list(
      table = result,
      links = links, 
      volcano = volcano,
      bar_chart = bar_chart
    )
  )
}
