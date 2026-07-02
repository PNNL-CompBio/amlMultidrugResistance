library(dplyr)
library(EnhancedVolcano)
library(forcats)
library(ggplot2)
library(here)
library(KSEAapp)
library(limma)
library(synapser)

here::i_am("src/r/ksea_mutations/run_ksea_mutations.R")
token <- readLines(here("auth_token.txt"), warn = FALSE)
synapser::synLogin(authToken = token)


MUTATION <- c("NPM1", "FLT3_ITD")

ks_data <- read.csv(
  synapser::synGet("syn73849653")$path
)

# Trim meta to patients with phosphoproteomic measurements
phospho <- read.csv(here("src/r/ksea_mutations/phospho.csv"), row.names = 1)
meta <- read.csv(here("src/r/ksea_mutations/meta.csv"), row.names = 1)
meta <- meta[meta$Race %in% c("Black", "White"),]

# Append both mutation column, relabel NPM1/FLT3_ITD for patients with
# mutations in both genes
mutant_meta <- meta[,c("Race", "Study")]
mutant_meta$Mutation <- "WT"
mutant_meta[
  meta$NPM1 == "Mutant",
  "Mutation"
] <- "NPM1"
mutant_meta[
  meta$FLT3_ITD == "Mutant",
  "Mutation"
] <- "FLT3"
mutant_meta[
  meta$NPM1 == "Mutant" & meta$FLT3_ITD == "Mutant",
  "Mutation"
] <- "Both"

mutant_meta[] <- lapply(mutant_meta, factor)
mutant_meta$Mutation <- interaction(mutant_meta$Mutation, mutant_meta$Race)

design <- model.matrix(
  ~0 + Mutation + Study,
  data = mutant_meta,
  contrasts.arg = lapply(
    mutant_meta,
    contrasts, 
    contrasts = FALSE
  )
)

# Trim to patients with meta-data, format to expression matrix
phospho <- phospho[row.names(design),]
phospho <- t(as.matrix(phospho))

# Fit LIMMA, get contrasts by race
data_fit <- lmFit(phospho, design)

contrast_list <- list(
  "Black_NPM1" = makeContrasts(
    MutationNPM1.Black - MutationWT.Black,
    levels = design
  ),
  "Black_FLT3" = makeContrasts(
    MutationFLT3.Black - MutationWT.Black,
    levels = design
  ),
  "Black_Both" = makeContrasts(
    MutationBoth.Black - MutationWT.Black,
    levels = design
  ),
  "White_NPM1" = makeContrasts(
    MutationNPM1.White - MutationWT.White,
    levels = design
  ),
  "White_FLT3" = makeContrasts(
    MutationFLT3.White - MutationWT.White,
    levels = design
  ),
  "White_Both" = makeContrasts(
    MutationBoth.White - MutationWT.White,
    levels = design
  )
)

table <- data.frame(
  matrix("", nrow = 4, ncol = 3)
)
rownames(table) <- c(
  "Black - Increased activity",
  "Black - Reduced activity",
  "White - Increased activity",
  "White - Reduced activity"
)
colnames(table) <- c("NPM1", "FLT3", "Both")

for (comparison in names(contrast_list)) {
  cont.matrix <- contrast_list[[comparison]]
  race_contrasts <- contrasts.fit(data_fit, cont.matrix)
  race_contrasts <- eBayes(race_contrasts)
  limma_res <- topTable(race_contrasts, number = Inf, adjust.method = "BH")
  
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
  # at least 5 measured substrates
  df <- result
  df <- df[df$FDR < 0.05,]
  df <- df[df$m > 5,]
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
  
  ggsave(
    sprintf("%s%s.png", here("src/r/ksea_mutations/results/"), comparison),
    bar_chart
  )
  write.csv(
    result,
    sprintf("%s%s.csv", here("src/r/ksea_mutations/results/"), comparison)
  )
  
  race <- substr(comparison, 1, 5)
  mutation <- substr(comparison, nchar(comparison) - 3, nchar(comparison))
  
  table[sprintf("%s - Increased activity", race), mutation] <- paste(
    df[df$FDR < 0.05 & df$z.score > 0, "Kinase.Gene"], 
    collapse=", "
  )
  table[sprintf("%s - Reduced activity", race), mutation] <- paste(
    df[df$FDR < 0.05 & df$z.score < 0, "Kinase.Gene"], 
    collapse=", "
  )
}

write.csv(table, here("src/r/ksea_mutations/results/overview.csv"))
