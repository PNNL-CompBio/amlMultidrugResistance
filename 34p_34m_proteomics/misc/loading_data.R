library(PlexedPiper)
library(dplyr)
library(PNNL.DMS.utils)
library(Biostrings)
library(tidyverse)
library(MSnSet.utils)
library(cowplot)


boxplot_colors <- c("none 0h" = "grey",
                    "DMSO 0.5h" = "#A4FF9C",
                    "DMSO 6h" = "#5C9058",
                    "Ven 0.5h" = "#7654FF",
                    "Ven 6h" = "#4B35A1")
subtype_colors = c("Subtype 1" = "#d85600",
                   "Subtype 2" = "#7671b1",
                   "Subtype 3" = "#5b9902",
                   "Subtype 4" = "#f2a300",
                   "34+" = "firebrick",
                   "34-" = "steelblue")


## Load Braun data
m_corrected_p <- readRDS(syn$get("syn72330920")$path)
m_corrected_g <- readRDS(syn$get("syn72330921")$path)
m_corrected_a <- readRDS(syn$get("syn72330922")$path)
contrasts <- c("Ven_0.5h_34m-DMSO_0.5h_34m", "Ven_6h_34m-DMSO_6h_34m",
               "Ven_0.5h_34p-DMSO_0.5h_34p", "Ven_6h_34p-DMSO_6h_34p",
               "none_0h_34m-none_0h_34p")

## Load BeatAML data
m_beatAML = readRDS(syn$get("syn72644178")$path)
beataml_meta1 = read.table(syn$get("syn72644194")$path, sep = "\t", header = T)
beataml_meta2 = read.table(syn$get("syn26534982")$path, sep = "\t", header = T)
beataml_meta = left_join(beataml_meta1, beataml_meta2, by = c("Sample" = "SampleID.abbrev"))
drug_data <- read.csv(syn$get("syn51674470")$path)
drug_data_df <- rbind(drug_data %>% filter(inhibitor == "Venetoclax")) %>%
   select(Barcode.ID = sample_id, inhibitor, auc)