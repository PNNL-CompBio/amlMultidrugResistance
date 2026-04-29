library(dplyr)
library(ggplot2)
library(openxlsx)


# load Limma result
df <- read.xlsx("PTRC_EXP28_Refractory VS relapse and no-relapse.xlsx", check.names=FALSE)
## load kinase substrate database
KSDB <- read.csv('PSP&NetworKIN_Kinase_Substrate_Dataset_July2016 1.csv', stringsAsFactors = FALSE)

## parse out 3 comparisons from Limma output file
NRRef <- df %>% filter(contrast == "Norelapse-Refractory")   
RRef  <- df %>% filter(contrast == "Relapse-Refractory")  
NRR   <- df %>% filter(contrast == "Relapse-Norelapse")  

## subject NRRef, RRef, NRR to following codes for KSEA analysis and plotting.
## prepare input for KSEA
fold_change <- NRRef$logFC
fold_change <- 2**fold_change

PhosInp <- data.frame(Protein = "NULL", Gene = NRRef$feature, Peptide = "NULL", 
                      Residue.Both = NRRef$feature, p = "NULL", FC = fold_change) %>%
  dplyr::mutate(Residue.Both = sub("^.*-", "", Residue.Both)) %>%
  dplyr::mutate(Gene = sub("^(.*)-[^-]*$", "\\1", Gene))

ksea_res_full <- KSEAapp::KSEA.Scores(KSDB, PhosInp, NetworKIN = TRUE, NetworKIN.cutoff = Inf) 

ksea_res <- ksea_res_full %>%
  dplyr::select(Kinase.Gene, m, p.value, FDR, z.score) %>%
  dplyr::rename(kinase = Kinase.Gene, z_score = z.score, p_value = p.value,
                adj_p_val = FDR, site_size = m)

# filter out kinases enriched by >=5 phosphosites and have p-value <0.05
kinase <- ksea_res %>% filter(site_size >=5 & p_value < 0.05)

kinase <- kinase %>%
  arrange(desc(z_score)) %>%
  mutate(kinase=factor(kinase, levels=kinase))


# plot KSEA results with lollipop plot
ggplot(kinase, aes(x = z_score, y = kinase)) +
  geom_segment(aes(x = 0, xend = z_score, y = kinase, yend = kinase),
               size=1.5) +
  # points, size by m, color by FDR
  geom_point(aes(size = site_size, colour = p_value)) +
  scale_size(range = c(7, 12), name = "phosphosite\nsubstrates") +
  scale_colour_viridis_c(option = "plasma", direction = -1, name = "p-value") +
  geom_vline(xintercept=0, color="black", size=1.5)+
  labs(x = "z-score", y = NULL, title = "") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_line(color="grey90", linewidth=0.3),
    panel.grid.minor = element_line(color="grey90", linewidth=0.3),
    axis.title.x=element_text(size=30, face="bold"),
    axis.text.y=element_text(size=27, face="bold", color="black"),
    axis.text.x=element_text(size=30, face="bold", color="black"),
    legend.title=element_text(size=22, face="bold"),
    legend.text=element_text(size=20)
  )

