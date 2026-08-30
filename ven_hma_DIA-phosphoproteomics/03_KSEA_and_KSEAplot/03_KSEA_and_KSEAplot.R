library(dplyr)
library(ggplot2)
library(openxlsx)
library(synapser)

## synapse login with .Renviron
synLogin()

## define folder ID and the target files
folder_id <-"syn68733653"

file1_id <- synFindEntityId("PTRC_EXP28_Phospho_Stats_Results.xlsx", parent = folder_id)
file1_entity <- synGet(file1_id)
df <- read.xlsx(file1_entity$path)

file2_id <- synFindEntityId("PTRC_EXP28_KSEA_Dataset_July2016 1.csv", parent = folder_id)
file2_entity <- synGet(file2_id)
KSDB <- read.csv(file2_entity$path)

## parse out 2 comparisons from Limma output file
NRRef <- df %>% filter(contrast == "Norelapse-Refractory")   
RRef  <- df %>% filter(contrast == "Relapse-Refractory")  

## subject NRRef to following codes for KSEA analysis and plotting.
## prepare input for KSEA
fold_change <- NRRef$logFC
fold_change <- 2**fold_change

PhosInp <- data.frame(Protein = "NULL", Gene = NRRef$feature, Peptide = "NULL", 
                      Residue.Both = NRRef$feature, p = "NULL", FC = fold_change) %>%
  dplyr::mutate(Residue.Both = sub("^.*-", "", Residue.Both)) %>%
  dplyr::mutate(Gene = sub("^(.*)-[^-]*$", "\\1", Gene))

ksea_res_full <- KSEAapp::KSEA.Scores(KSDB, PhosInp, NetworKIN = TRUE, 
                                      NetworKIN.cutoff = Inf) 

ksea_res <- ksea_res_full %>%
  dplyr::select(Kinase.Gene, m, p.value, FDR, z.score) %>%
  dplyr::rename(kinase = Kinase.Gene, z_score = z.score, p_value = p.value,
                adj_p_val = FDR, site_size = m)


# filter out kinases enriched by >=5 phosphosites and have p-value <0.05
kinase <- ksea_res %>% filter(site_size >= 5 & p_value < 0.05)

kinase <- kinase %>%
  arrange(desc(z_score)) %>%
  mutate(kinase=factor(kinase, levels=kinase))


###ADDED BY SARA
##now get complete scores
ksea_comp <- KSEAapp::KSEA.Complete(KSDB, PhosInp, NetworKIN = TRUE, NetworKIN.cutoff = Inf, m.cutoff = 5, p.cutoff = 0.05) 

links <- readr::read_csv('Kinase-Substrate Links.csv')  |>
  dplyr::rename(kinase = 'Kinase.Gene') |>
  right_join(kinase)

links <- links |> 
  rowwise() |>
  mutate(site = paste0(c(`Substrate.Gene`, `Substrate.Mod`), collapse = '-'))


links1 <- links |>
  mutate(comparison = 'Norelapse-Refractory')

##we can also look at one kinase of interest

links |> subset(kinase == 'AURKA') |>
  ggplot(aes(x=reorder(site, log2FC), y = log2FC, fill = p_value)) + geom_bar(stat='identity') +
  coord_flip()
ggsave('nr_ref_aurka.png',height=9)

links$site[abs(links$log2FC) < 1.5] <- ""
 
ggplot(links, aes(x = reorder(kinase,z_score), y = log2FC, col = log2FC)) + 
  geom_boxplot(outliers=FALSE) + 
  geom_jitter() + 
  ggrepel::geom_label_repel(aes(label = site)) +
  coord_flip()

ggsave('nr_ref_subs.png', height=9)

####end add


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

ggsave('nr_ref_kins.png',height=9)


## subject RRef to following codes for KSEA analysis and plotting.
## prepare input for KSEA
fold_change <- RRef$logFC
fold_change <- 2**fold_change

PhosInp <- data.frame(Protein = "NULL", Gene = RRef$feature, Peptide = "NULL", 
                      Residue.Both = RRef$feature, p = "NULL", FC = fold_change) %>%
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


###ADDED BY SARA
##now get complete scores
ksea_comp <- KSEAapp::KSEA.Complete(KSDB, PhosInp, NetworKIN = TRUE, NetworKIN.cutoff = Inf, m.cutoff = 5, p.cutoff = 0.05) 

links <- readr::read_csv('Kinase-Substrate Links.csv')  |>
  dplyr::rename(kinase = 'Kinase.Gene') |>
  right_join(kinase)

links <- links |> 
  rowwise() |>
  mutate(site = paste0(c(`Substrate.Gene`, `Substrate.Mod`), collapse = '-')) 

links2 <- links |>
  mutate(comparison = 'Relapse-Refractory')

links |> 
  subset(kinase == 'AURKB') |>
  ggplot(aes(x = reorder(site, log2FC), y = log2FC, fill = p_value)) + 
  geom_bar(stat = 'identity') +
  coord_flip()

ggsave('nr_ref_aurkb.png',height=9)

rbind(links1, links2) |>
  subset(kinase %in% c('AURKB','AURKA')) |>
  ggplot(aes(x=reorder(site, log2FC), y = log2FC, col = kinase, shape = comparison)) +
  geom_jitter() +
  coord_flip()

ggsave('aurk_test.png')

links$site[abs(links$log2FC) < 1.5] <- ""

ggplot(links, aes(x = reorder(kinase,z_score), y = log2FC, col = log2FC)) + 
  geom_boxplot(outliers=FALSE) + 
  geom_jitter() + 
  ggrepel::geom_label_repel(aes(label = site)) +
  coord_flip()

##

####end add
ggsave('rel_ref_subs.png',height=9)


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

ggsave('rel_ref_kins.png',height=8)
