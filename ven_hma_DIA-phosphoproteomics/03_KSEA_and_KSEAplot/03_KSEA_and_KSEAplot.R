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

plot_extra = FALSE ##adding this to remove a lot of the extra stuff that Le added and didn't use

## parse out 2 comparisons from Limma output file
NRRef <- df %>% filter(contrast == "Norelapse-Refractory")   
RRef  <- df %>% filter(contrast == "Relapse-Refractory")  

## subject NRRef to following codes for KSEA analysis and plotting.
## prepare input for KSEA
#fold_change <- df$logFC#NRRef$logFC
#fold_change <- 2**fold_change

PhosInp <- data.frame(Protein = "NULL", Gene = NRRef$feature, Peptide = "NULL", 
                      Residue.Both = NRRef$feature, p = "NULL", FC = 2**NRRef$logFC) %>%
  dplyr::mutate(Residue.Both = sub("^.*-", "", Residue.Both)) %>%
  dplyr::mutate(Gene = sub("^(.*)-[^-]*$", "\\1", Gene))

ksea_nrref <- KSEAapp::KSEA.Scores(KSDB, PhosInp, NetworKIN = TRUE, 
                                      NetworKIN.cutoff = 1) |>
  mutate(comparison = 'Norelapse-Refractory')

ksea_res_nrref <- ksea_nrref %>%
  dplyr::select(Kinase.Gene, m, p.value, FDR, z.score,comparison) %>%
  dplyr::rename(kinase = Kinase.Gene, z_score = z.score, p_value = p.value,
                adj_p_val = FDR, site_size = m)


# filter out kinases enriched by >=5 phosphosites and have p-value <0.05
kinase_nrref_sig <- ksea_res_nrref %>% filter(site_size >= 5 & adj_p_val < 0.05)

kinase_nrref_sig <- kinase_nrref_sig %>%
  arrange(desc(z_score)) %>%
  mutate(kinase=factor(kinase, levels=kinase))


##now get complete scores
links_nrref <- KSEAapp::KSEA.KS_table(KSDB, PhosInp, NetworKIN = TRUE, NetworKIN.cutoff = 1)|>
  dplyr::rename(kinase = 'Kinase.Gene') 


links_nrref <- links_nrref |> 
  rowwise() |>
  mutate(site = paste0(c(`Substrate.Gene`, `Substrate.Mod`), collapse = '-'))


links_nrref_sig <- links_nrref|>
  #left_join(dplyr::rename(df,site='feature'))|>
  right_join(kinase_nrref_sig)
#, m.cutoff = 5, p.cutoff = 0.05) 

  #KSEAapp::KSEA.Complete(KSDB, PhosInp, NetworKIN = TRUE, NetworKIN.cutoff = Inf, m.cutoff = 5, p.cutoff = 0.05) 

# links <- readr::read_csv('Kinase-Substrate Links.csv')  |>
#   dplyr::rename(kinase = 'Kinase.Gene') |>
#   right_join(kinase)




##we can also look at one kinase of interest
if(plot_extra){

  links_nrref_sig |> subset(kinase == 'AURKA') |>
    ggplot(aes(x=reorder(site, log2FC), y = log2FC, fill = p_value)) + geom_bar(stat='identity') +
    coord_flip()
  if(plot_extra){
    ggsave('nr_ref_aurka.png',height=9)
  }
  links_nrref_sig$site[abs(links_nrref_sig$log2FC) < 1.5] <- ""
  
  ggplot(links_nrref_sig, aes(x = reorder(kinase,z_score), y = log2FC, col = log2FC)) + 
    geom_boxplot(outliers=FALSE) + 
    geom_jitter() + 
    ggrepel::geom_label_repel(aes(label = site)) +
    coord_flip()

  ggsave('nr_ref_subs.png', height=9)

}
####end add


# plot KSEA results with lollipop plot
ggplot(kinase_nrref_sig, aes(x = z_score, y = kinase)) +
  geom_segment(aes(x = 0, xend = z_score, y = kinase, yend = kinase),
               size=1.5) +
  # points, size by m, color by FDR
  geom_point(aes(size = site_size, colour = adj_p_val)) +
  scale_size(range = c(7, 12), name = "phosphosite\nsubstrates") +
  scale_colour_viridis_c(option = "plasma", direction = -1, name = "adjusted p-value") +
  geom_vline(xintercept=0, color="black", size=1.5)+
  labs(x = "z-score", y = NULL, title = "") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_line(color="grey90", linewidth=0.3),
    panel.grid.minor = element_line(color="grey90", linewidth=0.3),
    axis.title.x=element_text(size=30, face="plain"),
    axis.text.y=element_text(size=27, face="plain", color="black"),
    axis.text.x=element_text(size=30, face="plain", color="black"),
    legend.title=element_text(size=22, face="plain"),
    legend.text=element_text(size=20)
  )

ggsave('fig_5a_nrel_ref_kins.pdf',width=10,height=8)


## subject RRef to following codes for KSEA analysis and plotting.
## prepare input for KSEA
#fold_change <- RRef$logFC
#fold_change <- 2**fold_change

PhosInp <- data.frame(Protein = "NULL", Gene = RRef$feature, Peptide = "NULL", 
                      Residue.Both = RRef$feature, p = "NULL", FC = 2**RRef$logFC) %>%
  dplyr::mutate(Residue.Both = sub("^.*-", "", Residue.Both)) %>%
  dplyr::mutate(Gene = sub("^(.*)-[^-]*$", "\\1", Gene))

ksea_res_rref <- KSEAapp::KSEA.Scores(KSDB, PhosInp, NetworKIN = TRUE, NetworKIN.cutoff = 1) |>
  mutate(comparison = 'Relapse-Refractory')

ksea_res_rref <- ksea_res_rref %>%
  dplyr::select(Kinase.Gene, m, p.value, FDR, z.score,comparison) %>%
  dplyr::rename(kinase = Kinase.Gene, z_score = z.score, p_value = p.value,
                adj_p_val = FDR, site_size = m)

# filter out kinases enriched by >=5 phosphosites and have p-value <0.05
kinase_rref_sig <- ksea_res_rref %>% filter(site_size >=5 & adj_p_val < 0.05)

kinase_rref_sig <- kinase_rref_sig %>%
  arrange(desc(z_score)) %>%
  mutate(kinase=factor(kinase, levels=kinase))


###ADDED BY SARA
##now get complete scores
#ksea_comp <- KSEAapp::KSEA.Complete(KSDB, PhosInp, NetworKIN = TRUE, NetworKIN.cutoff = Inf, m.cutoff = 5, p.cutoff = 0.05) 

# links <- readr::read_csv('Kinase-Substrate Links.csv')  |>
#   dplyr::rename(kinase = 'Kinase.Gene') |>
#   right_join(kinase)

links_rref <- KSEAapp::KSEA.KS_table(KSDB, PhosInp, NetworKIN = TRUE, NetworKIN.cutoff = 1)|>
  dplyr::rename(kinase = 'Kinase.Gene')

links_rref <- links_rref |> 
  rowwise() |>
  mutate(site = paste0(c(`Substrate.Gene`, `Substrate.Mod`), collapse = '-')) 

links_rref_sig <- links_rref |>
  right_join(kinase_rref_sig)


if(plot_extra){
  links_rref_sig |> 
    subset(kinase == 'AURKB') |>
    ggplot(aes(x = reorder(site, log2FC), y = log2FC, fill = p_value)) + 
    geom_bar(stat = 'identity') +
    coord_flip()

  ggsave('nr_ref_aurkb.png',height=9)

  rbind(links_nrref_sig, links_rref_sig) |>
    subset(kinase %in% c('AURKB','AURKA')) |>
    ggplot(aes(x=reorder(site, log2FC), y = log2FC, col = kinase, shape = comparison)) +
    geom_jitter() +
    coord_flip()

  ggsave('aurk_test.png')

  links_rref_sig$site[abs(links_rref_sig$log2FC) < 1.5] <- ""

  ggplot(links_rref_sig, aes(x = reorder(kinase,z_score), y = log2FC, col = log2FC)) + 
    geom_boxplot(outliers=FALSE) + 
    geom_jitter() + 
    ggrepel::geom_label_repel(aes(label = site)) +
    coord_flip()

  ##

  ####end add
  ggsave('rel_ref_subs.png',height=9)
}

# plot KSEA results with lollipop plot
ggplot(kinase_rref_sig, aes(x = z_score, y = kinase)) +
  geom_segment(aes(x = 0, xend = z_score, y = kinase, yend = kinase),
               size=1.5) +
  # points, size by m, color by FDR
  geom_point(aes(size = site_size, colour =adj_p_val)) +
  scale_size(range = c(7, 12), name = "phosphosite\nsubstrates") +
  scale_colour_viridis_c(option = "plasma", direction = -1, name = "adjusted p-value") +
  geom_vline(xintercept=0, color="black", size=1.5)+
  labs(x = "z-score", y = NULL, title = "") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_line(color="grey90", linewidth=0.3),
    panel.grid.minor = element_line(color="grey90", linewidth=0.3),
    axis.title.x=element_text(size=30, face="plain"),
    axis.text.y=element_text(size=27, face="plain", color="black"),
    axis.text.x=element_text(size=30, face="plain", color="black"),
    legend.title=element_text(size=22, face="plain"),
    legend.text=element_text(size=20) 
  )

ggsave('fig_5b_rel_ref_kins.pdf',width=10,height=8)

full_tab <- rbind(ksea_res_nrref,ksea_res_rref)|>
  arrange(p_value) #|>
#  dplyr::select(-c(Substrate.Gene,Substrate.Mod,log2FC,site_size,site))|>
#  distinct()

write.xlsx(full_tab,file='Table_S4_KSEA_results.xlsx')
synStore(File('Table_S4_KSEA_results.xlsx', parent='syn76730214'))

##now to plot panel C

##collect any features that are significant
sig_feats <- subset(df, contrast %in% c('Norelapse-Refractory','Relapse-Refractory')& P.Value < 0.05 & abs(logFC)>1)

merged <- rbind(merge(links_nrref, ksea_res_nrref), merge(links_rref,ksea_res_rref))

merged |>
  subset(kinase == 'CDK1') |>
  subset(site %in% sig_feats$feature) |>
  ggplot(aes(x=reorder(site,log2FC),y=log2FC,fill=comparison)) + geom_bar(stat='identity',position='dodge') +
  scale_fill_manual(values=c('slategrey','lightgrey'))+
  theme_minimal(base_size = 12) +
  coord_flip()+labs(y = "log2 fold change", x = 'CDK1 substrate', title = "")

ggsave('fig_5c_cdk1_subs.pdf',width =8,height=7)

sig_feats <- subset(df, contrast %in% c('Norelapse-Refractory','Relapse-Refractory')& P.Value < 0.05 & abs(logFC)>1)

merged |>
  subset(kinase == 'MAPK1') |>
  subset(site %in% sig_feats$feature) |>
  ggplot(aes(x=reorder(site,log2FC),y=log2FC,fill=comparison)) + geom_bar(stat='identity',position='dodge') +
  scale_fill_manual(values=c('slategrey','lightgrey'))+
  theme_minimal(base_size = 12) +
  coord_flip()+labs(y = "log2 fold change", x = 'MAPK1 substrate', title = "")

ggsave('fig_5d_mapk1_subs.pdf',width=8,height=7)

##can we hone in on substrates?

kdiffs <- merged |>
  dplyr::select(kinase,site,comparison,log2FC) |>
  distinct() |>
  subset(kinase == 'MAPK1') |> tidyr::pivot_wider(names_from='comparison', values_from=log2FC) |>
  mutate(diff = `Norelapse-Refractory` - `Relapse-Refractory`) |>
  arrange(desc(diff))
##now ot plot Panel E, F

###READ IN DAT FROM last script

file1_id <- synFindEntityId("PTRC_EXP28_InSilico_Cleaned_PhosphositeData.xlsx", parent = folder_id)
file1_entity <- synGet(file1_id)
ddf <- read.xlsx(file1_entity$path)

file2_id <- synFindEntityId("PTRC_metadata_Exp28_removesamples.xlsx", parent = folder_id)
file2_entity <- synGet(file2_id)
meta <- read.xlsx(file2_entity$path)

## build MSnSet
exp <- ddf %>% arrange(., SITE) %>% 
  tidyr::pivot_longer(cols = 2:ncol(ddf),names_to='SampleID',values_to='log2Expr')

#fData <- data.frame(rownames(exprs)) %>% tibble::column_to_rownames(var="rownames.exprs.")
#exprs <- exprs[ , paste(meta$Sample, sep = "")]
full_dat <- meta %>% mutate(SampleID = Sample) %>% 
  right_join(exp)#tibble::column_to_rownames(var="Sample")

full_dat$subcohort = factor(full_dat$subcohort,levels=c('Norelapse','Relapse','Refractory'))

library(ggpubr)


for(psite in c('HMGN5-T31','RBMX-S326','FOXK1-S257','PRKAG2-S143','MAP7-S209','MKI67-S2783','PEAR1-S953','PEAR1-S1029')){

pse <- subset(full_dat, SITE==psite)
mval <- max(pse$log2Expr,na.rm=T)

  stat.df <- subset(df, feature==psite)|> 
    dplyr::select(contrast,p=P.Value) |> 
    tidyr::separate(contrast, sep='-',into=c('group1','group2')) |>
    mutate(y.position = c(mval + 1.2, mval + 0.8, mval + .3))

  pse |>
    ggplot(aes(x=subcohort, y = log2Expr,col=subcohort)) + geom_boxplot() + geom_point() +
    scale_color_manual(values=c('slategrey','darkgrey','black'))+
    theme_minimal(base_size = 12) +
    stat_pvalue_manual(stat.df)+
  labs(y = paste("log2 expression of", psite), x = 'Patient cohort', title = "")

  ggsave(paste0('fig5_',psite,'expr.pdf'),width=6,height=4)

}
