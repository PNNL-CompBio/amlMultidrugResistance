library(dplyr)
library(MSnSet.utils)
library(ggplot2)
library(ggpubr)
library(patchwork)



subcohort_colors =  c("#e8991b", "#00798c", "#d1495b", "grey27")
subcohort_colors2 =  c("#e8991b", "#00798c", "#d1495b", "grey27")
names(subcohort_colors) = c("Response_no_relapse", "Refractory", "Relapse", "Paired_relapse_sample")
vital_colors = c('orange3', 'plum4', 'darkgrey')
names(vital_colors) = c("Alive", "Dead", "LTFU")

## Load exp28
mat_path = syn$get("syn71717828")$path
xx = read.csv2(mat_path, sep = "\t")
# xx_ = xx %>% mutate(N_missing = rowSums(xx[, c(7:47)] == "")) %>%
#    group_by(Genes) %>% mutate(min_ = min(N_missing), max_ = max(N_missing)) %>%
#    ungroup() %>% filter(N_missing == min_)
# yy = read.csv2("../Aggregation_report.pg_matrix.tsv", sep = "\t")

mat_df = xx[, c(3, 7:47)] %>% as.data.frame()  ## Take Gene name and samples
colnames(mat_df) = sub("^.*WorkDir1\\.", "", colnames(mat_df))
colnames(mat_df) = sub("_FAIMS_.*$", "", colnames(mat_df))
colnames(mat_df) = sub("Astral_", "", colnames(mat_df))
colnames(mat_df) = sub("_03Nov25_Monty_ES906_3258_50SPD.mzML", "", colnames(mat_df))
colnames(mat_df)
## Rolling up, since there are duplicate gene names. Taking log2 then mean to rollup.
mat_df = mat_df %>%
   tidyr::pivot_longer(-Genes, names_to = "sample_name", values_to = "value") %>%
   mutate(value = as.numeric(value),
          value = log2(value)) %>%
   group_by(Genes, sample_name) %>%
   mutate(value_rollup = mean(value)) %>% ungroup() %>%
   select(Genes, sample_name, value_rollup) %>% unique() %>%
   tidyr::pivot_wider(names_from = "sample_name", values_from = "value_rollup")
mat_df = mat_df %>% filter(Genes != "")
mat = mat_df[, -1] %>% as.matrix()
rownames(mat) = mat_df$Genes


meta_path = syn$get("syn69058992")$path
meta = openxlsx::read.xlsx(meta_path)
meta$survival_days = meta[[13]]
p_data = meta %>% mutate(sample_name = datasetNameGlobal,
                         subcohort = gsub(" - ", "_", `sub-cohort.PNNL`),
                         subcohort = gsub(" ", "_", subcohort),
                         group = subcohort,
                         group = factor(group, levels = c("Response_no_relapse", "Refractory", "Relapse", "Paired_relapse_sample")),
                         group_surv = case_when(subcohort == "Response_no_relapse" & survival_days >= 500 ~ "response_long",
                                                subcohort == "Response_no_relapse" & survival_days < 500 ~ "response_short",
                                                subcohort == "Refractory" ~ "Refractory",
                                                subcohort == "Relapse" & survival_days < 500 ~ "relapse_short",
                                                subcohort == "Relapse" & survival_days >= 500 ~ "relapse_long",
                                                TRUE ~ subcohort))
rownames(p_data) = p_data$sample_name
p_data = p_data[colnames(mat), ]


## This is the Ven + HMA dataset (41 samples)
m_exp28 = MSnSet(exprs = mat, pData = p_data)
# hist(apply(exprs(m_exp28), 2, mean, na.rm = T))
# hist(apply(exprs(m_exp28), 1, mean, na.rm = T))
exprs(m_exp28) = sweep(exprs(m_exp28), 2, apply(exprs(m_exp28), 2, mean, na.rm = T), FUN = '-')




## Load corrected exp25 dataset, subset to 210 features?

## This is the ex vivo dataset
m_4pat <- readRDS(syn$get("syn64605135")$path)
m_210 <- readRDS(syn$get("syn72644178")$path)
m_210 <- m_210[complete.cases(exprs(m_210)), ]
m_4pat <- m_4pat[featureNames(m_4pat) %in% featureNames(m_210), ]
drug_data <- read.csv(syn$get("syn51674470")$path)

# zz = read.csv("../210_cohort_drug_auc_syn51674470.csv")
# table(zz$inhibitor)

## Inhibitors of interest
drug_data_df <- rbind(drug_data %>% filter(inhibitor == "Venetoclax"),   ## 127 samples)
                      drug_data %>% filter(inhibitor == "Azacytidine - Venetoclax"),  ## 20 samples
                      drug_data %>% filter(inhibitor == "Bortezomib - Venetoclax"),  ## 19 samples
                      drug_data %>% filter(inhibitor == "Dasatinib - Venetoclax"),  ## 92 samples
                      drug_data %>% filter(inhibitor == "Doramapimod - Venetoclax"),  ## 110 samples
                      drug_data %>% filter(inhibitor == "GW-2580 - Venetoclax"),  ## 20 samples
                      drug_data %>% filter(inhibitor == "Idelalisib - Venetoclax"),  ## 107 samples
                      drug_data %>% filter(inhibitor == "Olaparib - Venetoclax"),  ## 18 samples
                      drug_data %>% filter(inhibitor == "Quizartinib - Venetoclax"),  ## 80 samples
                      drug_data %>% filter(inhibitor == "Ruxolitinib - Venetoclax"),  ## 111 samples
                      drug_data %>% filter(inhibitor == "Sorafenib - Venetoclax"),  ## 92 samples
                      drug_data %>% filter(inhibitor == "Trametinib - Venetoclax"),  ## 84 samples
                      drug_data %>% filter(inhibitor == "Venetoclax - Artemisinin"),  ## 76 samples
                      drug_data %>% filter(inhibitor == "Venetoclax - Ibrutinib"),  ## 119 samples
                      drug_data %>% filter(inhibitor == "Venetoclax - JQ1"),  ## 85 samples
                      drug_data %>% filter(inhibitor == "Venetoclax - Palbociclib"),  ## 85 samples
                      drug_data %>% filter(inhibitor == "Venetoclax - Panobinostat"),  ## 82 samples
                      drug_data %>% filter(inhibitor == "Venetoclax - PH797804")  ## 11 samples
)






