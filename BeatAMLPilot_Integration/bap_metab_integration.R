library(magrittr)
library(pmartR)
library(malbacR)
library(igraph)
library(here)
library(ggplot2)
library(plotly)
library(synapser)
synapser::synLogin()

here::i_am("BeatAMLPilot_Integration/bap_metab_integration.R")

# Here we explore integration of BeatAML and Pilot Metabolomics datasets via
# ComBat, which is an approach that has historically done best in this regard. 

# Also included are various statistical analyses based on different versions of the data,
# i.e. 1) Pilot data alone, 2) BeatAML data alone, and 3) Integrated data

# We also apply some pre-processing to these data as needed. 


### Data Import -------------------------------------------------------------

## BeatAML Data -----------------------------------------------------

#all data
all_ba_data <- readxl::read_xlsx(synapser::synGet('syn68710501')$path) 

#patient data without proteomics sample numbers
ba_pats <- all_ba_data |>
  dplyr::mutate(RNASeq=ifelse(is.na(dbgap_rnaseq_sample),FALSE,TRUE), 
                Mutation = ifelse(is.na(dbgap_dnaseq_sample),FALSE,TRUE)) |>
  dplyr::select(labId,consensus_sex,reportedRace,reportedEthnicity,inferred_ethnicity,RNASeq, Mutation) |>
  dplyr::distinct()


##now we have extra metadata that wasn't accounted for
ex_data <- readxl::read_xlsx(synGet('syn68731261')$path)|>
  dplyr::mutate(RNASeq=FALSE,Mutation=FALSE) |>
  dplyr::select(labId='Accession',consensus_sex='Sex',reportedRace,reportedEthnicity,
                inferred_ethnicity='inferredEthnicity',RNASeq, Mutation) |>
  dplyr::distinct()

ba_pats <- rbind(ba_pats, ex_data)

all_ba_meta <- readxl::read_xlsx(synapser::synGet('syn25796769')$path)|>
  dplyr::left_join(ba_pats,by = 'labId')|>
  dplyr::rename(CALR='Calreticulin', NPM1='NPM1calls')

ba_meta <- all_ba_meta |>
  #mutate(Accession = labId) |>
  tidyr::separate(`SampleID(full)`, into = c("PTRC",'exp10','SampleID'),sep = '_') |>
  mutate(SampleID=as.numeric(SampleID))|>
  mutate(Race = ifelse(is.na(reportedRace),inferred_ethnicity,reportedRace))|> ##focus on reported RACE first
  mutate(Race = stringr::str_remove(Race,'Admixed')) |> ##remove the admixed 
  dplyr::select(Accession='labId', SampleID,'FLT3-ITDcalls',Sex='gender',Age='ageAtDiagnosis',Race,RNASeq, Mutation)|>
  mutate(Proteomics=TRUE,Phospho=TRUE) |> ##all these patients have proteomics
  mutate(study='BeatAML', source='BeatAML') # Note (8/6/2025): defined based on guidance from SG

ba_metab <- readxl::read_xlsx(synapser::synGet('syn53678273')$path) |>
  tidyr::pivot_longer(cols = starts_with('BEAT_AML'), names_to = 'sample', values_to = 'value') |>
  tidyr::separate(sample, into = c('Beat','AML','PNL','SampleID.abbrev', 'm','hilic','posneg'), 
                  sep = '_') |>
  dplyr::mutate(SampleID = as.numeric(SampleID.abbrev)) |>
  dplyr::left_join(ba_meta, by = 'SampleID') |>
  dplyr::select('m/z','Name','Standardized name','Super class','Main class','SampleID',
                'Accession')

# BeatAML metabolomics is only HILIC Pos - there aren't other modes
# represented, e.g. HILIC Neg, RP Pos, RP Neg
ba_metab_hp <- readxl::read_xlsx(synapser::synGet('syn53678273')$path) %>%
  setNames(make.names(colnames(.))) %>%
  setNames(gsub("_HILIC_POS", "", colnames(.), ignore.case = TRUE)) %>%
  setNames(gsub("Blank_BEAT_AML_022", "Blank_BEAT_AML_02_2", colnames(.), ignore.case = TRUE))

# ba_metab2 <- readxl::read_xlsx(synapser::synGet('syn53678273')$path) |>
#   tidyr::pivot_longer(cols = starts_with('BEAT_AML'), names_to = 'sample', values_to = 'value') |>
#   tidyr::separate(sample, into = c('Beat','AML','PNL','SampleID.abbrev', 'm','hilic','posneg'), 
#                   sep = '_') |>
#   dplyr::mutate(SampleID = as.numeric(SampleID.abbrev)) |>
#   dplyr::left_join(ba_meta, by = 'SampleID') |>
#   dplyr::select('m/z','Name','Standardized name','Super class','Main class','SampleID',
#                 'Accession')

## Pilot Data -------------------------------------------------------

#new pilot metaadata
sample_exp_metadata <- readxl::read_xlsx(synGet('syn68522106')$path) |>
  dplyr::rename(CEBPA = 'CEBPA_BZIP',FLT3 = 'FLT3_TKD',Race='Race_label')|>
  mutate(Sex = ifelse(`Sex (1-male, 0-female)` == 1,'Male','Female'))|>
  tidyr::pivot_longer(cols=c(11:40),names_to='gene',values_to='mutation')|>
  mutate(mutation=ifelse(is.na(mutation),'Not measured',ifelse(mutation==0,'WT','Mutant')))

#lookup table for linking metadata
pi_metadat <- readxl::read_xlsx(synGet('syn68835814')$path) %>%
  dplyr::rename(Accession = `Accession/\r\nBarcode`) %>%
  setNames(colnames(.)) %>%
  dplyr::mutate(Accession_og = Accession) %>%
  dplyr::mutate(Accession = dplyr::case_when(
    Accession == "PS88-0050 (Ricky, P)" ~ "PS88-0050",
    Accession == "PS88-0140 (Yates, N)" ~ "PS88-0140",
    Accession == "PS89-0158 (Alston Marvin)" ~ "PS89-0158",
    TRUE ~ Accession
  )) %>%
  dplyr::mutate(SampleID = `Metabolomics ID`) %>%
  dplyr::left_join(sample_exp_metadata %>% 
                     dplyr::select(Accession, Age, Sex, Race) %>%
                     dplyr::distinct(), by = "Accession") %>%
  dplyr::mutate(SampleID = gsub("Exp26", "exp26", SampleID))

all(pi_metadat$Accession %in% sample_exp_metadata$Accession)

# Here, there are HILIC Pos, HILIC Neg, RP Pos, and RP Neg.
# However, we will only read in HILIC Pos since that is all
# that is available in the BeatAML cohort. 
mfile <- synGet('syn68710369')$path

pi_metab_hp <- readxl::read_xlsx(mfile,sheet = 'HILIC Positive') %>%
  setNames(make.names(colnames(.))) %>%
  setNames(gsub("_HILIC_POS", "", colnames(.), ignore.case = TRUE))

# -------------------------------------------------------------------------

### Data Processing ---------------------------------------------------------

## BeatAML Data -------------------------------------------------------------------------
# Duplicate Handling ------------------------------------------------------

# Check for duplicate metabolites and append _A, _B, _C, etc labels as needed
ba_metab_hp <- ba_metab_hp %>%
  dplyr::mutate(Name_og = Name, 
                Name = make.names(paste0("hipos_", Name))) %>%
  dplyr::group_by(Name) %>%
  dplyr::mutate(dupcount = dplyr::n(), # count total duplicates for a given Name
                dupnum = dplyr::row_number(), # assign numeric labels for each duplicate
                dupletter = toupper(letters[dupnum]), # convert numeric labels to capital letters
                # Rename metabolites by appending duplicate labels where appropriate
                Name = ifelse(dupcount > 1, paste0(Name, "__", dupletter), Name)) %>%
  dplyr::ungroup() %>%
  dplyr::relocate(Name, Name_og) %>%
  dplyr::relocate(dupcount, dupnum, dupletter, .after = `Name_og`)

# No duplicates found
ba_hp_dup_df <- ba_metab_hp %>% 
  dplyr::filter(dupcount > 1) %>% 
  dplyr::select(RT..min., m.z, Name_og, Reference.Ion, Name) %>% 
  dplyr::arrange(Name_og)

# pmartR formatting -------------------------------------------------------

# Create initial edata, emeta, and fdata objects
ba_metab_hp_edat <- ba_metab_hp %>%
  dplyr::select(Name, 
                dplyr::contains("BEAT_AML_PNL_"))
ba_metab_hp_emet <- ba_metab_hp %>%
  dplyr::select(-dplyr::contains("BEAT_AML_PNL_"))

ba_metab_hp_fdat <- data.frame(SampleID = colnames(ba_metab_hp_edat)[-1]) %>%
  dplyr::mutate(link_id = purrr::map_dbl(SampleID, function(x){
    as.numeric(strsplit(x, "_")[[1]][4]) # gets the sample number at the end of SampleID string
  })) %>%
  dplyr::left_join(ba_meta %>%
                     dplyr::rename(link_id = SampleID)) %>%
  dplyr::mutate(study = ifelse(is.na(study), "BeatAML", study)) %>%
  # For now, select subset of variables that are in common with pilot fdata
  dplyr::select(SampleID, Age, Sex, Race, study) %>%
  dplyr::mutate(Race_mod = dplyr::case_when(
    Race == "Unknown" ~ NA,
    Race == "NA" ~ NA,
    Race %in% c("Black", "HispNative", "Asian") ~ "Non-white",
    TRUE ~ Race
  ))

# Check consistency of edat and fdat names
all(ba_metab_hp_fdat$SampleID %in% colnames(ba_metab_hp_edat)[-1])

ba_pm_hp <- as.metabData(e_data = ba_metab_hp_edat,
                         f_data = ba_metab_hp_fdat,
                         e_meta = ba_metab_hp_emet,
                         edata_cname = "Name", 
                         fdata_cname = "SampleID",
                         emeta_cname = "Name",
                         data_scale = "abundance",
                         data_types = "HILIC Positive")

ba_pm_hp <- edata_transform(ba_pm_hp, data_scale = "log2")

# -------------------------------------------------------------------------

## Pilot Data -------------------------------------------------------------------------
# Duplicate Handling ------------------------------------------------------

# Check for duplicate metabolites and append _A, _B, _C, etc labels as needed
pi_metab_hp <- pi_metab_hp %>%
  dplyr::mutate(Name_og = Name, 
                Name = make.names(paste0("hipos_", Name))) %>%
  dplyr::group_by(Name) %>%
  dplyr::mutate(dupcount = dplyr::n(), # count total duplicates for a given Name
                dupnum = dplyr::row_number(), # assign numeric labels for each duplicate
                dupletter = toupper(letters[dupnum]), # convert numeric labels to capital letters
                # Rename metabolites by appending duplicate labels where appropriate
                Name = ifelse(dupcount > 1, paste0(Name, "__", dupletter), Name)) %>%
  dplyr::ungroup() %>%
  dplyr::relocate(Name, Name_og) %>%
  dplyr::relocate(dupcount, dupnum, dupletter, .after = `Name_og`)

# No duplicates
pi_hp_dup_df <- pi_metab_hp %>% 
  dplyr::filter(dupcount > 1) %>% 
  dplyr::select(RT..min., m.z, Name_og, Reference.Ion, Name) %>% 
  dplyr::arrange(Name_og)

# pmartR formatting -------------------------------------------------------

# Create initial edata, emeta, and fdata objects
pi_metab_hp_edat <- pi_metab_hp %>%
  dplyr::select(Name, 
                dplyr::contains("PTRC_exp26_Met")) %>%
  dplyr::select(-dplyr::contains("Blank")) %>%
  dplyr::select(-dplyr::contains("QC_Pool"))
pi_metab_hp_emet <- pi_metab_hp %>%
  dplyr::select(Name:PTRC_exp26_Met_QC_Pool_07)

pi_metab_hp_fdat <- data.frame(SampleID = colnames(pi_metab_hp_edat)[-1]) %>%
  dplyr::mutate(link_id = purrr::map_dbl(SampleID, function(x){
    as.numeric(strsplit(x, "_")[[1]][4]) # gets the sample number at the end of SampleID string
  })) %>%
  dplyr::mutate(study = "pilot") %>%
  dplyr::left_join(pi_metadat) %>%
  # For now, select subset of variables that are in common with BeatAML fdata
  dplyr::select(SampleID, Age, Sex, Race, study) %>%
  dplyr::mutate(Race_mod = dplyr::case_when(
    Race == "Unknown" ~ NA,
    Race == "NA" ~ NA,
    Race %in% c("Black", "HispNative", "Asian") ~ "Non-white",
    TRUE ~ Race
  ))

# Check consistency of edat and fdat names
all(pi_metab_hp_fdat$SampleID %in% colnames(pi_metab_hp_edat)[-1])

pi_pm_hp <- as.metabData(e_data = pi_metab_hp_edat,
                         f_data = pi_metab_hp_fdat,
                         e_meta = pi_metab_hp_emet,
                         edata_cname = "Name", 
                         fdata_cname = "SampleID",
                         emeta_cname = "Name",
                         data_scale = "abundance",
                         data_types = "HILIC Positive")

pi_pm_hp <- edata_transform(pi_pm_hp, data_scale = "log2")

# -------------------------------------------------------------------------


### Flagging Common Metabolites ---------------------------------------------

# Flag based on:
# 1) shared name
# 2) m.z.
# 3) R.T.
# However, this matching may need to be done manually, as is the case for
# lipids.

temp_ba <- ba_pm_hp$e_meta %>%
  dplyr::select(Name, Name_og, 
                Standardized.name,
                m.z, RT..min.) %>%
  dplyr::mutate(cohort = 1) %>%
  dplyr::mutate(round_RT = round(RT..min., 2),
                round_mz = round(m.z, 2)) %>%
  dplyr::mutate(match_Name = "", 
                match_Name_og = "",
                match_Standardized.name = "",
                match_m.z = NA,
                match_RT..min. = NA) %>%
  dplyr::relocate(match_Name, .after = Name) %>%
  dplyr::relocate(match_Name_og, .after = Name_og) %>%
  dplyr::relocate(match_Standardized.name, .after = Standardized.name) %>%
  dplyr::relocate(match_m.z, .after = m.z) %>%
  dplyr::relocate(match_RT..min., .after = RT..min.)

temp_pi <- pi_pm_hp$e_meta %>%
  dplyr::select(Name, Name_og, 
                Standardized.name,
                m.z, RT..min.) %>%
  dplyr::mutate(cohort = 2) %>%
  dplyr::mutate(round_RT = round(RT..min., 2),
                round_mz = round(m.z, 2))

for(i in 1:nrow(temp_ba)){
  
  # Compute the euclidean distance between the mz/rt values of
  # each lipid sharing the same name as the query, and then take 
  # the top match
  cand_pi_matches <- temp_pi %>%
    dplyr::filter(Name_og == temp_ba$Name_og[i]) %>%
    dplyr::mutate(pt_dist = sqrt((RT..min.-temp_ba$RT..min.[i])^2 + 
                                   (m.z-temp_ba$m.z[i])^2)) %>%
    dplyr::arrange(pt_dist)
  
  if(nrow(cand_pi_matches) == 0){
    next
  } else{
    
    cand_pi_matches <- cand_pi_matches %>%
      dplyr::slice_head(n = 1)
  }
  
  temp_ba$match_Name[i] <- cand_pi_matches$Name
  temp_ba$match_Name_og[i] <- cand_pi_matches$Name_og
  temp_ba$match_Standardized.name[i] <- cand_pi_matches$Standardized.name
  temp_ba$match_m.z[i] <- cand_pi_matches$m.z
  temp_ba$match_RT..min.[i] <- cand_pi_matches$RT..min.
}

bapi_matchkey_hp <- temp_ba %>%
  dplyr::filter(match_Name != "")

rm(cand_pi_matches, temp_ba, temp_pi)

# -------------------------------------------------------------------------


### Normalization + ComBat Only ----------------------------------------------------------

## Normalization -----------------------------------------------------------
# BeatAML -------------------------------------------------------------------

# Normalize data using median centering
ba_pm_hp_norm <- normalize_global(omicsData = ba_pm_hp,
                                  subset_fn = "all",
                                  norm_fn = "median",
                                  apply_norm = TRUE,
                                  backtransform = TRUE)

# Pilot -------------------------------------------------------------------

# Normalize data using median centering
pi_pm_hp_norm <- normalize_global(omicsData = pi_pm_hp,
                                  subset_fn = "all",
                                  norm_fn = "median",
                                  apply_norm = TRUE,
                                  backtransform = TRUE)


## Integrate ---------------------------------------------------

# Define subsets containing only common lipids

# BeatAML
my_filt <- custom_filter(ba_pm_hp_norm, 
                         e_data_keep = bapi_matchkey_hp$Name)
ba_pm_hp_norm <- applyFilt(filter_object = my_filt,
                           omicsData = ba_pm_hp_norm)


# Pilot
my_filt <- custom_filter(pi_pm_hp_norm, 
                         e_data_keep = bapi_matchkey_hp$match_Name)
pi_pm_hp_norm <- applyFilt(filter_object = my_filt,
                           omicsData = pi_pm_hp_norm)


# Rename the formatted_name in Pilot to match that of BeatAML
for(i in 1:nrow(bapi_matchkey_hp)){
  tempname <- bapi_matchkey_hp$Name[i]
  
  temp_repname <- bapi_matchkey_hp$match_Name[i]
  repidx <- which(pi_pm_hp_norm$e_meta$Name == temp_repname)
  
  pi_pm_hp_norm$e_meta$Name[repidx] <- tempname
  pi_pm_hp_norm$e_data$Name[repidx] <- tempname
}
all(pi_pm_hp_norm$e_data$Name %in% ba_pm_hp_norm$e_data$Name) &
  all(ba_pm_hp_norm$e_data$Name %in% pi_pm_hp_norm$e_data$Name)


# Define new pmartR objects based on a concatenation of e_data, f_data
# Neg Mode
edat_common_hp <- ba_pm_hp_norm$e_data %>%
  dplyr::left_join(pi_pm_hp_norm$e_data, by = "Name")

fdat_common_hp <- rbind.data.frame(ba_pm_hp_norm$f_data %>%
                                      dplyr::mutate(batchid = 1),
                                   pi_pm_hp_norm$f_data %>%
                                      dplyr::mutate(batchid = 2))

emeta_common_hp <- ba_pm_hp_norm$e_meta %>%
  dplyr::left_join(pi_pm_hp_norm$e_meta, by = "Name") %>%
  setNames(gsub("\\.x", "\\.beataml", colnames(.))) %>%
  setNames(gsub("\\.y", "\\.pilot", colnames(.)))

# all(names(edat_common_hp)[-1] %in% fdat_common_hp$SampleID)
# setdiff(names(edat_common_hp)[-1], fdat_common_hp$SampleID)

pm_common_hp <- as.lipidData(e_data = edat_common_hp, 
                             f_data = fdat_common_hp,
                             e_meta = emeta_common_hp,
                             emeta_cname = "Name",
                             edata_cname = "Name", 
                             fdata_cname = "SampleID",
                             data_scale = "log2",
                             data_types = "HILIC Pos",
                             is_normalized = TRUE)

# Several samples (16 specifically: 4 Beat AML, 12 Pilot) have been removed due to
# missing race information
# pm_common_hp <- group_designation(pm_common_hp, main_effects = c("Race_mod"), batch_id = "batchid")
pm_common_hp <- group_designation(pm_common_hp, main_effects = c("study"), batch_id = "batchid")

pm_common_hp_combat <- bc_combat(omicsData = pm_common_hp, use_groups = FALSE)

## RMD Outlier Assessment --------------------------------------------------

# BeatAML -----------------------------------------------------------------

ba_pm_hp_norm <- group_designation(ba_pm_hp_norm, main_effects = c("study"))

myrmd <- rmd_filter(ba_pm_hp_norm)
png(here("BeatAMLPilot_Integration", "figures", "outlier_beataml_hp_rmd.png"), units="in",
    width=10, height=7, res=300)
plot(myrmd, pvalue_threshold = 0.00001)
dev.off()

# Pilot -------------------------------------------------------------------

pi_pm_hp_norm <- group_designation(pi_pm_hp_norm, main_effects = c("study"))

myrmd <- rmd_filter(pi_pm_hp_norm)
png(here("BeatAMLPilot_Integration", "figures", "outlier_pilot_hp_rmd.png"), units="in",
    width=10, height=7, res=300)
plot(myrmd, pvalue_threshold = 0.00001)
dev.off()

# Combined, No BC ---------------------------------------------------------

pm_common_hp <- group_designation(pm_common_hp, main_effects = c("study"), batch_id = "batchid")

myrmd <- rmd_filter(pm_common_hp)
png(here("BeatAMLPilot_Integration", "figures", "outlier_nobc_hp_rmd.png"), units="in",
    width=10, height=7, res=300)
plot(myrmd, pvalue_threshold = 0.00001)
dev.off()

# Combined, BC ------------------------------------------------------------

pm_common_hp_combat <- group_designation(pm_common_hp_combat, main_effects = c("study"), batch_id = "batchid")

myrmd <- rmd_filter(pm_common_hp_combat)
png(here("BeatAMLPilot_Integration", "figures", "outlier_combat_hp_rmd.png"), units="in",
    width=10, height=7, res=300)
plot(myrmd, pvalue_threshold = 0.00001)
dev.off()

# -------------------------------------------------------------------------

## PCA --------------------------------------------------------------------

# BeatAML -----------------------------------------------------------------

ba_pm_hp_norm <- group_designation(ba_pm_hp_norm, main_effects = c("study"))

mypca <- dim_reduction(ba_pm_hp_norm)
png(here("BeatAMLPilot_Integration", "figures", "beataml_hp_pca.png"), units="in",
    width=10, height=7, res=300)
plot(mypca, pvalue_threshold = 0.00001, title_lab = "HP: BeatAML")
dev.off()

# Pilot -------------------------------------------------------------------

pi_pm_hp_norm <- group_designation(pi_pm_hp_norm, main_effects = c("study"))

mypca <- dim_reduction(pi_pm_hp_norm)
png(here("BeatAMLPilot_Integration", "figures", "pilot_hp_pca.png"), units="in",
    width=10, height=7, res=300)
plot(mypca, pvalue_threshold = 0.00001, title_lab = "HP: Pilot")
dev.off()

# Combined, No BC ---------------------------------------------------------

mypca <- dim_reduction(pm_common_hp)

png(here("BeatAMLPilot_Integration", "figures", "nobc_hp_pca.png"), units="in",
    width=10, height=7, res=300)
plot(mypca, omicsData = pm_common_hp, color_by = "study", 
     title_lab = "HP: No Batch-Correction")+theme(text=element_text(size=21))
dev.off()

# Combined, BC ------------------------------------------------------------

mypca <- dim_reduction(pm_common_hp_combat)

png(here("BeatAMLPilot_Integration", "figures", "combat_hp_pca.png"), units="in",
    width=10, height=7, res=300)
plot(mypca, omicsData = pm_common_hp_combat, color_by = "study", 
     title_lab = "HP: ComBat")+theme(text=element_text(size=21))
dev.off()

# -------------------------------------------------------------------------

## ANOVA -------------------------------------------------------------------

# BeatAML -----------------------------------------------------------------

ba_pm_hp_norm_temp <- group_designation(ba_pm_hp_norm, main_effects = c("Race_mod"))

ba_hptest_race <- imd_anova(omicsData = ba_pm_hp_norm_temp, 
                            test_method = "anova", 
                            pval_adjust_a_multcomp = "holm", 
                            pval_thresh = 0.05)
summary(ba_hptest_race)

# ba_pm_hp_norm_temp <- group_designation(ba_pm_hp_norm, main_effects = c("study"))
# 
# ba_hptest_study <- imd_anova(omicsData = ba_pm_hp_norm_temp, 
#                              test_method = "anova", 
#                              pval_adjust_a_multcomp = "holm", 
#                              pval_thresh = 0.05)
# summary(ba_hptest_study)

# Pilot -------------------------------------------------------------------

# Not enough data in each group for a comparison
# pi_pm_hp_norm_temp <- group_designation(pi_pm_hp_norm, main_effects = c("Race_mod"))
# 
# pi_hptest_race <- imd_anova(omicsData = pi_pm_hp_norm_temp, 
#                             test_method = "anova", 
#                             pval_adjust_a_multcomp = "holm", 
#                             pval_thresh = 0.05)
# summary(pi_hptest_race)

# pi_pm_hp_norm_temp <- group_designation(pi_pm_hp_norm, main_effects = c("study"))
# 
# pi_hptest_study <- imd_anova(omicsData = pi_pm_hp_norm_temp, 
#                              test_method = "anova", 
#                              pval_adjust_a_multcomp = "holm", 
#                              pval_thresh = 0.05)
# summary(pi_hptest_study)

# Combined, No BC ---------------------------------------------------------

pm_common_hp_temp <- group_designation(pm_common_hp, main_effects = c("Race_mod"))

nobc_hptest_race <- imd_anova(omicsData = pm_common_hp_temp,
                              test_method = "anova",
                              pval_adjust_a_multcomp = "holm",
                              pval_thresh = 0.05)
summary(nobc_hptest_race)

pm_common_hp_temp <- group_designation(pm_common_hp, main_effects = c("study"))

nobc_hptest_study <- imd_anova(omicsData = pm_common_hp_temp,
                               test_method = "anova",
                               pval_adjust_a_multcomp = "holm",
                               pval_thresh = 0.05)
summary(nobc_hptest_study)

# Combined, BC ------------------------------------------------------------

pm_common_hp_combat_temp <- group_designation(pm_common_hp_combat, main_effects = c("Race_mod"))

combat_hptest_race <- imd_anova(omicsData = pm_common_hp_combat_temp,
                                test_method = "anova",
                                pval_adjust_a_multcomp = "holm",
                                pval_thresh = 0.05)
summary(combat_hptest_race)

pm_common_hp_combat_temp <- group_designation(pm_common_hp_combat, main_effects = c("study"))

combat_hptest_study <- imd_anova(omicsData = pm_common_hp_combat_temp,
                                 test_method = "anova",
                                 pval_adjust_a_multcomp = "holm",
                                 pval_thresh = 0.05)
summary(combat_hptest_study)

# -------------------------------------------------------------------------

