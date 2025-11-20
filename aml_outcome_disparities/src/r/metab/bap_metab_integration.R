library(magrittr)
library(pmartR)
library(malbacR)
library(igraph)
library(here)
library(ggplot2)
library(plotly)
library(synapser)
synapser::synLogin()

here::i_am("bap_metab_integration.R")

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

ba_metab_hp <- readxl::read_xlsx(synapser::synGet('syn53678273')$path, sheet = "HILIC Positive") %>%
  setNames(make.names(colnames(.))) %>%
  setNames(gsub("_HILIC_POS", "", colnames(.), ignore.case = TRUE)) %>%
  setNames(gsub("Blank_BEAT_AML_022", "Blank_BEAT_AML_02_2", colnames(.), ignore.case = TRUE))

ba_metab_hn <- readxl::read_xlsx(synapser::synGet('syn53678273')$path, sheet = "HILIC Negative") %>%
setNames(make.names(colnames(.))) %>%
  setNames(gsub("_HILIC_NEG", "", colnames(.), ignore.case = TRUE)) %>%
  setNames(gsub("Blank_BEAT_AML_012", "Blank_BEAT_AML_01_2", colnames(.), ignore.case = TRUE)) %>%
  setNames(gsub("Blank_BEAT_AML_022", "Blank_BEAT_AML_02_2", colnames(.), ignore.case = TRUE)) %>%
  setNames(gsub("Blank_BEAT_AML_023", "Blank_BEAT_AML_02_3", colnames(.), ignore.case = TRUE))

ba_metab_rp <- readxl::read_xlsx(synapser::synGet('syn53678273')$path, sheet = "RP Positive") %>%
  setNames(make.names(colnames(.))) %>%
  setNames(gsub("_RP_POS", "", colnames(.), ignore.case = TRUE)) %>%
  setNames(gsub("Blank_BEAT_AML_013", "Blank_BEAT_AML_01_3", colnames(.), ignore.case = TRUE)) %>%
  setNames(gsub("Blank_BEAT_AML_012", "Blank_BEAT_AML_01_2", colnames(.), ignore.case = TRUE))

ba_metab_rn <- readxl::read_xlsx(synapser::synGet('syn53678273')$path, sheet = "RP Negative") %>%
  setNames(make.names(colnames(.))) %>%
  setNames(gsub("_RP_Neg", "", colnames(.), ignore.case = TRUE)) %>%
  setNames(gsub("Blank_BEAT_AML_02_", "Blank_BEAT_AML_02_2", colnames(.), ignore.case = TRUE))

## Pilot Data -------------------------------------------------------

#new pilot metaadata
sample_exp_metadata <- readxl::read_xlsx(synGet('syn68522106')$path) |>
  dplyr::rename(CEBPA = 'CEBPA_BZIP',FLT3 = 'FLT3_TKD',Race='Race_label')|>
  mutate(Sex = ifelse(`Sex (1-male, 0-female)` == 1,'Male','Female'))|>
  tidyr::pivot_longer(cols=c(11:40),names_to='gene',values_to='mutation')|>
  mutate(mutation=ifelse(is.na(mutation),'Not measured',ifelse(mutation==0,'WT','Mutant')))

mfile <- synGet('syn68710369')$path

pi_metab_hp <- readxl::read_xlsx(mfile, sheet = 'HILIC Positive') %>%
  setNames(make.names(colnames(.))) %>%
  setNames(gsub("_HILIC_POS", "", colnames(.), ignore.case = TRUE))

pi_metab_hn <- readxl::read_xlsx(mfile, sheet = 'HILIC Negative') %>%
  setNames(make.names(colnames(.))) %>%
  setNames(gsub("_HILIC_Neg", "", colnames(.), ignore.case = TRUE))

pi_metab_rp <- readxl::read_xlsx(mfile, sheet = 'RP Positive') %>%
  setNames(make.names(colnames(.))) %>%
  setNames(gsub("_RP_Pos", "", colnames(.), ignore.case = TRUE))

pi_metab_rn <- readxl::read_xlsx(mfile, sheet = 'RP Negative') %>%
  setNames(make.names(colnames(.))) %>%
  setNames(gsub("_RP_Neg", "", colnames(.), ignore.case = TRUE))

# -------------------------------------------------------------------------

# Update Pilot metadata ------------------------------------------------------

##filter so that we only have new samples, no beatAML
new_samp <- subset(sample_exp_metadata, !Accession %in% ba_meta$Accession)  |>
  dplyr::select(c('Accession','Age','Sex',Race,gene,mutation)) |>
  distinct()

#filter for betaaml data that is in the new cohort
ba <- subset(ba_meta,Accession %in% sample_exp_metadata$Accession)  |>
  dplyr::select(Accession, Age,Sex, Race,source) |>
  mutate(Race = ifelse(is.na(Race),'White',Race)) 

##how many samples by race and sex are there in each experiment? 
comb <- new_samp |>
  dplyr::select(Accession,Age,Sex,Race) |>
  mutate(source = 'newSamples')
##what is the age distribution of each experiment?

full_meta_pilot <- rbind(comb,ba) |>
  mutate(Study = 'pilotStudy')

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
  dplyr::left_join(full_meta_pilot %>% 
                     dplyr::select(Accession, Age, Sex, Race) %>%
                     dplyr::distinct(), by = "Accession") %>%
  dplyr::mutate(SampleID = gsub("Exp26", "exp26", SampleID))

all(pi_metadat$Accession %in% sample_exp_metadata$Accession)

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

# Check for duplicate metabolites and append _A, _B, _C, etc labels as needed
ba_metab_hn <- ba_metab_hn %>%
  dplyr::mutate(Name_og = Name, 
                Name = make.names(paste0("hineg_", Name))) %>%
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
ba_hn_dup_df <- ba_metab_hn %>% 
  dplyr::filter(dupcount > 1) %>% 
  dplyr::select(RT..min., m.z, Name_og, Reference.Ion, Name) %>% 
  dplyr::arrange(Name_og)


# Check for duplicate metabolites and append _A, _B, _C, etc labels as needed
ba_metab_rp <- ba_metab_rp %>%
  dplyr::mutate(Name_og = Name, 
                Name = make.names(paste0("rppos_", Name))) %>%
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
ba_rp_dup_df <- ba_metab_rp %>% 
  dplyr::filter(dupcount > 1) %>% 
  dplyr::select(RT..min., m.z, Name_og, Reference.Ion, Name) %>% 
  dplyr::arrange(Name_og)


# Check for duplicate metabolites and append _A, _B, _C, etc labels as needed
ba_metab_rn <- ba_metab_rn %>%
  dplyr::mutate(Name_og = Name, 
                Name = make.names(paste0("rpneg_", Name))) %>%
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
ba_rn_dup_df <- ba_metab_rn %>% 
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
  ), Race_mod = factor(Race_mod, levels = c("White", "Non-white")))

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



# Create initial edata, emeta, and fdata objects
ba_metab_hn_edat <- ba_metab_hn %>%
  dplyr::select(Name, 
                dplyr::contains("BEAT_AML_PNL_"))
ba_metab_hn_emet <- ba_metab_hn %>%
  dplyr::select(-dplyr::contains("BEAT_AML_PNL_"))

ba_metab_hn_fdat <- data.frame(SampleID = colnames(ba_metab_hn_edat)[-1]) %>%
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
  ), Race_mod = factor(Race_mod, levels = c("White", "Non-white")))

# Check consistency of edat and fdat names
all(ba_metab_hn_fdat$SampleID %in% colnames(ba_metab_hn_edat)[-1])

ba_pm_hn <- as.metabData(e_data = ba_metab_hn_edat,
                         f_data = ba_metab_hn_fdat,
                         e_meta = ba_metab_hn_emet,
                         edata_cname = "Name", 
                         fdata_cname = "SampleID",
                         emeta_cname = "Name",
                         data_scale = "abundance",
                         data_types = "HILIC Negative")

ba_pm_hn <- edata_transform(ba_pm_hn, data_scale = "log2")



# Create initial edata, emeta, and fdata objects
ba_metab_rp_edat <- ba_metab_rp %>%
  dplyr::select(Name, 
                dplyr::contains("BEAT_AML_PNL_"))
ba_metab_rp_emet <- ba_metab_rp %>%
  dplyr::select(-dplyr::contains("BEAT_AML_PNL_"))

ba_metab_rp_fdat <- data.frame(SampleID = colnames(ba_metab_rp_edat)[-1]) %>%
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
  ), Race_mod = factor(Race_mod, levels = c("White", "Non-white")))

# Check consistency of edat and fdat names
all(ba_metab_rp_fdat$SampleID %in% colnames(ba_metab_rp_edat)[-1])

ba_pm_rp <- as.metabData(e_data = ba_metab_rp_edat,
                         f_data = ba_metab_rp_fdat,
                         e_meta = ba_metab_rp_emet,
                         edata_cname = "Name", 
                         fdata_cname = "SampleID",
                         emeta_cname = "Name",
                         data_scale = "abundance",
                         data_types = "RP Positive")

ba_pm_rp <- edata_transform(ba_pm_rp, data_scale = "log2")


# Create initial edata, emeta, and fdata objects
ba_metab_rn_edat <- ba_metab_rn %>%
  dplyr::select(Name, 
                dplyr::contains("BEAT_AML_PNL_"))
ba_metab_rn_emet <- ba_metab_rn %>%
  dplyr::select(-dplyr::contains("BEAT_AML_PNL_"))

ba_metab_rn_fdat <- data.frame(SampleID = colnames(ba_metab_rn_edat)[-1]) %>%
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
  ), Race_mod = factor(Race_mod, levels = c("White", "Non-white")))

# Check consistency of edat and fdat names
all(ba_metab_rn_fdat$SampleID %in% colnames(ba_metab_rn_edat)[-1])

ba_pm_rn <- as.metabData(e_data = ba_metab_rn_edat,
                         f_data = ba_metab_rn_fdat,
                         e_meta = ba_metab_rn_emet,
                         edata_cname = "Name", 
                         fdata_cname = "SampleID",
                         emeta_cname = "Name",
                         data_scale = "abundance",
                         data_types = "RP Negative")

ba_pm_rn <- edata_transform(ba_pm_rn, data_scale = "log2")

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


# Check for duplicate metabolites and append _A, _B, _C, etc labels as needed
pi_metab_hn <- pi_metab_hn %>%
  dplyr::mutate(Name_og = Name, 
                Name = make.names(paste0("hineg_", Name))) %>%
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
pi_hn_dup_df <- pi_metab_hn %>% 
  dplyr::filter(dupcount > 1) %>% 
  dplyr::select(RT..min., m.z, Name_og, Reference.Ion, Name) %>% 
  dplyr::arrange(Name_og)


# Check for duplicate metabolites and append _A, _B, _C, etc labels as needed
pi_metab_rp <- pi_metab_rp %>%
  dplyr::mutate(Name_og = Name, 
                Name = make.names(paste0("rppos_", Name))) %>%
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
pi_rp_dup_df <- pi_metab_rp %>% 
  dplyr::filter(dupcount > 1) %>% 
  dplyr::select(RT..min., m.z, Name_og, Reference.Ion, Name) %>% 
  dplyr::arrange(Name_og)


# Check for duplicate metabolites and append _A, _B, _C, etc labels as needed
pi_metab_rn <- pi_metab_rn %>%
  dplyr::mutate(Name_og = Name, 
                Name = make.names(paste0("rpneg_", Name))) %>%
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
pi_rn_dup_df <- pi_metab_rn %>% 
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
  ), Race_mod = factor(Race_mod, levels = c("White", "Non-white")))

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



# Create initial edata, emeta, and fdata objects
pi_metab_hn_edat <- pi_metab_hn %>%
  dplyr::select(Name, 
                dplyr::contains("PTRC_exp26_Met")) %>%
  dplyr::select(-dplyr::contains("Blank")) %>%
  dplyr::select(-dplyr::contains("QC_Pool"))
pi_metab_hn_emet <- pi_metab_hn %>%
  dplyr::select(Name:PTRC_exp26_Met_QC_Pool_07)

pi_metab_hn_fdat <- data.frame(SampleID = colnames(pi_metab_hn_edat)[-1]) %>%
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
  ), Race_mod = factor(Race_mod, levels = c("White", "Non-white")))

# Check consistency of edat and fdat names
all(pi_metab_hn_fdat$SampleID %in% colnames(pi_metab_hn_edat)[-1])

pi_pm_hn <- as.metabData(e_data = pi_metab_hn_edat,
                         f_data = pi_metab_hn_fdat,
                         e_meta = pi_metab_hn_emet,
                         edata_cname = "Name", 
                         fdata_cname = "SampleID",
                         emeta_cname = "Name",
                         data_scale = "abundance",
                         data_types = "HILIC Negative")

pi_pm_hn <- edata_transform(pi_pm_hn, data_scale = "log2")




# Create initial edata, emeta, and fdata objects
pi_metab_rp_edat <- pi_metab_rp %>%
  dplyr::select(Name, 
                dplyr::contains("PTRC_exp26_Met")) %>%
  dplyr::select(-dplyr::contains("Blank")) %>%
  dplyr::select(-dplyr::contains("QC_Pool"))
pi_metab_rp_emet <- pi_metab_rp %>%
  dplyr::select(Name:PTRC_exp26_Met_QC_Pool_07)

pi_metab_rp_fdat <- data.frame(SampleID = colnames(pi_metab_rp_edat)[-1]) %>%
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
  ), Race_mod = factor(Race_mod, levels = c("White", "Non-white")))

# Check consistency of edat and fdat names
all(pi_metab_rp_fdat$SampleID %in% colnames(pi_metab_rp_edat)[-1])

pi_pm_rp <- as.metabData(e_data = pi_metab_rp_edat,
                         f_data = pi_metab_rp_fdat,
                         e_meta = pi_metab_rp_emet,
                         edata_cname = "Name", 
                         fdata_cname = "SampleID",
                         emeta_cname = "Name",
                         data_scale = "abundance",
                         data_types = "RP Positive")

pi_pm_rp <- edata_transform(pi_pm_rp, data_scale = "log2")




# Create initial edata, emeta, and fdata objects
pi_metab_rn_edat <- pi_metab_rn %>%
  dplyr::select(Name, 
                dplyr::contains("PTRC_exp26_Met")) %>%
  dplyr::select(-dplyr::contains("Blank")) %>%
  dplyr::select(-dplyr::contains("QC_Pool"))
pi_metab_rn_emet <- pi_metab_rn %>%
  dplyr::select(Name:PTRC_exp26_Met_QC_Pool_07)

pi_metab_rn_fdat <- data.frame(SampleID = colnames(pi_metab_rn_edat)[-1]) %>%
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
  ), Race_mod = factor(Race_mod, levels = c("White", "Non-white")))

# Check consistency of edat and fdat names
all(pi_metab_rn_fdat$SampleID %in% colnames(pi_metab_rn_edat)[-1])

pi_pm_rn <- as.metabData(e_data = pi_metab_rn_edat,
                         f_data = pi_metab_rn_fdat,
                         e_meta = pi_metab_rn_emet,
                         edata_cname = "Name", 
                         fdata_cname = "SampleID",
                         emeta_cname = "Name",
                         data_scale = "abundance",
                         data_types = "RP Negative")

pi_pm_rn <- edata_transform(pi_pm_rn, data_scale = "log2")

# -------------------------------------------------------------------------


### Flagging Common Metabolites ---------------------------------------------

# Hilic Positive ----------------------------------------------------------

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

# Hilic Negative ----------------------------------------------------------

# Flag based on:
# 1) shared name
# 2) m.z.
# 3) R.T.
# However, this matching may need to be done manually, as is the case for
# lipids.

temp_ba <- ba_pm_hn$e_meta %>%
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

temp_pi <- pi_pm_hn$e_meta %>%
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

bapi_matchkey_hn <- temp_ba %>%
  dplyr::filter(match_Name != "")

rm(cand_pi_matches, temp_ba, temp_pi)

# RP Positive ----------------------------------------------------------

# Flag based on:
# 1) shared name
# 2) m.z.
# 3) R.T.
# However, this matching may need to be done manually, as is the case for
# lipids.

temp_ba <- ba_pm_rp$e_meta %>%
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

temp_pi <- pi_pm_rp$e_meta %>%
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

bapi_matchkey_rp <- temp_ba %>%
  dplyr::filter(match_Name != "")

rm(cand_pi_matches, temp_ba, temp_pi)

# RP Negative ----------------------------------------------------------

# Flag based on:
# 1) shared name
# 2) m.z.
# 3) R.T.
# However, this matching may need to be done manually, as is the case for
# lipids.

temp_ba <- ba_pm_rn$e_meta %>%
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

temp_pi <- pi_pm_rn$e_meta %>%
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

bapi_matchkey_rn <- temp_ba %>%
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

# Normalize data using median centering
ba_pm_hn_norm <- normalize_global(omicsData = ba_pm_hn,
                                  subset_fn = "all",
                                  norm_fn = "median",
                                  apply_norm = TRUE,
                                  backtransform = TRUE)

# Normalize data using median centering
ba_pm_rp_norm <- normalize_global(omicsData = ba_pm_rp,
                                  subset_fn = "all",
                                  norm_fn = "median",
                                  apply_norm = TRUE,
                                  backtransform = TRUE)

# Normalize data using median centering
ba_pm_rn_norm <- normalize_global(omicsData = ba_pm_rn,
                                  subset_fn = "all",
                                  norm_fn = "median",
                                  apply_norm = TRUE,
                                  backtransform = TRUE)

# Combine modes

ba_pm_hp_norm_temp <- as.lipidData(e_data = ba_pm_hp_norm$e_data, 
                                   f_data = ba_pm_hp_norm$f_data, 
                                   e_meta = ba_pm_hp_norm$e_meta,
                                   edata_cname = "Name", 
                                   fdata_cname = "SampleID",
                                   emeta_cname = "Name",
                                   data_scale = "log2",
                                   is_normalized = TRUE,
                                   data_types = "HILIC Positive")

ba_pm_hn_norm_temp <- as.lipidData(e_data = ba_pm_hn_norm$e_data, 
                                   f_data = ba_pm_hn_norm$f_data, 
                                   e_meta = ba_pm_hn_norm$e_meta,
                                   edata_cname = "Name", 
                                   fdata_cname = "SampleID",
                                   emeta_cname = "Name",
                                   data_scale = "log2",
                                   is_normalized = TRUE,
                                   data_types = "HILIC Negative")

ba_pm_hcomb <- combine_omicsData(ba_pm_hp_norm_temp, ba_pm_hn_norm_temp,
                              retain_filters = TRUE)
rm(ba_pm_hp_norm_temp, ba_pm_hn_norm_temp)



ba_pm_rp_norm_temp <- as.lipidData(e_data = ba_pm_rp_norm$e_data, 
                                   f_data = ba_pm_rp_norm$f_data, 
                                   e_meta = ba_pm_rp_norm$e_meta,
                                   edata_cname = "Name", 
                                   fdata_cname = "SampleID",
                                   emeta_cname = "Name",
                                   data_scale = "log2",
                                   is_normalized = TRUE,
                                   data_types = "RP Positive")

ba_pm_rn_norm_temp <- as.lipidData(e_data = ba_pm_rn_norm$e_data, 
                                   f_data = ba_pm_rn_norm$f_data, 
                                   e_meta = ba_pm_rn_norm$e_meta,
                                   edata_cname = "Name", 
                                   fdata_cname = "SampleID",
                                   emeta_cname = "Name",
                                   data_scale = "log2",
                                   is_normalized = TRUE,
                                   data_types = "RP Negative")

ba_pm_rcomb <- combine_omicsData(ba_pm_rp_norm_temp, ba_pm_rn_norm_temp,
                              retain_filters = TRUE)
rm(ba_pm_rp_norm_temp, ba_pm_rn_norm_temp)

# Pilot -------------------------------------------------------------------

# Normalize data using median centering
pi_pm_hp_norm <- normalize_global(omicsData = pi_pm_hp,
                                  subset_fn = "all",
                                  norm_fn = "median",
                                  apply_norm = TRUE,
                                  backtransform = TRUE)

# Normalize data using median centering
pi_pm_hn_norm <- normalize_global(omicsData = pi_pm_hn,
                                  subset_fn = "all",
                                  norm_fn = "median",
                                  apply_norm = TRUE,
                                  backtransform = TRUE)

# Normalize data using median centering
pi_pm_rp_norm <- normalize_global(omicsData = pi_pm_rp,
                                  subset_fn = "all",
                                  norm_fn = "median",
                                  apply_norm = TRUE,
                                  backtransform = TRUE)

# Normalize data using median centering
pi_pm_rn_norm <- normalize_global(omicsData = pi_pm_rn,
                                  subset_fn = "all",
                                  norm_fn = "median",
                                  apply_norm = TRUE,
                                  backtransform = TRUE)

# Combine modes

pi_pm_hp_norm_temp <- as.lipidData(e_data = pi_pm_hp_norm$e_data, 
                                   f_data = pi_pm_hp_norm$f_data, 
                                   e_meta = pi_pm_hp_norm$e_meta,
                                   edata_cname = "Name", 
                                   fdata_cname = "SampleID",
                                   emeta_cname = "Name",
                                   data_scale = "log2",
                                   is_normalized = TRUE,
                                   data_types = "HILIC Positive")

pi_pm_hn_norm_temp <- as.lipidData(e_data = pi_pm_hn_norm$e_data, 
                                   f_data = pi_pm_hn_norm$f_data, 
                                   e_meta = pi_pm_hn_norm$e_meta,
                                   edata_cname = "Name", 
                                   fdata_cname = "SampleID",
                                   emeta_cname = "Name",
                                   data_scale = "log2",
                                   is_normalized = TRUE,
                                   data_types = "HILIC Negative")

pi_pm_hcomb <- combine_omicsData(pi_pm_hp_norm_temp, pi_pm_hn_norm_temp,
                                 retain_filters = TRUE)
rm(pi_pm_hp_norm_temp, pi_pm_hn_norm_temp)



pi_pm_rp_norm_temp <- as.lipidData(e_data = pi_pm_rp_norm$e_data, 
                                   f_data = pi_pm_rp_norm$f_data, 
                                   e_meta = pi_pm_rp_norm$e_meta,
                                   edata_cname = "Name", 
                                   fdata_cname = "SampleID",
                                   emeta_cname = "Name",
                                   data_scale = "log2",
                                   is_normalized = TRUE,
                                   data_types = "RP Positive")

pi_pm_rn_norm_temp <- as.lipidData(e_data = pi_pm_rn_norm$e_data, 
                                   f_data = pi_pm_rn_norm$f_data, 
                                   e_meta = pi_pm_rn_norm$e_meta,
                                   edata_cname = "Name", 
                                   fdata_cname = "SampleID",
                                   emeta_cname = "Name",
                                   data_scale = "log2",
                                   is_normalized = TRUE,
                                   data_types = "RP Negative")

pi_pm_rcomb <- combine_omicsData(pi_pm_rp_norm_temp, pi_pm_rn_norm_temp,
                                 retain_filters = TRUE)
rm(pi_pm_rp_norm_temp, pi_pm_rn_norm_temp)

## Integrate ---------------------------------------------------


# Combined matchkeys according to their modes.
bapi_matchkey_hcomb <- rbind.data.frame(bapi_matchkey_hp,
                                        bapi_matchkey_hn)

bapi_matchkey_rcomb <- rbind.data.frame(bapi_matchkey_rp,
                                        bapi_matchkey_rn)


# HILIC -------------------------------------------------------------------

# Define subsets containing only common lipids

# BeatAML
my_filt <- custom_filter(ba_pm_hcomb, 
                         e_data_keep = bapi_matchkey_hcomb$Name)
ba_pm_hcomb <- applyFilt(filter_object = my_filt,
                           omicsData = ba_pm_hcomb)


# Pilot
my_filt <- custom_filter(pi_pm_hcomb, 
                         e_data_keep = bapi_matchkey_hcomb$match_Name)
pi_pm_hcomb <- applyFilt(filter_object = my_filt,
                           omicsData = pi_pm_hcomb)


# Rename the formatted_name in Pilot to match that of BeatAML
for(i in 1:nrow(bapi_matchkey_hcomb)){
  tempname <- bapi_matchkey_hcomb$Name[i]
  
  temp_repname <- bapi_matchkey_hcomb$match_Name[i]
  repidx <- which(pi_pm_hcomb$e_meta$Name == temp_repname)
  
  pi_pm_hcomb$e_meta$Name[repidx] <- tempname
  pi_pm_hcomb$e_data$Name[repidx] <- tempname
}
all(pi_pm_hcomb$e_data$Name %in% ba_pm_hcomb$e_data$Name) &
  all(ba_pm_hcomb$e_data$Name %in% pi_pm_hcomb$e_data$Name)




# Define new pmartR objects based on a concatenation of e_data, f_data
edat_common_hcomb <- ba_pm_hcomb$e_data %>%
  dplyr::left_join(pi_pm_hcomb$e_data, by = "Name")

fdat_common_hcomb <- rbind.data.frame(ba_pm_hcomb$f_data %>%
                                      dplyr::mutate(batchid = 1),
                                   pi_pm_hcomb$f_data %>%
                                      dplyr::mutate(batchid = 2))

emeta_common_hcomb <- ba_pm_hcomb$e_meta %>%
  dplyr::left_join(pi_pm_hcomb$e_meta, by = "Name") %>%
  setNames(gsub("\\.x", "\\.beataml", colnames(.))) %>%
  setNames(gsub("\\.y", "\\.pilot", colnames(.)))

# all(names(edat_common_hcomb)[-1] %in% fdat_common_hcomb$SampleID)
# setdiff(names(edat_common_hcomb)[-1], fdat_common_hcomb$SampleID)

pm_common_hcomb <- as.lipidData(e_data = edat_common_hcomb, 
                             f_data = fdat_common_hcomb,
                             e_meta = emeta_common_hcomb,
                             emeta_cname = "Name",
                             edata_cname = "Name", 
                             fdata_cname = "SampleID",
                             data_scale = "log2",
                             data_types = "HILIC",
                             is_normalized = TRUE)

# Several samples (13 specifically: 12 Beat AML, 1 Pilot) have been removed due to
# missing race information
# pm_common_hcomb <- group_designation(pm_common_hcomb, main_effects = c("Race_mod"), batch_id = "batchid")
pm_common_hcomb <- group_designation(pm_common_hcomb, main_effects = c("study"), batch_id = "batchid")

pm_common_hcomb_combat <- bc_combat(omicsData = pm_common_hcomb, use_groups = FALSE)

# RP -------------------------------------------------------------------

# Define subsets containing only common lipids

# BeatAML
my_filt <- custom_filter(ba_pm_rcomb, 
                         e_data_keep = bapi_matchkey_rcomb$Name)
ba_pm_rcomb <- applyFilt(filter_object = my_filt,
                         omicsData = ba_pm_rcomb)


# Pilot
my_filt <- custom_filter(pi_pm_rcomb, 
                         e_data_keep = bapi_matchkey_rcomb$match_Name)
pi_pm_rcomb <- applyFilt(filter_object = my_filt,
                         omicsData = pi_pm_rcomb)


# Rename the formatted_name in Pilot to match that of BeatAML
for(i in 1:nrow(bapi_matchkey_rcomb)){
  tempname <- bapi_matchkey_rcomb$Name[i]
  
  temp_repname <- bapi_matchkey_rcomb$match_Name[i]
  repidx <- which(pi_pm_rcomb$e_meta$Name == temp_repname)
  
  pi_pm_rcomb$e_meta$Name[repidx] <- tempname
  pi_pm_rcomb$e_data$Name[repidx] <- tempname
}
all(pi_pm_rcomb$e_data$Name %in% ba_pm_rcomb$e_data$Name) &
  all(ba_pm_rcomb$e_data$Name %in% pi_pm_rcomb$e_data$Name)




# Define new pmartR objects based on a concatenation of e_data, f_data
edat_common_rcomb <- ba_pm_rcomb$e_data %>%
  dplyr::left_join(pi_pm_rcomb$e_data, by = "Name")

fdat_common_rcomb <- rbind.data.frame(ba_pm_rcomb$f_data %>%
                                        dplyr::mutate(batchid = 1),
                                      pi_pm_rcomb$f_data %>%
                                        dplyr::mutate(batchid = 2))

emeta_common_rcomb <- ba_pm_rcomb$e_meta %>%
  dplyr::left_join(pi_pm_rcomb$e_meta, by = "Name") %>%
  setNames(gsub("\\.x", "\\.beataml", colnames(.))) %>%
  setNames(gsub("\\.y", "\\.pilot", colnames(.)))

# all(names(edat_common_rcomb)[-1] %in% fdat_common_rcomb$SampleID)
# setdiff(names(edat_common_rcomb)[-1], fdat_common_rcomb$SampleID)

pm_common_rcomb <- as.lipidData(e_data = edat_common_rcomb, 
                                f_data = fdat_common_rcomb,
                                e_meta = emeta_common_rcomb,
                                emeta_cname = "Name",
                                edata_cname = "Name", 
                                fdata_cname = "SampleID",
                                data_scale = "log2",
                                data_types = "RP",
                                is_normalized = TRUE)

# Several samples (13 specifically: 12 Beat AML, 1 Pilot) have been removed due to
# missing race information
# pm_common_rcomb <- group_designation(pm_common_rcomb, main_effects = c("Race_mod"), batch_id = "batchid")
pm_common_rcomb <- group_designation(pm_common_rcomb, main_effects = c("study"), batch_id = "batchid")

pm_common_rcomb_combat <- bc_combat(omicsData = pm_common_rcomb, use_groups = FALSE)


# -------------------------------------------------------------------------


### HILIC -------------------------------------------------------------------

## RMD Outlier Assessment --------------------------------------------------

# BeatAML -----------------------------------------------------------------

ba_pm_hcomb <- group_designation(ba_pm_hcomb, main_effects = c("study"))

myrmd <- rmd_filter(ba_pm_hcomb)
png(here("figures", "outlier_beataml_hi_rmd.png"), units="in",
    width=10, height=7, res=300)
plot(myrmd, pvalue_threshold = 0.00001)
dev.off()

# Pilot -------------------------------------------------------------------

pi_pm_hcomb <- group_designation(pi_pm_hcomb, main_effects = c("study"))

myrmd <- rmd_filter(pi_pm_hcomb)
png(here("figures", "outlier_pilot_hi_rmd.png"), units="in",
    width=10, height=7, res=300)
plot(myrmd, pvalue_threshold = 0.00001)
dev.off()

# Combined, No BC ---------------------------------------------------------

pm_common_hcomb <- group_designation(pm_common_hcomb, main_effects = c("study"), batch_id = "batchid")

myrmd <- rmd_filter(pm_common_hcomb)
png(here("figures", "outlier_nobc_hi_rmd.png"), units="in",
    width=10, height=7, res=300)
plot(myrmd, pvalue_threshold = 0.00001)
dev.off()

# Combined, BC ------------------------------------------------------------

pm_common_hcomb_combat <- group_designation(pm_common_hcomb_combat, main_effects = c("study"), batch_id = "batchid")

myrmd <- rmd_filter(pm_common_hcomb_combat)
png(here("figures", "outlier_combat_hi_rmd.png"), units="in",
    width=10, height=7, res=300)
plot(myrmd, pvalue_threshold = 0.00001)
dev.off()

# -------------------------------------------------------------------------

## PCA --------------------------------------------------------------------

# BeatAML -----------------------------------------------------------------

ba_pm_hcomb <- group_designation(ba_pm_hcomb, main_effects = c("study"))

mypca <- dim_reduction(ba_pm_hcomb)
png(here("figures", "beataml_hi_pca.png"), units="in",
    width=10, height=7, res=300)
plot(mypca, pvalue_threshold = 0.00001, title_lab = "HP: BeatAML")
dev.off()

# Pilot -------------------------------------------------------------------

pi_pm_hcomb <- group_designation(pi_pm_hcomb, main_effects = c("study"))

mypca <- dim_reduction(pi_pm_hcomb)
png(here("figures", "pilot_hi_pca.png"), units="in",
    width=10, height=7, res=300)
plot(mypca, pvalue_threshold = 0.00001, title_lab = "HP: Pilot")
dev.off()

# Combined, No BC ---------------------------------------------------------

mypca <- dim_reduction(pm_common_hcomb)

png(here("figures", "nobc_hi_pca.png"), units="in",
    width=10, height=7, res=300)
plot(mypca, omicsData = pm_common_hcomb, color_by = "study", 
     title_lab = "HILIC: No Batch-Correction")+theme(text=element_text(size=21))
dev.off()

# Combined, BC ------------------------------------------------------------

mypca <- dim_reduction(pm_common_hcomb_combat)

png(here("figures", "combat_hi_pca.png"), units="in",
    width=10, height=7, res=300)
plot(mypca, omicsData = pm_common_hcomb_combat, color_by = "study", 
     title_lab = "HILIC: ComBat")+theme(text=element_text(size=21))
dev.off()

# -------------------------------------------------------------------------

## ANOVA -------------------------------------------------------------------

# BeatAML -----------------------------------------------------------------

ba_pm_hcomb_temp <- group_designation(ba_pm_hcomb, main_effects = c("Race_mod"))
gdf <- attr(ba_pm_hcomb_temp, "group_DF") %>% 
  dplyr::select(-SampleID) %>% # Remove replicate identifier; only retain design information
  dplyr::distinct()
gdf_comps <- data.frame(Control = c(gdf %>% 
                                      dplyr::filter(Group == "Non-white") %>% .$Group),
                        Test    = c(gdf %>% 
                                      dplyr::filter(Group == "White") %>% .$Group))
ba_hitest_race <- imd_anova(omicsData = ba_pm_hcomb_temp, 
                            comparisons = gdf_comps,
                            test_method = "anova", 
                            pval_adjust_a_multcomp = "holm", 
                            pval_thresh = 0.05)

summary(ba_hitest_race)$sig_table %>%
  dplyr::mutate(Comparison = gsub("_vs_", " vs ", Comparison)) %>%
  dplyr::mutate(dplyr::across(dplyr::matches(":"), ~ paste0(.x, " (",
                                                            round(.x/dim(ba_hitest_race)[1]*100,2), "%)"))) %>%
  dplyr::mutate(Comparison = gsub("X", "", Comparison))

whichsig_ba_hitest_race <- ba_hitest_race %>%
  dplyr::select(Name, `P_value_A_White_vs_Non-white`, `Flag_A_White_vs_Non-white`) %>%
  dplyr::filter(`Flag_A_White_vs_Non-white` != 0)

# ba_pm_hcomb_temp <- group_designation(ba_pm_hcomb, main_effects = c("study"))
# 
# ba_hitest_study <- imd_anova(omicsData = ba_pm_hcomb_temp, 
#                              test_method = "anova", 
#                              pval_adjust_a_multcomp = "holm", 
#                              pval_thresh = 0.05)
# summary(ba_hitest_study)

# Pilot -------------------------------------------------------------------

pi_pm_hcomb_temp <- group_designation(pi_pm_hcomb, main_effects = c("Race_mod"))
gdf <- attr(pi_pm_hcomb_temp, "group_DF") %>% 
  dplyr::select(-SampleID) %>% # Remove replicate identifier; only retain design information
  dplyr::distinct()
gdf_comps <- data.frame(Control = c(gdf %>% 
                                      dplyr::filter(Group == "Non-white") %>% .$Group),
                        Test    = c(gdf %>% 
                                      dplyr::filter(Group == "White") %>% .$Group))
pi_hitest_race <- imd_anova(omicsData = pi_pm_hcomb_temp,
                            comparisons = gdf_comps,
                            test_method = "anova",
                            pval_adjust_a_multcomp = "holm",
                            pval_thresh = 0.05)

summary(pi_hitest_race)$sig_table %>%
  dplyr::mutate(Comparison = gsub("_vs_", " vs ", Comparison)) %>%
  dplyr::mutate(dplyr::across(dplyr::matches(":"), ~ paste0(.x, " (",
                                                            round(.x/dim(pi_hitest_race)[1]*100,2), "%)"))) %>%
  dplyr::mutate(Comparison = gsub("X", "", Comparison))

whichsig_pi_hitest_race <- pi_hitest_race %>%
  dplyr::select(Name, `P_value_A_White_vs_Non-white`, `Flag_A_White_vs_Non-white`) %>%
  dplyr::filter(`Flag_A_White_vs_Non-white` != 0)

# pi_pm_hcomb_temp <- group_designation(pi_pm_hcomb, main_effects = c("study"))
# 
# pi_hitest_study <- imd_anova(omicsData = pi_pm_hcomb_temp, 
#                              test_method = "anova", 
#                              pval_adjust_a_multcomp = "holm", 
#                              pval_thresh = 0.05)
# summary(pi_hitest_study)

# Combined, No BC ---------------------------------------------------------

pm_common_hcomb_temp <- group_designation(pm_common_hcomb, main_effects = c("Race_mod"))
gdf <- attr(pm_common_hcomb_temp, "group_DF") %>% 
  dplyr::select(-SampleID) %>% # Remove replicate identifier; only retain design information
  dplyr::distinct()
gdf_comps <- data.frame(Control = c(gdf %>% 
                                      dplyr::filter(Group == "Non-white") %>% .$Group),
                        Test    = c(gdf %>% 
                                      dplyr::filter(Group == "White") %>% .$Group))
nobc_hitest_race <- imd_anova(omicsData = pm_common_hcomb_temp,
                              comparisons = gdf_comps,
                              test_method = "anova",
                              pval_adjust_a_multcomp = "holm",
                              pval_thresh = 0.05)
summary(nobc_hitest_race)$sig_table %>%
  dplyr::mutate(Comparison = gsub("_vs_", " vs ", Comparison)) %>%
  dplyr::mutate(dplyr::across(dplyr::matches(":"), ~ paste0(.x, " (",
                                                            round(.x/dim(nobc_hitest_race)[1]*100,2), "%)"))) %>%
  dplyr::mutate(Comparison = gsub("X", "", Comparison))
whichsig_nobc_hitest_race <- nobc_hitest_race %>%
  dplyr::select(Name, `P_value_A_White_vs_Non-white`, `Flag_A_White_vs_Non-white`) %>%
  dplyr::filter(`Flag_A_White_vs_Non-white` != 0)

pm_common_hcomb_temp <- group_designation(pm_common_hcomb, main_effects = c("study"))
gdf <- attr(pm_common_hcomb_temp, "group_DF") %>% 
  dplyr::select(-SampleID) %>% # Remove replicate identifier; only retain design information
  dplyr::distinct()
gdf_comps <- data.frame(Control = c(gdf %>% 
                                      dplyr::filter(Group == "BeatAML") %>% .$Group),
                        Test    = c(gdf %>% 
                                      dplyr::filter(Group == "pilot") %>% .$Group))
nobc_hitest_study <- imd_anova(omicsData = pm_common_hcomb_temp,
                               comparisons = gdf_comps,
                               test_method = "anova",
                               pval_adjust_a_multcomp = "holm",
                               pval_thresh = 0.05)
summary(nobc_hitest_study)$sig_table %>%
  dplyr::mutate(Comparison = gsub("_vs_", " vs ", Comparison)) %>%
  dplyr::mutate(dplyr::across(dplyr::matches(":"), ~ paste0(.x, " (",
                                                            round(.x/dim(nobc_hitest_study)[1]*100,2), "%)"))) %>%
  dplyr::mutate(Comparison = gsub("X", "", Comparison))
whichsig_nobc_hitest_study <- nobc_hitest_study %>%
  dplyr::select(Name, `P_value_A_pilot_vs_BeatAML`, `Flag_A_pilot_vs_BeatAML`) %>%
  dplyr::filter(`Flag_A_pilot_vs_BeatAML` != 0)

# Combined, BC ------------------------------------------------------------

temp <- pm_common_hcomb_combat
sampkeep <- temp$f_data %>%
  dplyr::filter(Race %in% c("White", "Black")) %>%
  .$SampleID
myfilt <- custom_filter(temp, f_data_keep = sampkeep)
temp <- applyFilt(myfilt, temp)

temp <- group_designation(temp, main_effects = c("Race"))
gdf <- attr(temp, "group_DF") %>% 
  dplyr::select(-SampleID) %>% # Remove replicate identifier; only retain design information
  dplyr::distinct()
gdf_comps <- data.frame(Control = c(gdf %>% 
                                      dplyr::filter(Group == "White") %>% .$Group),
                        Test    = c(gdf %>% 
                                      dplyr::filter(Group == "Black") %>% .$Group))
combat_hitest_race <- imd_anova(omicsData = temp,
                                comparisons = gdf_comps,
                                test_method = "anova",
                                pval_adjust_a_multcomp = "none",
                                pval_thresh = 0.05)

writexl::write_xlsx(x = list(`Black vs White` = combat_hitest_race %>%
                               dplyr::left_join(temp$e_meta %>%
                                                  dplyr::select(Name, Name_og.beataml) %>%
                                                  dplyr::rename(unformatted_name = Name_og.beataml)) %>%
                               dplyr::relocate(unformatted_name, .after = Name)),
                    path = here::here("beataml_pilot_bvsw_hilic_metabolites.xlsx"))


pm_common_hcomb_combat_temp <- group_designation(pm_common_hcomb_combat, main_effects = c("Race_mod"))
gdf <- attr(pm_common_hcomb_combat_temp, "group_DF") %>% 
  dplyr::select(-SampleID) %>% # Remove replicate identifier; only retain design information
  dplyr::distinct()
gdf_comps <- data.frame(Control = c(gdf %>% 
                                      dplyr::filter(Group == "Non-white") %>% .$Group),
                        Test    = c(gdf %>% 
                                      dplyr::filter(Group == "White") %>% .$Group))
combat_hitest_race <- imd_anova(omicsData = pm_common_hcomb_combat_temp,
                                comparisons = gdf_comps,
                                test_method = "anova",
                                pval_adjust_a_multcomp = "holm",
                                pval_thresh = 0.05)
summary(combat_hitest_race)$sig_table %>%
  dplyr::mutate(Comparison = gsub("_vs_", " vs ", Comparison)) %>%
  dplyr::mutate(dplyr::across(dplyr::matches(":"), ~ paste0(.x, " (",
                                                            round(.x/dim(combat_hitest_race)[1]*100,2), "%)"))) %>%
  dplyr::mutate(Comparison = gsub("X", "", Comparison))
whichsig_combat_hitest_race <- combat_hitest_race %>%
  dplyr::select(Name, `P_value_A_White_vs_Non-white`, `Flag_A_White_vs_Non-white`) %>%
  dplyr::filter(`Flag_A_White_vs_Non-white` != 0)

pm_common_hcomb_combat_temp <- group_designation(pm_common_hcomb_combat, main_effects = c("study"))
gdf <- attr(pm_common_hcomb_combat_temp, "group_DF") %>% 
  dplyr::select(-SampleID) %>% # Remove replicate identifier; only retain design information
  dplyr::distinct()
gdf_comps <- data.frame(Control = c(gdf %>% 
                                      dplyr::filter(Group == "BeatAML") %>% .$Group),
                        Test    = c(gdf %>% 
                                      dplyr::filter(Group == "pilot") %>% .$Group))
combat_hitest_study <- imd_anova(omicsData = pm_common_hcomb_combat_temp,
                                 comparisons = gdf_comps,
                                 test_method = "anova",
                                 pval_adjust_a_multcomp = "holm",
                                 pval_thresh = 0.05)
summary(combat_hitest_study)$sig_table %>%
  dplyr::mutate(Comparison = gsub("_vs_", " vs ", Comparison)) %>%
  dplyr::mutate(dplyr::across(dplyr::matches(":"), ~ paste0(.x, " (",
                                                            round(.x/dim(combat_hitest_study)[1]*100,2), "%)"))) %>%
  dplyr::mutate(Comparison = gsub("X", "", Comparison))
whichsig_combat_hitest_study <- combat_hitest_study %>%
  dplyr::select(Name, `P_value_A_pilot_vs_BeatAML`, `Flag_A_pilot_vs_BeatAML`) %>%
  dplyr::filter(`Flag_A_pilot_vs_BeatAML` != 0)

# -------------------------------------------------------------------------

## Significant Set Comparisons ---------------------------------------------

# BeatAML vs Pilot --------------------------------------------------------

# Commonly significant, Race Comparison
intersect(whichsig_ba_hitest_race$Name, whichsig_pi_hitest_race$Name)

# Uniquely significant, BeatAML, Race Comparison
whichsig_ba_hitest_race %>%
  dplyr::filter(!(Name %in% intersect(whichsig_ba_hitest_race$Name, whichsig_pi_hitest_race$Name))) %>%
  .$Name

# Uniquely significant, Pilot, Race Comparison
whichsig_pi_hitest_race %>%
  dplyr::filter(!(Name %in% intersect(whichsig_ba_hitest_race$Name, whichsig_pi_hitest_race$Name))) %>%
  .$Name

# BeatAML vs Combined, No BC ----------------------------------------------

# Commonly significant, Race Comparison
intersect(whichsig_ba_hitest_race$Name, whichsig_nobc_hitest_race$Name)

# Uniquely significant, BeatAML, Race Comparison
whichsig_ba_hitest_race %>%
  dplyr::filter(!(Name %in% intersect(whichsig_ba_hitest_race$Name, whichsig_pi_hitest_race$Name))) %>%
  .$Name

# Uniquely significant, Combined No BC, Race Comparison
whichsig_nobc_hitest_race %>%
  dplyr::filter(!(Name %in% intersect(whichsig_ba_hitest_race$Name, whichsig_nobc_hitest_race$Name))) %>%
  .$Name

# BeatAML vs Combined, BC -------------------------------------------------

# Commonly significant, Race Comparison
intersect(whichsig_ba_hitest_race$Name, whichsig_combat_hitest_race$Name)

# Uniquely significant, BeatAML, Race Comparison
whichsig_ba_hitest_race %>%
  dplyr::filter(!(Name %in% intersect(whichsig_ba_hitest_race$Name, whichsig_combat_hitest_race$Name))) %>%
  .$Name

# Uniquely significant, Combined BC, Race Comparison
whichsig_combat_hitest_race %>%
  dplyr::filter(!(Name %in% intersect(whichsig_ba_hitest_race$Name, whichsig_combat_hitest_race$Name))) %>%
  .$Name

# Pilot vs Combined, No BC ----------------------------------------------

# Commonly significant, Race Comparison
intersect(whichsig_pi_hitest_race$Name, whichsig_nobc_hitest_race$Name)

# Uniquely significant, Pilot, Race Comparison
whichsig_pi_hitest_race %>%
  dplyr::filter(!(Name %in% intersect(whichsig_pi_hitest_race$Name, whichsig_pi_hitest_race$Name))) %>%
  .$Name

# Uniquely significant, Combined No BC, Race Comparison
whichsig_nobc_hitest_race %>%
  dplyr::filter(!(Name %in% intersect(whichsig_pi_hitest_race$Name, whichsig_nobc_hitest_race$Name))) %>%
  .$Name

# Pilot vs Combined, BC -------------------------------------------------

# Commonly significant, Race Comparison
intersect(whichsig_pi_hitest_race$Name, whichsig_combat_hitest_race$Name)

# Uniquely significant, Pilot, Race Comparison
whichsig_pi_hitest_race %>%
  dplyr::filter(!(Name %in% intersect(whichsig_pi_hitest_race$Name, whichsig_combat_hitest_race$Name))) %>%
  .$Name

# Uniquely significant, Combined BC, Race Comparison
whichsig_combat_hitest_race %>%
  dplyr::filter(!(Name %in% intersect(whichsig_pi_hitest_race$Name, whichsig_combat_hitest_race$Name))) %>%
  .$Name

# Combined, No BC vs Combined, BC -------------------------------------------------

# Commonly significant, Race Comparison
intersect(whichsig_nobc_hitest_race$Name, whichsig_combat_hitest_race$Name)

# Uniquely significant, Combined No BC, Race Comparison
whichsig_nobc_hitest_race %>%
  dplyr::filter(!(Name %in% intersect(whichsig_nobc_hitest_race$Name, whichsig_combat_hitest_race$Name))) %>%
  .$Name

# Uniquely significant, Combined BC, Race Comparison
whichsig_combat_hitest_race %>%
  dplyr::filter(!(Name %in% intersect(whichsig_nobc_hitest_race$Name, whichsig_combat_hitest_race$Name))) %>%
  .$Name

# -------------------------------------------------------------------------


### RP -------------------------------------------------------------------

## RMD Outlier Assessment --------------------------------------------------

# BeatAML -----------------------------------------------------------------

ba_pm_rcomb <- group_designation(ba_pm_rcomb, main_effects = c("study"))

myrmd <- rmd_filter(ba_pm_rcomb)
png(here("figures", "outlier_beataml_r_rmd.png"), units="in",
    width=10, height=7, res=300)
plot(myrmd, pvalue_threshold = 0.00001)
dev.off()

# Pilot -------------------------------------------------------------------

pi_pm_rcomb <- group_designation(pi_pm_rcomb, main_effects = c("study"))

myrmd <- rmd_filter(pi_pm_rcomb)
png(here("figures", "outlier_pilot_r_rmd.png"), units="in",
    width=10, height=7, res=300)
plot(myrmd, pvalue_threshold = 0.00001)
dev.off()

# Combined, No BC ---------------------------------------------------------

pm_common_rcomb <- group_designation(pm_common_rcomb, main_effects = c("study"), batch_id = "batchid")

myrmd <- rmd_filter(pm_common_rcomb)
png(here("figures", "outlier_nobc_r_rmd.png"), units="in",
    width=10, height=7, res=300)
plot(myrmd, pvalue_threshold = 0.00001)
dev.off()

# Combined, BC ------------------------------------------------------------

pm_common_rcomb_combat <- group_designation(pm_common_rcomb_combat, main_effects = c("study"), batch_id = "batchid")

myrmd <- rmd_filter(pm_common_rcomb_combat)
png(here("figures", "outlier_combat_r_rmd.png"), units="in",
    width=10, height=7, res=300)
plot(myrmd, pvalue_threshold = 0.00001)
dev.off()

# -------------------------------------------------------------------------

## PCA --------------------------------------------------------------------

# BeatAML -----------------------------------------------------------------

ba_pm_rcomb <- group_designation(ba_pm_rcomb, main_effects = c("study"))

mypca <- dim_reduction(ba_pm_rcomb)
png(here("figures", "beataml_r_pca.png"), units="in",
    width=10, height=7, res=300)
plot(mypca, pvalue_threshold = 0.00001, title_lab = "HP: BeatAML")
dev.off()

# Pilot -------------------------------------------------------------------

pi_pm_rcomb <- group_designation(pi_pm_rcomb, main_effects = c("study"))

mypca <- dim_reduction(pi_pm_rcomb)
png(here("figures", "pilot_r_pca.png"), units="in",
    width=10, height=7, res=300)
plot(mypca, pvalue_threshold = 0.00001, title_lab = "HP: Pilot")
dev.off()

# Combined, No BC ---------------------------------------------------------

mypca <- dim_reduction(pm_common_rcomb)

png(here("figures", "nobc_r_pca.png"), units="in",
    width=10, height=7, res=300)
plot(mypca, omicsData = pm_common_rcomb, color_by = "study", 
     title_lab = "RP: No Batch-Correction")+theme(text=element_text(size=21))
dev.off()

# Combined, BC ------------------------------------------------------------

mypca <- dim_reduction(pm_common_rcomb_combat)

png(here("figures", "combat_r_pca.png"), units="in",
    width=10, height=7, res=300)
plot(mypca, omicsData = pm_common_rcomb_combat, color_by = "study", 
     title_lab = "RP: ComBat")+theme(text=element_text(size=21))
dev.off()

# -------------------------------------------------------------------------

## ANOVA -------------------------------------------------------------------

# BeatAML -----------------------------------------------------------------

ba_pm_rcomb_temp <- group_designation(ba_pm_rcomb, main_effects = c("Race_mod"))
gdf <- attr(ba_pm_rcomb_temp, "group_DF") %>% 
  dplyr::select(-SampleID) %>% # Remove replicate identifier; only retain design information
  dplyr::distinct()
gdf_comps <- data.frame(Control = c(gdf %>% 
                                      dplyr::filter(Group == "Non-white") %>% .$Group),
                        Test    = c(gdf %>% 
                                      dplyr::filter(Group == "White") %>% .$Group))
ba_rtest_race <- imd_anova(omicsData = ba_pm_rcomb_temp, 
                            comparisons = gdf_comps,
                            test_method = "anova", 
                            pval_adjust_a_multcomp = "holm", 
                            pval_thresh = 0.05)

summary(ba_rtest_race)$sig_table %>%
  dplyr::mutate(Comparison = gsub("_vs_", " vs ", Comparison)) %>%
  dplyr::mutate(dplyr::across(dplyr::matches(":"), ~ paste0(.x, " (",
                                                            round(.x/dim(ba_rtest_race)[1]*100,2), "%)"))) %>%
  dplyr::mutate(Comparison = gsub("X", "", Comparison))

whichsig_ba_rtest_race <- ba_rtest_race %>%
  dplyr::select(Name, `P_value_A_White_vs_Non-white`, `Flag_A_White_vs_Non-white`) %>%
  dplyr::filter(`Flag_A_White_vs_Non-white` != 0)

# ba_pm_rcomb_temp <- group_designation(ba_pm_rcomb, main_effects = c("study"))
# 
# ba_rtest_study <- imd_anova(omicsData = ba_pm_rcomb_temp, 
#                              test_method = "anova", 
#                              pval_adjust_a_multcomp = "holm", 
#                              pval_thresh = 0.05)
# summary(ba_rtest_study)

# Pilot -------------------------------------------------------------------

pi_pm_rcomb_temp <- group_designation(pi_pm_rcomb, main_effects = c("Race_mod"))
gdf <- attr(pi_pm_rcomb_temp, "group_DF") %>% 
  dplyr::select(-SampleID) %>% # Remove replicate identifier; only retain design information
  dplyr::distinct()
gdf_comps <- data.frame(Control = c(gdf %>% 
                                      dplyr::filter(Group == "Non-white") %>% .$Group),
                        Test    = c(gdf %>% 
                                      dplyr::filter(Group == "White") %>% .$Group))
pi_rtest_race <- imd_anova(omicsData = pi_pm_rcomb_temp,
                            comparisons = gdf_comps,
                            test_method = "anova",
                            pval_adjust_a_multcomp = "holm",
                            pval_thresh = 0.05)

summary(pi_rtest_race)$sig_table %>%
  dplyr::mutate(Comparison = gsub("_vs_", " vs ", Comparison)) %>%
  dplyr::mutate(dplyr::across(dplyr::matches(":"), ~ paste0(.x, " (",
                                                            round(.x/dim(pi_rtest_race)[1]*100,2), "%)"))) %>%
  dplyr::mutate(Comparison = gsub("X", "", Comparison))

whichsig_pi_rtest_race <- pi_rtest_race %>%
  dplyr::select(Name, `P_value_A_White_vs_Non-white`, `Flag_A_White_vs_Non-white`) %>%
  dplyr::filter(`Flag_A_White_vs_Non-white` != 0)

# pi_pm_rcomb_temp <- group_designation(pi_pm_rcomb, main_effects = c("study"))
# 
# pi_rtest_study <- imd_anova(omicsData = pi_pm_rcomb_temp, 
#                              test_method = "anova", 
#                              pval_adjust_a_multcomp = "holm", 
#                              pval_thresh = 0.05)
# summary(pi_rtest_study)

# Combined, No BC ---------------------------------------------------------

pm_common_rcomb_temp <- group_designation(pm_common_rcomb, main_effects = c("Race_mod"))
gdf <- attr(pm_common_rcomb_temp, "group_DF") %>% 
  dplyr::select(-SampleID) %>% # Remove replicate identifier; only retain design information
  dplyr::distinct()
gdf_comps <- data.frame(Control = c(gdf %>% 
                                      dplyr::filter(Group == "Non-white") %>% .$Group),
                        Test    = c(gdf %>% 
                                      dplyr::filter(Group == "White") %>% .$Group))
nobc_rtest_race <- imd_anova(omicsData = pm_common_rcomb_temp,
                              comparisons = gdf_comps,
                              test_method = "anova",
                              pval_adjust_a_multcomp = "holm",
                              pval_thresh = 0.05)
summary(nobc_rtest_race)$sig_table %>%
  dplyr::mutate(Comparison = gsub("_vs_", " vs ", Comparison)) %>%
  dplyr::mutate(dplyr::across(dplyr::matches(":"), ~ paste0(.x, " (",
                                                            round(.x/dim(nobc_rtest_race)[1]*100,2), "%)"))) %>%
  dplyr::mutate(Comparison = gsub("X", "", Comparison))
whichsig_nobc_rtest_race <- nobc_rtest_race %>%
  dplyr::select(Name, `P_value_A_White_vs_Non-white`, `Flag_A_White_vs_Non-white`) %>%
  dplyr::filter(`Flag_A_White_vs_Non-white` != 0)

pm_common_rcomb_temp <- group_designation(pm_common_rcomb, main_effects = c("study"))
gdf <- attr(pm_common_rcomb_temp, "group_DF") %>% 
  dplyr::select(-SampleID) %>% # Remove replicate identifier; only retain design information
  dplyr::distinct()
gdf_comps <- data.frame(Control = c(gdf %>% 
                                      dplyr::filter(Group == "BeatAML") %>% .$Group),
                        Test    = c(gdf %>% 
                                      dplyr::filter(Group == "pilot") %>% .$Group))
nobc_rtest_study <- imd_anova(omicsData = pm_common_rcomb_temp,
                               comparisons = gdf_comps,
                               test_method = "anova",
                               pval_adjust_a_multcomp = "holm",
                               pval_thresh = 0.05)
summary(nobc_rtest_study)$sig_table %>%
  dplyr::mutate(Comparison = gsub("_vs_", " vs ", Comparison)) %>%
  dplyr::mutate(dplyr::across(dplyr::matches(":"), ~ paste0(.x, " (",
                                                            round(.x/dim(nobc_rtest_study)[1]*100,2), "%)"))) %>%
  dplyr::mutate(Comparison = gsub("X", "", Comparison))
whichsig_nobc_rtest_study <- nobc_rtest_study %>%
  dplyr::select(Name, `P_value_A_pilot_vs_BeatAML`, `Flag_A_pilot_vs_BeatAML`) %>%
  dplyr::filter(`Flag_A_pilot_vs_BeatAML` != 0)

# Combined, BC ------------------------------------------------------------

temp <- pm_common_rcomb_combat
sampkeep <- temp$f_data %>%
  dplyr::filter(Race %in% c("White", "Black")) %>%
  .$SampleID
myfilt <- custom_filter(temp, f_data_keep = sampkeep)
temp <- applyFilt(myfilt, temp)

temp <- group_designation(temp, main_effects = c("Race"))
gdf <- attr(temp, "group_DF") %>% 
  dplyr::select(-SampleID) %>% # Remove replicate identifier; only retain design information
  dplyr::distinct()
gdf_comps <- data.frame(Control = c(gdf %>% 
                                      dplyr::filter(Group == "White") %>% .$Group),
                        Test    = c(gdf %>% 
                                      dplyr::filter(Group == "Black") %>% .$Group))
combat_rtest_race <- imd_anova(omicsData = temp,
                                comparisons = gdf_comps,
                                test_method = "anova",
                                pval_adjust_a_multcomp = "none",
                                pval_thresh = 0.05)

writexl::write_xlsx(x = list(`Black vs White` = combat_rtest_race %>%
                               dplyr::left_join(temp$e_meta %>%
                                                  dplyr::select(Name, Name_og.beataml) %>%
                                                  dplyr::rename(unformatted_name = Name_og.beataml)) %>%
                               dplyr::relocate(unformatted_name, .after = Name)),
                    path = here::here("beataml_pilot_bvsw_rp_metabolites.xlsx"))


pm_common_rcomb_combat_temp <- group_designation(pm_common_rcomb_combat, main_effects = c("Race_mod"))
gdf <- attr(pm_common_rcomb_combat_temp, "group_DF") %>% 
  dplyr::select(-SampleID) %>% # Remove replicate identifier; only retain design information
  dplyr::distinct()
gdf_comps <- data.frame(Control = c(gdf %>% 
                                      dplyr::filter(Group == "Non-white") %>% .$Group),
                        Test    = c(gdf %>% 
                                      dplyr::filter(Group == "White") %>% .$Group))
combat_rtest_race <- imd_anova(omicsData = pm_common_rcomb_combat_temp,
                                comparisons = gdf_comps,
                                test_method = "anova",
                                pval_adjust_a_multcomp = "holm",
                                pval_thresh = 0.05)
summary(combat_rtest_race)$sig_table %>%
  dplyr::mutate(Comparison = gsub("_vs_", " vs ", Comparison)) %>%
  dplyr::mutate(dplyr::across(dplyr::matches(":"), ~ paste0(.x, " (",
                                                            round(.x/dim(combat_rtest_race)[1]*100,2), "%)"))) %>%
  dplyr::mutate(Comparison = gsub("X", "", Comparison))
whichsig_combat_rtest_race <- combat_rtest_race %>%
  dplyr::select(Name, `P_value_A_White_vs_Non-white`, `Flag_A_White_vs_Non-white`) %>%
  dplyr::filter(`Flag_A_White_vs_Non-white` != 0)

pm_common_rcomb_combat_temp <- group_designation(pm_common_rcomb_combat, main_effects = c("study"))
gdf <- attr(pm_common_rcomb_combat_temp, "group_DF") %>% 
  dplyr::select(-SampleID) %>% # Remove replicate identifier; only retain design information
  dplyr::distinct()
gdf_comps <- data.frame(Control = c(gdf %>% 
                                      dplyr::filter(Group == "BeatAML") %>% .$Group),
                        Test    = c(gdf %>% 
                                      dplyr::filter(Group == "pilot") %>% .$Group))
combat_rtest_study <- imd_anova(omicsData = pm_common_rcomb_combat_temp,
                                 comparisons = gdf_comps,
                                 test_method = "anova",
                                 pval_adjust_a_multcomp = "holm",
                                 pval_thresh = 0.05)
summary(combat_rtest_study)$sig_table %>%
  dplyr::mutate(Comparison = gsub("_vs_", " vs ", Comparison)) %>%
  dplyr::mutate(dplyr::across(dplyr::matches(":"), ~ paste0(.x, " (",
                                                            round(.x/dim(combat_rtest_study)[1]*100,2), "%)"))) %>%
  dplyr::mutate(Comparison = gsub("X", "", Comparison))
whichsig_combat_rtest_study <- combat_rtest_study %>%
  dplyr::select(Name, `P_value_A_pilot_vs_BeatAML`, `Flag_A_pilot_vs_BeatAML`) %>%
  dplyr::filter(`Flag_A_pilot_vs_BeatAML` != 0)

# -------------------------------------------------------------------------

## Significant Set Comparisons ---------------------------------------------

# BeatAML vs Pilot --------------------------------------------------------

# Commonly significant, Race Comparison
intersect(whichsig_ba_rtest_race$Name, whichsig_pi_rtest_race$Name)

# Uniquely significant, BeatAML, Race Comparison
whichsig_ba_rtest_race %>%
  dplyr::filter(!(Name %in% intersect(whichsig_ba_rtest_race$Name, whichsig_pi_rtest_race$Name))) %>%
  .$Name

# Uniquely significant, Pilot, Race Comparison
whichsig_pi_rtest_race %>%
  dplyr::filter(!(Name %in% intersect(whichsig_ba_rtest_race$Name, whichsig_pi_rtest_race$Name))) %>%
  .$Name

# BeatAML vs Combined, No BC ----------------------------------------------

# Commonly significant, Race Comparison
intersect(whichsig_ba_rtest_race$Name, whichsig_nobc_rtest_race$Name)

# Uniquely significant, BeatAML, Race Comparison
whichsig_ba_rtest_race %>%
  dplyr::filter(!(Name %in% intersect(whichsig_ba_rtest_race$Name, whichsig_pi_rtest_race$Name))) %>%
  .$Name

# Uniquely significant, Combined No BC, Race Comparison
whichsig_nobc_rtest_race %>%
  dplyr::filter(!(Name %in% intersect(whichsig_ba_rtest_race$Name, whichsig_nobc_rtest_race$Name))) %>%
  .$Name

# BeatAML vs Combined, BC -------------------------------------------------

# Commonly significant, Race Comparison
intersect(whichsig_ba_rtest_race$Name, whichsig_combat_rtest_race$Name)

# Uniquely significant, BeatAML, Race Comparison
whichsig_ba_rtest_race %>%
  dplyr::filter(!(Name %in% intersect(whichsig_ba_rtest_race$Name, whichsig_combat_rtest_race$Name))) %>%
  .$Name

# Uniquely significant, Combined BC, Race Comparison
whichsig_combat_rtest_race %>%
  dplyr::filter(!(Name %in% intersect(whichsig_ba_rtest_race$Name, whichsig_combat_rtest_race$Name))) %>%
  .$Name

# Pilot vs Combined, No BC ----------------------------------------------

# Commonly significant, Race Comparison
intersect(whichsig_pi_rtest_race$Name, whichsig_nobc_rtest_race$Name)

# Uniquely significant, Pilot, Race Comparison
whichsig_pi_rtest_race %>%
  dplyr::filter(!(Name %in% intersect(whichsig_pi_rtest_race$Name, whichsig_pi_rtest_race$Name))) %>%
  .$Name

# Uniquely significant, Combined No BC, Race Comparison
whichsig_nobc_rtest_race %>%
  dplyr::filter(!(Name %in% intersect(whichsig_pi_rtest_race$Name, whichsig_nobc_rtest_race$Name))) %>%
  .$Name

# Pilot vs Combined, BC -------------------------------------------------

# Commonly significant, Race Comparison
intersect(whichsig_pi_rtest_race$Name, whichsig_combat_rtest_race$Name)

# Uniquely significant, Pilot, Race Comparison
whichsig_pi_rtest_race %>%
  dplyr::filter(!(Name %in% intersect(whichsig_pi_rtest_race$Name, whichsig_combat_rtest_race$Name))) %>%
  .$Name

# Uniquely significant, Combined BC, Race Comparison
whichsig_combat_rtest_race %>%
  dplyr::filter(!(Name %in% intersect(whichsig_pi_rtest_race$Name, whichsig_combat_rtest_race$Name))) %>%
  .$Name

# Combined, No BC vs Combined, BC -------------------------------------------------

# Commonly significant, Race Comparison
intersect(whichsig_nobc_rtest_race$Name, whichsig_combat_rtest_race$Name)

# Uniquely significant, Combined No BC, Race Comparison
whichsig_nobc_rtest_race %>%
  dplyr::filter(!(Name %in% intersect(whichsig_nobc_rtest_race$Name, whichsig_combat_rtest_race$Name))) %>%
  .$Name

# Uniquely significant, Combined BC, Race Comparison
whichsig_combat_rtest_race %>%
  dplyr::filter(!(Name %in% intersect(whichsig_nobc_rtest_race$Name, whichsig_combat_rtest_race$Name))) %>%
  .$Name

# -------------------------------------------------------------------------


