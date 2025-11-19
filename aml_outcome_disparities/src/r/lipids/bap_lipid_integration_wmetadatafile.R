library(magrittr)
library(pmartR)
library(malbacR)
library(igraph)
library(here)
library(ggplot2)
library(plotly)
library(synapser)
#here::i_am("BeatAMLPilot_Integration/bap_lipid_integration_wmetadatafile.R")

synapser::synLogin()

# Here we explore integration of BeatAML and Pilot Lipidoimcs datasets via
# ComBat, which is an approach that has historically done best in this regard. 

# Also included are various statistical analyses based on different versions of the data,
# i.e. 1) Pilot data alone, 2) BeatAML data alone, and 3) Integrated data

# We also apply some pre-processing to these data as needed. 


### Data Import -------------------------------------------------------------


## Metadata for BeatAML+Pilot -------------------------------------------------------------------------

beatpilot_metadata <- read.csv(synapser::synGet('syn69692583')$path) %>%
  dplyr::select(Accession:Study) %>%
  dplyr::distinct() %>%
  dplyr::rename(Age2 = Age,
                Sex2 = Sex,
                Race2 = Race,
                source2 = source,
                Study2 = Study)

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

#new pilot metaadata
sample_exp_metadata <- readxl::read_xlsx(synGet('syn68522106')$path) |>
  dplyr::rename(CEBPA = 'CEBPA_BZIP',FLT3 = 'FLT3_TKD',Race='Race_label')|>
  mutate(Sex = ifelse(`Sex (1-male, 0-female)` == 1,'Male','Female'))|>
  tidyr::pivot_longer(cols=c(11:40),names_to='gene',values_to='mutation')|>
  mutate(mutation=ifelse(is.na(mutation),'Not measured',ifelse(mutation==0,'WT','Mutant')))

## BeatAML Data -----------------------------------------------------

# See beat_lipids_proc for why local file is used. Essentially, more extensive preprocessing
# for lipidomics data was required.
ba_something <- readxl::read_xlsx(synapser::synGet('syn71210151')$path,
  #here("BeatAMLPilot_Integration", "data", "processed", "beataml_lipids_log2_sum.xlsx"),
                             sheet = "Negative") %>%
  dplyr::select(-dplyr::contains("CPTAC4")) %>% # removes QC columns
  tidyr::pivot_longer(cols = starts_with('BEAT_AML'), names_to = 'sample', values_to = 'value') %>%
  tidyr::separate(sample, into = c('Beat','AML','PNL','SampleID.abbrev'), 
                  sep = '_') |>
  dplyr::mutate(SampleID = as.numeric(SampleID.abbrev)) |>
  dplyr::left_join(ba_meta, by = 'SampleID') 

ba_lneg <- readxl::read_xlsx(synapser::synGet('syn71210151')$path,
#                             here("BeatAMLPilot_Integration", "data", "processed", "beataml_lipids_log2_sum.xlsx"),
                             sheet = "Negative") %>%
  dplyr::select(-dplyr::contains("CPTAC4"))

ba_lpos <- readxl::read_xlsx(synapser::synGet('syn71210151')$path,
                             #here("BeatAMLPilot_Integration", "data", "processed", "beataml_lipids_log2_sum.xlsx"),
                                 sheet = "Positive") %>%
  dplyr::select(-dplyr::contains("CPTAC4"))


## Pilot Data -------------------------------------------------------


pi_lneg <- readxl::read_xlsx(synapser::synGet('syn71210151')$path,
                             #here("BeatAMLPilot_Integration", "data", "processed", "ptrc_lipids_log2_sum.xlsx"),
                             sheet = 'Negative')

pi_lpos <- readxl::read_xlsx(synapser::synGet('syn71210151')$path,
                             #here("BeatAMLPilot_Integration", "data", "processed", "ptrc_lipids_log2_sum.xlsx"),
                             sheet = 'Positive')

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


# # Update beatpilot_metadata with sample IDs -------------------------------------------------------------------------
# 
# beatpilot_metadata <- beatpilot_metadata %>%
#   dplyr::left_join(ba_meta %>%
#                      dplyr::select(Accession, SampleID) %>%
#                      dplyr::mutate(SampleID = as.character(SampleID))) %>%
#   dplyr::left_join(pi_metadat %>%
#                      dplyr::select(Accession, SampleID), by = "Accession")

### Data Processing ---------------------------------------------------------

## BeatAML Data -------------------------------------------------------------------------

# pmartR formatting -------------------------------------------------------

# Create initial edata, emeta, and fdata objects
ba_lneg_edat <- ba_lneg %>%
  dplyr::select(formatted_name, 
                dplyr::contains("BEAT_AML_PNL_"))
ba_lneg_emet <- ba_lneg %>%
  dplyr::select(-dplyr::contains("BEAT_AML_PNL_"))

ba_lneg_fdat <- data.frame(SampleID = colnames(ba_lneg_edat)[-1]) %>%
  dplyr::mutate(link_id = purrr::map_dbl(SampleID, function(x){
    as.numeric(strsplit(x, "_")[[1]][4]) # gets the sample number at the end of SampleID string
  })) %>%
  dplyr::left_join(ba_meta %>%
                     dplyr::rename(link_id = SampleID)) %>%
  dplyr::left_join(beatpilot_metadata %>% 
                     dplyr::left_join(ba_meta %>%
                                        dplyr::select(Accession, SampleID) %>%
                                        dplyr::rename(link_id = SampleID))) %>%
  dplyr::mutate(study = ifelse(is.na(study), "BeatAML", study)) %>%
  # For now, select subset of variables that are in common with pilot fdata
  dplyr::select(SampleID, Age, Age2, Sex, Sex2, Race, Race2, study) %>%
  dplyr::mutate(Race_mod = dplyr::case_when(
    Race == "Unknown" ~ NA,
    Race == "NA" ~ NA,
    Race %in% c("Black", "HispNative", "Asian") ~ "Non-white",
    TRUE ~ Race
  ), Race_mod = factor(Race_mod, levels = c("White", "Non-white"))) %>%
  dplyr::distinct()

# Check consistency of edat and fdat names
all(ba_lneg_fdat$SampleID %in% colnames(ba_lneg_edat)[-1])

ba_pm_lneg <- as.metabData(e_data = ba_lneg_edat,
                           f_data = ba_lneg_fdat,
                           e_meta = ba_lneg_emet,
                           edata_cname = "formatted_name", 
                           fdata_cname = "SampleID",
                           emeta_cname = "formatted_name",
                           data_scale = "log2",
                           data_types = "Negative")


# Create initial edata, emeta, and fdata objects
ba_lpos_edat <- ba_lpos %>%
  dplyr::select(formatted_name, 
                dplyr::contains("BEAT_AML_PNL_"))
ba_lpos_emet <- ba_lpos %>%
  dplyr::select(-dplyr::contains("BEAT_AML_PNL_"))

ba_lpos_fdat <- data.frame(SampleID = colnames(ba_lpos_edat)[-1]) %>%
  dplyr::mutate(link_id = purrr::map_dbl(SampleID, function(x){
    as.numeric(strsplit(x, "_")[[1]][4]) # gets the sample number at the end of SampleID string
  })) %>%
  dplyr::left_join(ba_meta %>%
                     dplyr::rename(link_id = SampleID)) %>%
  dplyr::left_join(beatpilot_metadata %>% 
                     dplyr::left_join(ba_meta %>%
                                        dplyr::select(Accession, SampleID) %>%
                                        dplyr::rename(link_id = SampleID))) %>%
  dplyr::mutate(study = ifelse(is.na(study), "BeatAML", study)) %>%
  # For now, select subset of variables that are in common with pilot fdata
  dplyr::select(SampleID, Age, Age2, Sex, Sex2, Race, Race2, study) %>%
  dplyr::mutate(Race_mod = dplyr::case_when(
    Race == "Unknown" ~ NA,
    Race == "NA" ~ NA,
    Race %in% c("Black", "HispNative", "Asian") ~ "Non-white",
    TRUE ~ Race
  ), Race_mod = factor(Race_mod, levels = c("White", "Non-white")))  %>%
  dplyr::distinct()

# Check consistency of edat and fdat names
all(ba_lpos_fdat$SampleID %in% colnames(ba_lpos_edat)[-1])

ba_pm_lpos <- as.metabData(e_data = ba_lpos_edat,
                           f_data = ba_lpos_fdat,
                           e_meta = ba_lpos_emet,
                           edata_cname = "formatted_name", 
                           fdata_cname = "SampleID",
                           emeta_cname = "formatted_name",
                           data_scale = "log2",
                           data_types = "Negative")

# -------------------------------------------------------------------------

## Pilot Data -------------------------------------------------------------------------

# pmartR formatting -------------------------------------------------------

# Create initial edata, emeta, and fdata objects
pi_lneg_edat <- pi_lneg %>%
  dplyr::select(-dplyr::contains("Lip_PB")) %>%
  dplyr::select(-dplyr::contains("Blank")) %>%
  dplyr::select(-dplyr::contains("QC_Pool")) %>%
  dplyr::select(formatted_name, dplyr::contains("Lip"))
pi_lneg_emet <- pi_lneg %>%
  dplyr::select(formatted_name:Lip_PB_03, dplyr::contains("QC_Pool"))

pi_lneg_fdat <- data.frame(SampleID = colnames(pi_lneg_edat)[-1]) %>%
  dplyr::mutate(link_id = purrr::map_dbl(SampleID, function(x){
    as.numeric(strsplit(x, "_")[[1]][2]) # gets the sample number at the end of SampleID string
  })) %>%
  dplyr::mutate(study = "pilot") %>%
  dplyr::left_join(pi_metadat %>%
                     dplyr::mutate(link_id = purrr::map_dbl(`Lipidomics ID`, function(x){
                       as.numeric(strsplit(x, "_")[[1]][4]) # gets the sample number at the end of SampleID string
                     })) %>%
                     dplyr::select(-SampleID)) %>%
  dplyr::left_join(beatpilot_metadata %>% 
                     dplyr::left_join(pi_metadat %>%
                                        dplyr::mutate(link_id = purrr::map_dbl(`Lipidomics ID`, function(x){
                                          as.numeric(strsplit(x, "_")[[1]][4]) # gets the sample number at the end of SampleID string
                                        })) %>%
                                        dplyr::select(Accession, link_id))) %>%
  # For now, select subset of variables that are in common with BeatAML fdata
  dplyr::select(SampleID, Age, Age2, Sex, Sex2, Race, Race2, study) %>%
  dplyr::mutate(Race_mod = dplyr::case_when(
    Race == "Unknown" ~ NA,
    Race == "NA" ~ NA,
    Race %in% c("Black", "HispNative", "Asian") ~ "Non-white",
    TRUE ~ Race
  ), Race_mod = factor(Race_mod, levels = c("White", "Non-white"))) %>%
  dplyr::distinct()

# Check consistency of edat and fdat names
all(pi_lneg_fdat$SampleID %in% colnames(pi_lneg_edat)[-1])

pi_pm_lneg <- as.metabData(e_data = pi_lneg_edat,
                           f_data = pi_lneg_fdat,
                           e_meta = pi_lneg_emet,
                           edata_cname = "formatted_name", 
                           fdata_cname = "SampleID",
                           emeta_cname = "formatted_name",
                           data_scale = "log2",
                           data_types = "Negative")



# Create initial edata, emeta, and fdata objects
pi_lpos_edat <- pi_lpos %>%
  dplyr::select(-dplyr::contains("Lip_PB")) %>%
  dplyr::select(-dplyr::contains("Blank")) %>%
  dplyr::select(-dplyr::contains("QC_Pool")) %>%
  dplyr::select(formatted_name, dplyr::contains("Lip"))
pi_lpos_emet <- pi_lpos %>%
  dplyr::select(formatted_name:Lip_PB_03, dplyr::contains("QC_Pool"))

pi_lpos_fdat <- data.frame(SampleID = colnames(pi_lpos_edat)[-1]) %>%
  dplyr::mutate(link_id = purrr::map_dbl(SampleID, function(x){
    as.numeric(strsplit(x, "_")[[1]][2]) # gets the sample number at the end of SampleID string
  })) %>%
  dplyr::mutate(study = "pilot") %>%
  dplyr::left_join(pi_metadat %>%
                     dplyr::mutate(link_id = purrr::map_dbl(`Lipidomics ID`, function(x){
                       as.numeric(strsplit(x, "_")[[1]][4]) # gets the sample number at the end of SampleID string
                     })) %>%
                     dplyr::select(-SampleID)) %>%
  dplyr::left_join(beatpilot_metadata %>% 
                     dplyr::left_join(pi_metadat %>%
                                        dplyr::mutate(link_id = purrr::map_dbl(`Lipidomics ID`, function(x){
                                          as.numeric(strsplit(x, "_")[[1]][4]) # gets the sample number at the end of SampleID string
                                        })) %>%
                                        dplyr::select(Accession, link_id))) %>%
  # For now, select subset of variables that are in common with BeatAML fdata
  dplyr::select(SampleID, Age, Age2, Sex, Sex2, Race, Race2, study) %>%
  dplyr::mutate(Race_mod = dplyr::case_when(
    Race == "Unknown" ~ NA,
    Race == "NA" ~ NA,
    Race %in% c("Black", "HispNative", "Asian") ~ "Non-white",
    TRUE ~ Race
  ), Race_mod = factor(Race_mod, levels = c("White", "Non-white"))) %>%
  dplyr::distinct()

# Check consistency of edat and fdat names
all(pi_lpos_fdat$SampleID %in% colnames(pi_lpos_edat)[-1])

pi_pm_lpos <- as.metabData(e_data = pi_lpos_edat,
                           f_data = pi_lpos_fdat,
                           e_meta = pi_lpos_emet,
                           edata_cname = "formatted_name", 
                           fdata_cname = "SampleID",
                           emeta_cname = "formatted_name",
                           data_scale = "log2",
                           data_types = "Positive")

# -------------------------------------------------------------------------


### Flagging Common Lipids ---------------------------------------------

# Positive ----------------------------------------------------------

pos_pi_match_cands <- vector("list", length = nrow(pi_pm_lpos$e_meta))
for(i in 1:nrow(pi_pm_lpos$e_meta)){
  
  if(is.na(pi_pm_lpos$e_meta$BEAT.AML.Alignment.ID[i])){
    next
  }
  
  piba_matched_alignid <- strsplit(pi_pm_lpos$e_meta$BEAT.AML.Alignment.ID[i], split = ",")[[1]]
  piba_matched_alignid <- gsub(" - not 100% sure", "", piba_matched_alignid)
  piba_matched_alignid <- gsub("\\?", "", piba_matched_alignid)
  piba_matched_alignid <- trimws(piba_matched_alignid)
  
  cands <- vector("list", length = nrow(ba_pm_lpos$e_meta))
  for(j in 1:nrow(ba_pm_lpos$e_meta)){
    ba_alignid <- strsplit(ba_pm_lpos$e_meta$Alignment.ID[j], split = ",")[[1]]
    ba_alignid <- trimws(ba_alignid)
    if(any(ba_alignid %in% piba_matched_alignid)){
      cands[[j]] <- data.frame(pi_pm_lpos$e_meta[i,]) %>%
        dplyr::select(formatted_name:BEAT.AML.RT) %>%
        dplyr::mutate(formatted_name_beat = ba_pm_lpos$e_meta$formatted_name[j],
                      Alignment.ID_beat = ba_pm_lpos$e_meta$Alignment.ID[j],
                      Average.Rt.min._beat = ba_pm_lpos$e_meta$Average.Rt.min.[j],
                      Average.Mz_beat = ba_pm_lpos$e_meta$Average.Mz[j],
                      original_name_beat = ba_pm_lpos$e_meta$original_name[j])
      
    }
  }
  cands <- Reduce("rbind", cands[!sapply(cands, is.null)])
  
  
  pos_pi_match_cands[[i]] <- cands
  
}
pos_pi_match_cands <- pos_pi_match_cands[!sapply(pos_pi_match_cands, is.null)]
dim_vec <- Reduce("c", lapply(pos_pi_match_cands, nrow))
table(dim_vec) # all match for match - no multi-matching. Good
pos_pi_match_cands_df <- Reduce("rbind", pos_pi_match_cands)



any(duplicated(pos_pi_match_cands_df$formatted_name_beat))
pos_dups <- pos_pi_match_cands_df$formatted_name_beat[which(duplicated(pos_pi_match_cands_df$formatted_name_beat))]
# View(pos_pi_match_cands_df %>%
#        dplyr::filter(formatted_name_beat %in% pos_dups))
# There appears to be a mistake in one of the matches that JK made. Will need to
# confirm with her at some point. For now, we will assume it is a mistake and correct 
# for it. 
pos_pi_match_cands_df <- pos_pi_match_cands_df %>%
  dplyr::filter(!(Alignment.ID == "4083"))



# Negative ----------------------------------------------------------

neg_align_data <- readxl::read_xlsx(synapser::synGet('syn71210151')$path)
  #here("BeatAMLPilot_Integration", "data", "PTRC_lipids_NEG_aligned_with_BEAT.xlsx"))
colnames(neg_align_data) <- make.names(colnames(neg_align_data))

pi_pm_lneg$e_meta <- pi_pm_lneg$e_meta %>%
  dplyr::left_join(neg_align_data %>%
                     dplyr::select(formatted_name, BEAT.NEG.Alignment.ID, BEAT.RT)) %>%
  dplyr::relocate(BEAT.NEG.Alignment.ID, BEAT.RT, .after = original_name) %>%
  # Some alignment IDs appear to be retention times instead. This is corrected based on
  # manual inspection
  dplyr::mutate(BEAT.NEG.Alignment.ID =  dplyr::case_when(
    BEAT.NEG.Alignment.ID == "19.443000000000001" ~ "10759",
    BEAT.NEG.Alignment.ID == "19.826000000000001" ~ "10756",
    TRUE ~ BEAT.NEG.Alignment.ID
  ))


neg_pi_match_cands <- vector("list", length = nrow(pi_pm_lneg$e_meta))
for(i in 1:nrow(pi_pm_lneg$e_meta)){
  
  if(is.na(pi_pm_lneg$e_meta$BEAT.NEG.Alignment.ID[i])){
    next
  }
  
  piba_matched_alignid <- strsplit(pi_pm_lneg$e_meta$BEAT.NEG.Alignment.ID[i], split = ",")[[1]]
  piba_matched_alignid <- gsub(" - not 100% sure", "", piba_matched_alignid)
  piba_matched_alignid <- gsub("\\?", "", piba_matched_alignid)
  piba_matched_alignid <- trimws(piba_matched_alignid)
  
  cands <- vector("list", length = nrow(ba_pm_lneg$e_meta))
  for(j in 1:nrow(ba_pm_lneg$e_meta)){
    ba_alignid <- strsplit(ba_pm_lneg$e_meta$Alignment.ID[j], split = ",")[[1]]
    ba_alignid <- trimws(ba_alignid)
    if(any(ba_alignid %in% piba_matched_alignid)){
      cands[[j]] <- data.frame(pi_pm_lneg$e_meta[i,]) %>%
        dplyr::select(formatted_name:BEAT.RT) %>%
        dplyr::mutate(formatted_name_beat = ba_pm_lneg$e_meta$formatted_name[j],
                      Alignment.ID_beat = ba_pm_lneg$e_meta$Alignment.ID[j],
                      Average.Rt.min._beat = ba_pm_lneg$e_meta$Average.Rt.min.[j],
                      Average.Mz_beat = ba_pm_lneg$e_meta$Average.Mz[j],
                      original_name_beat = ba_pm_lneg$e_meta$original_name[j])
      
    }
  }
  cands <- Reduce("rbind", cands[!sapply(cands, is.null)])
  
  
  neg_pi_match_cands[[i]] <- cands
  
}
neg_pi_match_cands <- neg_pi_match_cands[!sapply(neg_pi_match_cands, is.null)]
dim_vec <- Reduce("c", lapply(neg_pi_match_cands, nrow))
table(dim_vec) # all match for match - no multi-matching. Good
neg_pi_match_cands_df <- Reduce("rbind", neg_pi_match_cands)



any(duplicated(neg_pi_match_cands_df$formatted_name_beat))
neg_dups <- neg_pi_match_cands_df$formatted_name_beat[which(duplicated(neg_pi_match_cands_df$formatted_name_beat))]
# View(neg_pi_match_cands_df %>%
#        dplyr::filter(formatted_name_beat %in% neg_dups))
# No duplicates

# -------------------------------------------------------------------------


### Normalization + ComBat Only ----------------------------------------------------------

## Normalization -----------------------------------------------------------
# BeatAML -------------------------------------------------------------------

# Normalize data using median centering
ba_pm_lneg_norm <- normalize_global(omicsData = ba_pm_lneg,
                                    subset_fn = "all",
                                    norm_fn = "median",
                                    apply_norm = TRUE,
                                    backtransform = TRUE)

# Normalize data using median centering
ba_pm_lpos_norm <- normalize_global(omicsData = ba_pm_lpos,
                                    subset_fn = "all",
                                    norm_fn = "median",
                                    apply_norm = TRUE,
                                    backtransform = TRUE)

# Combine modes

ba_pm_lneg_norm_temp <- as.lipidData(e_data = ba_pm_lneg_norm$e_data %>%
                                       dplyr::mutate(formatted_name = paste0("neg_", formatted_name)), 
                                     f_data = ba_pm_lneg_norm$f_data, 
                                     e_meta = ba_pm_lneg_norm$e_meta %>%
                                       dplyr::mutate(formatted_name = paste0("neg_", formatted_name)),
                                     edata_cname = "formatted_name", 
                                     fdata_cname = "SampleID",
                                     emeta_cname = "formatted_name",
                                     data_scale = "log2",
                                     is_normalized = TRUE,
                                     data_types = "Negative")

ba_pm_lpos_norm_temp <- as.lipidData(e_data = ba_pm_lpos_norm$e_data %>%
                                       dplyr::mutate(formatted_name = paste0("pos_", formatted_name)), 
                                     f_data = ba_pm_lpos_norm$f_data, 
                                     e_meta = ba_pm_lpos_norm$e_meta %>%
                                       dplyr::mutate(formatted_name = paste0("pos_", formatted_name)),
                                     edata_cname = "formatted_name", 
                                     fdata_cname = "SampleID",
                                     emeta_cname = "formatted_name",
                                     data_scale = "log2",
                                     is_normalized = TRUE,
                                     data_types = "Positive")

ba_pm_comb <- combine_lipidData(ba_pm_lneg_norm_temp, ba_pm_lpos_norm_temp,
                                retain_filters = TRUE)
rm(ba_pm_lneg_norm_temp, ba_pm_lpos_norm_temp)

# Pilot -------------------------------------------------------------------

# Normalize data using median centering
pi_pm_lneg_norm <- normalize_global(omicsData = pi_pm_lneg,
                                    subset_fn = "all",
                                    norm_fn = "median",
                                    apply_norm = TRUE,
                                    backtransform = TRUE)

# Normalize data using median centering
pi_pm_lpos_norm <- normalize_global(omicsData = pi_pm_lpos,
                                    subset_fn = "all",
                                    norm_fn = "median",
                                    apply_norm = TRUE,
                                    backtransform = TRUE)

# Combine modes

pi_pm_lneg_norm_temp <- as.lipidData(e_data = pi_pm_lneg_norm$e_data %>%
                                       dplyr::mutate(formatted_name = paste0("neg_", formatted_name)), 
                                     f_data = pi_pm_lneg_norm$f_data, 
                                     e_meta = pi_pm_lneg_norm$e_meta %>%
                                       dplyr::mutate(formatted_name = paste0("neg_", formatted_name)),
                                     edata_cname = "formatted_name", 
                                     fdata_cname = "SampleID",
                                     emeta_cname = "formatted_name",
                                     data_scale = "log2",
                                     is_normalized = TRUE,
                                     data_types = "Negative")

pi_pm_lpos_norm_temp <- as.lipidData(e_data = pi_pm_lpos_norm$e_data %>%
                                       dplyr::mutate(formatted_name = paste0("pos_", formatted_name)), 
                                     f_data = pi_pm_lpos_norm$f_data, 
                                     e_meta = pi_pm_lpos_norm$e_meta %>%
                                       dplyr::mutate(formatted_name = paste0("pos_", formatted_name)),
                                     edata_cname = "formatted_name", 
                                     fdata_cname = "SampleID",
                                     emeta_cname = "formatted_name",
                                     data_scale = "log2",
                                     is_normalized = TRUE,
                                     data_types = "Positive")

pi_pm_comb <- combine_lipidData(pi_pm_lneg_norm_temp, pi_pm_lpos_norm_temp,
                                retain_filters = TRUE)
rm(pi_pm_lneg_norm_temp, pi_pm_lpos_norm_temp)

## Integrate ---------------------------------------------------

# Define subsets containing only common lipids

# BeatAML
my_filt <- custom_filter(ba_pm_comb, 
                         e_data_keep = c(neg_pi_match_cands_df %>%
                                           dplyr::mutate(formatted_name_beat = paste0("neg_", formatted_name_beat)) %>%
                                           .$formatted_name_beat,
                                         pos_pi_match_cands_df %>%
                                           dplyr::mutate(formatted_name_beat = paste0("pos_", formatted_name_beat)) %>%
                                           .$formatted_name_beat))
ba_pm_comb <- applyFilt(filter_object = my_filt,
                        omicsData = ba_pm_comb)


# Pilot
my_filt <- custom_filter(pi_pm_comb, 
                         e_data_keep = c(neg_pi_match_cands_df %>%
                                           dplyr::mutate(formatted_name = paste0("neg_", formatted_name)) %>%
                                           .$formatted_name,
                                         pos_pi_match_cands_df %>%
                                           dplyr::mutate(formatted_name = paste0("pos_", formatted_name)) %>%
                                           .$formatted_name))
pi_pm_comb <- applyFilt(filter_object = my_filt,
                         omicsData = pi_pm_comb)


# Rename the formatted_name in Pilot to match that of BeatAML
# First neg mode
for(i in 1:nrow(neg_pi_match_cands_df)){
  tempname <- paste0("neg_", neg_pi_match_cands_df$formatted_name_beat[i])
  
  temp_repname <- paste0("neg_", neg_pi_match_cands_df$formatted_name[i])
  temp_alignid <- neg_pi_match_cands_df$Alignment.ID[i]
  repidx <- which(pi_pm_comb$e_meta$formatted_name == temp_repname &
                    pi_pm_comb$e_meta$Alignment.ID == temp_alignid)
  
  pi_pm_comb$e_meta$formatted_name[repidx] <- tempname
  pi_pm_comb$e_data$formatted_name[repidx] <- tempname
}
# Then pos mode
for(i in 1:nrow(pos_pi_match_cands_df)){
  tempname <- paste0("pos_", pos_pi_match_cands_df$formatted_name_beat[i])
  
  temp_repname <- paste0("pos_", pos_pi_match_cands_df$formatted_name[i])
  temp_alignid <- pos_pi_match_cands_df$Alignment.ID[i]
  repidx <- which(pi_pm_comb$e_meta$formatted_name == temp_repname &
                    pi_pm_comb$e_meta$Alignment.ID == temp_alignid)
  
  pi_pm_comb$e_meta$formatted_name[repidx] <- tempname
  pi_pm_comb$e_data$formatted_name[repidx] <- tempname
}
all(pi_pm_comb$e_data$formatted_name %in% ba_pm_comb$e_data$formatted_name) &
  all(ba_pm_comb$e_data$formatted_name %in% pi_pm_comb$e_data$formatted_name)


# Define new pmartR objects based on a concatenation of e_data, f_data
edat_common_comb <- ba_pm_comb$e_data %>%
  dplyr::left_join(pi_pm_comb$e_data, by = "formatted_name")

fdat_common_comb <- rbind.data.frame(ba_pm_comb$f_data %>%
                                       dplyr::mutate(batchid = 1,
                                                     study = "beataml"),
                                     pi_pm_comb$f_data %>%
                                       dplyr::mutate(batchid = 2,
                                                     study = "pilot"))

emeta_common_comb <- ba_pm_comb$e_meta %>%
  dplyr::left_join(pi_pm_comb$e_meta, by = "formatted_name") %>%
  setNames(gsub("\\.x", "\\.beataml", colnames(.))) %>%
  setNames(gsub("\\.y", "\\.pilot", colnames(.)))

# all(names(edat_common_comb)[-1] %in% fdat_common_comb$SampleID)
# setdiff(names(edat_common_comb)[-1], fdat_common_comb$SampleID)

pm_common_comb <- as.lipidData(e_data = edat_common_comb, 
                                f_data = fdat_common_comb,
                                e_meta = emeta_common_comb,
                                emeta_cname = "formatted_name",
                                edata_cname = "formatted_name", 
                                fdata_cname = "SampleID",
                                data_scale = "log2",
                                is_normalized = TRUE)

# Several samples (13 specifically: 12 Beat AML, 1 Pilot) have been removed due to
# missing race information
# pm_common_comb <- group_designation(pm_common_comb, main_effects = c("Race_mod"), batch_id = "batchid")
pm_common_comb <- group_designation(pm_common_comb, main_effects = c("study"), batch_id = "batchid")

pm_common_comb_combat <- bc_combat(omicsData = pm_common_comb, use_groups = FALSE)

# -------------------------------------------------------------------------


### Analyses -------------------------------------------------------------------

## RMD Outlier Assessment --------------------------------------------------

# BeatAML -----------------------------------------------------------------

ba_pm_comb <- group_designation(ba_pm_comb, main_effects = c("study"))

myrmd <- rmd_filter(ba_pm_comb)
png(here("BeatAMLPilot_Integration", "figures", "outlier_beataml_lipid_rmd.png"), units="in",
    width=10, height=7, res=300)
plot(myrmd, pvalue_threshold = 0.00001)
dev.off()

# Pilot -------------------------------------------------------------------

pi_pm_comb <- group_designation(pi_pm_comb, main_effects = c("study"))

myrmd <- rmd_filter(pi_pm_comb)
png(here("BeatAMLPilot_Integration", "figures", "outlier_pilot_lipid_rmd.png"), units="in",
    width=10, height=7, res=300)
plot(myrmd, pvalue_threshold = 0.00001)
dev.off()

# Combined, No BC ---------------------------------------------------------

pm_common_comb <- group_designation(pm_common_comb, main_effects = c("study"), batch_id = "batchid")

myrmd <- rmd_filter(pm_common_comb)
png(here("BeatAMLPilot_Integration", "figures", "outlier_nobc_lipid_rmd.png"), units="in",
    width=10, height=7, res=300)
plot(myrmd, pvalue_threshold = 0.00001)
dev.off()

# Combined, BC ------------------------------------------------------------

pm_common_comb_combat <- group_designation(pm_common_comb_combat, main_effects = c("study"), batch_id = "batchid")

myrmd <- rmd_filter(pm_common_comb_combat)
png(here("BeatAMLPilot_Integration", "figures", "outlier_combat_lipid_rmd.png"), units="in",
    width=10, height=7, res=300)
plot(myrmd, pvalue_threshold = 0.00001)
dev.off()

# -------------------------------------------------------------------------

## PCA --------------------------------------------------------------------

# BeatAML -----------------------------------------------------------------

ba_pm_comb <- group_designation(ba_pm_comb, main_effects = c("study"))

mypca <- dim_reduction(ba_pm_comb)
png(here("BeatAMLPilot_Integration", "figures", "beataml_lipid_pca.png"), units="in",
    width=10, height=7, res=300)
plot(mypca, pvalue_threshold = 0.00001, title_lab = "Lipids: BeatAML")
dev.off()

# Pilot -------------------------------------------------------------------

pi_pm_comb <- group_designation(pi_pm_comb, main_effects = c("study"))

mypca <- dim_reduction(pi_pm_comb)
png(here("BeatAMLPilot_Integration", "figures", "pilot_lipid_pca.png"), units="in",
    width=10, height=7, res=300)
plot(mypca, pvalue_threshold = 0.00001, title_lab = "Lipids: Pilot")
dev.off()

# Combined, No BC ---------------------------------------------------------

mypca <- dim_reduction(pm_common_comb)

png(here("BeatAMLPilot_Integration", "figures", "nobc_lipid_pca.png"), units="in",
    width=10, height=7, res=300)
plot(mypca, omicsData = pm_common_comb, color_by = "study", 
     title_lab = "Lipids: No Batch-Correction")+theme(text=element_text(size=21))
dev.off()

# Combined, BC ------------------------------------------------------------

mypca <- dim_reduction(pm_common_comb_combat)

png(here("BeatAMLPilot_Integration", "figures", "combat_lipid_pca.png"), units="in",
    width=10, height=7, res=300)
plot(mypca, omicsData = pm_common_comb_combat, color_by = "study", 
     title_lab = "Lipids: ComBat")+theme(text=element_text(size=21))
dev.off()

# -------------------------------------------------------------------------

## ANOVA -------------------------------------------------------------------

# BeatAML -----------------------------------------------------------------

ba_pm_comb_temp <- group_designation(ba_pm_comb, main_effects = c("Race_mod"))
gdf <- attr(ba_pm_comb_temp, "group_DF") %>% 
  dplyr::select(-SampleID) %>% # Remove replicate identifier; only retain design information
  dplyr::distinct()
gdf_comps <- data.frame(Control = c(gdf %>% 
                                      dplyr::filter(Group == "Non-white") %>% .$Group),
                        Test    = c(gdf %>% 
                                      dplyr::filter(Group == "White") %>% .$Group))
ba_hitest_race <- imd_anova(omicsData = ba_pm_comb_temp, 
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
  dplyr::select(formatted_name, `P_value_A_White_vs_Non-white`, `Flag_A_White_vs_Non-white`) %>%
  dplyr::filter(`Flag_A_White_vs_Non-white` != 0)

# ba_pm_comb_temp <- group_designation(ba_pm_comb, main_effects = c("study"))
# 
# ba_hitest_study <- imd_anova(omicsData = ba_pm_comb_temp, 
#                              test_method = "anova", 
#                              pval_adjust_a_multcomp = "holm", 
#                              pval_thresh = 0.05)
# summary(ba_hitest_study)

# Pilot -------------------------------------------------------------------

pi_pm_comb_temp <- group_designation(pi_pm_comb, main_effects = c("Race_mod"))
gdf <- attr(pi_pm_comb_temp, "group_DF") %>% 
  dplyr::select(-SampleID) %>% # Remove replicate identifier; only retain design information
  dplyr::distinct()
gdf_comps <- data.frame(Control = c(gdf %>% 
                                      dplyr::filter(Group == "Non-white") %>% .$Group),
                        Test    = c(gdf %>% 
                                      dplyr::filter(Group == "White") %>% .$Group))
pi_hitest_race <- imd_anova(omicsData = pi_pm_comb_temp,
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
  dplyr::select(formatted_name, `P_value_A_White_vs_Non-white`, `Flag_A_White_vs_Non-white`) %>%
  dplyr::filter(`Flag_A_White_vs_Non-white` != 0)

# pi_pm_comb_temp <- group_designation(pi_pm_comb, main_effects = c("study"))
# 
# pi_hitest_study <- imd_anova(omicsData = pi_pm_comb_temp, 
#                              test_method = "anova", 
#                              pval_adjust_a_multcomp = "holm", 
#                              pval_thresh = 0.05)
# summary(pi_hitest_study)

# Combined, No BC ---------------------------------------------------------

pm_common_comb_temp <- group_designation(pm_common_comb, main_effects = c("Race_mod"))
gdf <- attr(pm_common_comb_temp, "group_DF") %>% 
  dplyr::select(-SampleID) %>% # Remove replicate identifier; only retain design information
  dplyr::distinct()
gdf_comps <- data.frame(Control = c(gdf %>% 
                                      dplyr::filter(Group == "Non-white") %>% .$Group),
                        Test    = c(gdf %>% 
                                      dplyr::filter(Group == "White") %>% .$Group))
nobc_hitest_race <- imd_anova(omicsData = pm_common_comb_temp,
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
  dplyr::select(formatted_name, `P_value_A_White_vs_Non-white`, `Flag_A_White_vs_Non-white`) %>%
  dplyr::filter(`Flag_A_White_vs_Non-white` != 0)

pm_common_comb_temp <- group_designation(pm_common_comb, main_effects = c("study"))
gdf <- attr(pm_common_comb_temp, "group_DF") %>% 
  dplyr::select(-SampleID) %>% # Remove replicate identifier; only retain design information
  dplyr::distinct()
gdf_comps <- data.frame(Control = c(gdf %>% 
                                      dplyr::filter(Group == "beataml") %>% .$Group),
                        Test    = c(gdf %>% 
                                      dplyr::filter(Group == "pilot") %>% .$Group))
nobc_hitest_study <- imd_anova(omicsData = pm_common_comb_temp,
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
  dplyr::select(formatted_name, `P_value_A_pilot_vs_beataml`, `Flag_A_pilot_vs_beataml`) %>%
  dplyr::filter(`Flag_A_pilot_vs_beataml` != 0)

# Combined, BC ------------------------------------------------------------

pm_common_comb_combat_temp <- group_designation(pm_common_comb_combat, main_effects = c("Race_mod"))
gdf <- attr(pm_common_comb_combat_temp, "group_DF") %>% 
  dplyr::select(-SampleID) %>% # Remove replicate identifier; only retain design information
  dplyr::distinct()
gdf_comps <- data.frame(Control = c(gdf %>% 
                                      dplyr::filter(Group == "Non-white") %>% .$Group),
                        Test    = c(gdf %>% 
                                      dplyr::filter(Group == "White") %>% .$Group))
combat_hitest_race <- imd_anova(omicsData = pm_common_comb_combat_temp,
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
  dplyr::select(formatted_name, `P_value_A_White_vs_Non-white`, `Flag_A_White_vs_Non-white`) %>%
  dplyr::filter(`Flag_A_White_vs_Non-white` != 0)

pm_common_comb_combat_temp <- group_designation(pm_common_comb_combat, main_effects = c("study"))
gdf <- attr(pm_common_comb_combat_temp, "group_DF") %>% 
  dplyr::select(-SampleID) %>% # Remove replicate identifier; only retain design information
  dplyr::distinct()
gdf_comps <- data.frame(Control = c(gdf %>% 
                                      dplyr::filter(Group == "beataml") %>% .$Group),
                        Test    = c(gdf %>% 
                                      dplyr::filter(Group == "pilot") %>% .$Group))
combat_hitest_study <- imd_anova(omicsData = pm_common_comb_combat_temp,
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
  dplyr::select(formatted_name, `P_value_A_pilot_vs_beataml`, `Flag_A_pilot_vs_beataml`) %>%
  dplyr::filter(`Flag_A_pilot_vs_beataml` != 0)

# -------------------------------------------------------------------------

## Significant Set Comparisons ---------------------------------------------

# BeatAML vs Pilot --------------------------------------------------------

# Commonly significant, Race Comparison
intersect(whichsig_ba_hitest_race$formatted_name, whichsig_pi_hitest_race$formatted_name)

# Uniquely significant, BeatAML, Race Comparison
whichsig_ba_hitest_race %>%
  dplyr::filter(!(formatted_name %in% intersect(whichsig_ba_hitest_race$formatted_name, whichsig_pi_hitest_race$formatted_name))) %>%
  .$formatted_name

# Uniquely significant, Pilot, Race Comparison
whichsig_pi_hitest_race %>%
  dplyr::filter(!(formatted_name %in% intersect(whichsig_ba_hitest_race$formatted_name, whichsig_pi_hitest_race$formatted_name))) %>%
  .$formatted_name

# BeatAML vs Combined, No BC ----------------------------------------------

# Commonly significant, Race Comparison
intersect(whichsig_ba_hitest_race$formatted_name, whichsig_nobc_hitest_race$formatted_name)

# Uniquely significant, BeatAML, Race Comparison
whichsig_ba_hitest_race %>%
  dplyr::filter(!(formatted_name %in% intersect(whichsig_ba_hitest_race$formatted_name, whichsig_pi_hitest_race$formatted_name))) %>%
  .$formatted_name

# Uniquely significant, Combined No BC, Race Comparison
whichsig_nobc_hitest_race %>%
  dplyr::filter(!(formatted_name %in% intersect(whichsig_ba_hitest_race$formatted_name, whichsig_nobc_hitest_race$formatted_name))) %>%
  .$formatted_name

# BeatAML vs Combined, BC -------------------------------------------------

# Commonly significant, Race Comparison
intersect(whichsig_ba_hitest_race$formatted_name, whichsig_combat_hitest_race$formatted_name)

# Uniquely significant, BeatAML, Race Comparison
whichsig_ba_hitest_race %>%
  dplyr::filter(!(formatted_name %in% intersect(whichsig_ba_hitest_race$formatted_name, whichsig_combat_hitest_race$formatted_name))) %>%
  .$formatted_name

# Uniquely significant, Combined BC, Race Comparison
whichsig_combat_hitest_race %>%
  dplyr::filter(!(formatted_name %in% intersect(whichsig_ba_hitest_race$formatted_name, whichsig_combat_hitest_race$formatted_name))) %>%
  .$formatted_name

# Pilot vs Combined, No BC ----------------------------------------------

# Commonly significant, Race Comparison
intersect(whichsig_pi_hitest_race$formatted_name, whichsig_nobc_hitest_race$formatted_name)

# Uniquely significant, Pilot, Race Comparison
whichsig_pi_hitest_race %>%
  dplyr::filter(!(formatted_name %in% intersect(whichsig_pi_hitest_race$formatted_name, whichsig_pi_hitest_race$formatted_name))) %>%
  .$formatted_name

# Uniquely significant, Combined No BC, Race Comparison
whichsig_nobc_hitest_race %>%
  dplyr::filter(!(formatted_name %in% intersect(whichsig_pi_hitest_race$formatted_name, whichsig_nobc_hitest_race$formatted_name))) %>%
  .$formatted_name

# Pilot vs Combined, BC -------------------------------------------------

# Commonly significant, Race Comparison
intersect(whichsig_pi_hitest_race$formatted_name, whichsig_combat_hitest_race$formatted_name)

# Uniquely significant, Pilot, Race Comparison
whichsig_pi_hitest_race %>%
  dplyr::filter(!(formatted_name %in% intersect(whichsig_pi_hitest_race$formatted_name, whichsig_combat_hitest_race$formatted_name))) %>%
  .$formatted_name

# Uniquely significant, Combined BC, Race Comparison
whichsig_combat_hitest_race %>%
  dplyr::filter(!(formatted_name %in% intersect(whichsig_pi_hitest_race$formatted_name, whichsig_combat_hitest_race$formatted_name))) %>%
  .$formatted_name

# Combined, No BC vs Combined, BC -------------------------------------------------

# Commonly significant, Race Comparison
intersect(whichsig_nobc_hitest_race$formatted_name, whichsig_combat_hitest_race$formatted_name)

# Uniquely significant, Combined No BC, Race Comparison
whichsig_nobc_hitest_race %>%
  dplyr::filter(!(formatted_name %in% intersect(whichsig_nobc_hitest_race$formatted_name, whichsig_combat_hitest_race$formatted_name))) %>%
  .$formatted_name

# Uniquely significant, Combined BC, Race Comparison
whichsig_combat_hitest_race %>%
  dplyr::filter(!(formatted_name %in% intersect(whichsig_nobc_hitest_race$formatted_name, whichsig_combat_hitest_race$formatted_name))) %>%
  .$formatted_name

# -------------------------------------------------------------------------


