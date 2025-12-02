library(magrittr)
library(pmartR)
library(ggplot2)
library(synapser)
library(here)

i_am("src/r/lipids/beat_lipids_proc.R")

synapser::synLogin()

##### Data Processing ---------------------------------------------------------
#### Import ------------------------------------------------------------------

#data_filenames <- list.files(here("BeatAMLPilot_Integration","data"))

# Neg ---------------------------------------------------------------------

## Skip first row since contains fdata
lneg_dat <- readxl::read_excel(path = synapser::synGet('syn71718883')$path,
                               sheet = "NEG_with_QCs_v2")
colnames(lneg_dat) <- make.names(colnames(lneg_dat))

# Pos ---------------------------------------------------------------------

# Note: Had to use local datafile here because JK had modified it to include additional
# lipid classes. (See email from JK on 8/14/2025 with subject line "PTRC aligned with BEAT - POS")
lpos_dat <- readxl::read_excel(path = synapser::synGet('syn71718883')$path,
                                 # here("BeatAMLPilot_Integration", "data", 
                                      #     data_filenames[grepl("BEAT_AML_", data_filenames, ignore.case = TRUE)]),
                               sheet = "POS_with_QCs")
colnames(lpos_dat) <- make.names(colnames(lpos_dat))

# -------------------------------------------------------------------------

#### Filtering ---------------------------------------------------------------
### S/N Filter (insufficient info to do for these data) ---------------------------------------------------------------------
# Neg ---------------------------------------------------------------------

# filt_count_neg <- nrow(lneg_dat %>%
#                          dplyr::filter(S.N.average < 10))
# 
# lneg_dat <- lneg_dat %>%
#   dplyr::filter(S.N.average >= 10)

# Pos ---------------------------------------------------------------------

# filt_count_pos <- nrow(lpos_dat %>%
#                          dplyr::filter(S.N.average < 10))
# 
# lpos_dat <- lpos_dat %>%
#   dplyr::filter(S.N.average >= 10)

# -------------------------------------------------------------------------

### Process Blank / Instrument Blank Filter (insufficient info to do for these data) ---------------------------------
# Neg ---------------------------------------------------------------------

# Filter out if average process blank value for lipid is greater than 10% or more of the
# values in the samples, remove it.
# Compare between using the average of the process blanks to the average of the instrument
# blanks
# pb_means <- apply(lneg_dat %>% 
#                     dplyr::select(dplyr::contains("PB")), 1, function(x){mean(as.numeric(x), na.rm = TRUE)})
# ib_means <- apply(lneg_dat %>% 
#                     dplyr::select(dplyr::contains("Blank_Co")), 1, function(x){mean(as.numeric(x), na.rm = TRUE)})
# blank_comparison <- vector("list", length = length(pb_means))
# for(i in 1:length(pb_means)){
#   tempdat <- as.numeric(lneg_dat %>%
#                           dplyr::select(dplyr::contains("PTRC_exp26_Lip"),
#                                         -dplyr::contains("PB"),
#                                         -dplyr::contains("Blank_Co")) %>%
#                           .[i,])
#   
#   pb_prop <- sum(pb_means[i] > tempdat)/length(tempdat)*100
#   ib_prop <- sum(ib_means[i] > tempdat)/length(tempdat)*100
#   blank_comparison[[i]] <- data.frame(Metabolite.name = lneg_dat$Metabolite.name[i],
#                                       Average.Rt.min. = lneg_dat$Average.Rt.min.[i],
#                                       Adduct.type = lneg_dat$Adduct.type[i],
#                                       pb_prop = pb_prop,
#                                       ib_prop = ib_prop,
#                                       pb_flag = ifelse(pb_prop >= 10, 1, 0),
#                                       ib_flag = ifelse(ib_prop >= 10, 1, 0))
# }
# blank_comparison <- Reduce("rbind", blank_comparison) %>%
#   dplyr::mutate(join_flag = ifelse(pb_flag == 1 & ib_flag == 1,1,0),
#                 different_flag = ifelse(pb_flag == ib_flag, 0, 1))
# 
# blank_comparison_neg <- blank_comparison %>%
#   dplyr::mutate(blank_remove = ifelse(pb_flag == 1 | ib_flag == 1, 1, 0)) %>%
#   dplyr::select(-pb_prop, -ib_prop)
# 
# blank_comparison <- blank_comparison %>%
#   dplyr::mutate(blank_remove = ifelse(pb_flag == 1 | ib_flag == 1, 1, 0)) %>%
#   dplyr::select(Metabolite.name, Average.Rt.min., Adduct.type, blank_remove)
# 
# # As per JK 4/25/24...remove based on both pb and ib
# lneg_dat <- lneg_dat %>%
#   dplyr::left_join(blank_comparison) %>%
#   dplyr::filter(blank_remove == 0) %>%
#   dplyr::select(-blank_remove)
# nrow(blank_comparison_neg) - nrow(lneg_dat)
# 
# rm(i, ib_means, ib_prop, pb_means, pb_prop, tempdat, blank_comparison)

# Pos ---------------------------------------------------------------------

# Filter out if average process blank value for lipid is greater than 10% or more of the
# values in the samples, remove it.
# Compare between using the average of the process blanks to the average of the instrument
# blanks
# pb_means <- apply(lpos_dat %>% 
#                     dplyr::select(dplyr::contains("PB")), 1, function(x){mean(as.numeric(x), na.rm = TRUE)})
# ib_means <- apply(lpos_dat %>% 
#                     dplyr::select(dplyr::contains("Blank_Co")), 1, function(x){mean(as.numeric(x), na.rm = TRUE)})
# blank_comparison <- vector("list", length = length(pb_means))
# for(i in 1:length(pb_means)){
#   tempdat <- as.numeric(lpos_dat %>%
#                           dplyr::select(dplyr::contains("PTRC_exp26_Lip"),
#                                         -dplyr::contains("PB"),
#                                         -dplyr::contains("Blank_Co")) %>%
#                           .[i,])
#   pb_prop <- sum(pb_means[i] > tempdat)/length(tempdat)*100
#   ib_prop <- sum(ib_means[i] > tempdat)/length(tempdat)*100
#   blank_comparison[[i]] <- data.frame(Metabolite.name = lpos_dat$Metabolite.name[i],
#                                       Average.Rt.min. = lpos_dat$Average.Rt.min.[i],
#                                       Adduct.type = lpos_dat$Adduct.type[i],
#                                       pb_prop = pb_prop,
#                                       ib_prop = ib_prop,
#                                       pb_flag = ifelse(pb_prop >= 10, 1, 0),
#                                       ib_flag = ifelse(ib_prop >= 10, 1, 0))
# }
# blank_comparison <- Reduce("rbind", blank_comparison) %>%
#   dplyr::mutate(join_flag = ifelse(pb_flag == 1 & ib_flag == 1,1,0),
#                 different_flag = ifelse(pb_flag == ib_flag, 0, 1))
# 
# blank_comparison_pos <- blank_comparison %>%
#   dplyr::mutate(blank_remove = ifelse(pb_flag == 1 | ib_flag == 1, 1, 0)) %>%
#   dplyr::select(-pb_prop, -ib_prop)
# 
# blank_comparison <- blank_comparison %>%
#   dplyr::mutate(blank_remove = ifelse(pb_flag == 1 | ib_flag == 1, 1, 0)) %>%
#   dplyr::select(Metabolite.name, Average.Rt.min., Adduct.type, blank_remove)
# 
# # As per JK 4/25/24...remove based on both pb and ib
# lpos_dat <- lpos_dat %>%
#   dplyr::left_join(blank_comparison) %>%
#   dplyr::filter(blank_remove == 0) %>%
#   dplyr::select(-blank_remove)
# nrow(blank_comparison_pos) - nrow(lpos_dat)
# 
# rm(i, ib_means, ib_prop, pb_means, pb_prop, tempdat, blank_comparison)
# 
# writexl::write_xlsx(x = list(`Neg` = blank_comparison_neg,
#                              `Pos` = blank_comparison_pos),
#                     path = here("data", "processed", "ptrc_lipids_blank_flags.xlsx"))

# -------------------------------------------------------------------------

# rm(blank_comparison_neg, blank_comparison_pos)

#### Processing --------------------------------------------------------------
### Duplicate Labeling ------------------------------------------------------
# Neg ---------------------------------------------------------------------

# Check for duplicate metabolites (in name and adduct), and append _A, _B, _C, etc labels as needed
lneg_dat <- lneg_dat %>%
  # Create new metabolite name to include adduct information
  dplyr::mutate(Metabolite.name = paste0("neg_", Metabolite.name),
                metname.with.adduct = paste0(Metabolite.name, "_", Adduct.type)) %>%
  dplyr::group_by(metname.with.adduct) %>%
  dplyr::arrange(Average.Rt.min.) %>%
  dplyr::mutate(dupcount = dplyr::n(), # count total duplicates for a given metname.with.adduct
                dupnum = dplyr::row_number(), # assign numeric labels for each duplicate
                dupletter = toupper(letters[dupnum]), # convert numeric labels to capital letters
                # Rename metabolites by appending duplicate labels where appropriate
                # Metabolite.name = ifelse(dupcount > 1, paste0(Metabolite.name, "__", dupletter), Metabolite.name),
                metname.with.adduct = ifelse(dupcount > 1, paste0(metname.with.adduct, "__", dupletter), metname.with.adduct)) %>%
  dplyr::ungroup() %>%
  dplyr::relocate(metname.with.adduct, .after = `Adduct.type`) %>%
  dplyr::relocate(dupcount, dupnum, dupletter, .after = `MS.MS.spectrum`)

neg_dup_df <- lneg_dat %>% 
  dplyr::filter(dupcount > 1) %>% 
  dplyr::select(Average.Rt.min.:metname.with.adduct) %>% 
  dplyr::arrange(metname.with.adduct)

# Pos ---------------------------------------------------------------------

# Check for duplicate metabolites (in name and adduct), and append _A, _B, _C, etc labels as needed
lpos_dat <- lpos_dat %>%
  # Create new metabolite name to include adduct information
  dplyr::mutate(Metabolite.name = paste0("pos_", Metabolite.name),
                metname.with.adduct = paste0(Metabolite.name, "_", Adduct.type)) %>%
  dplyr::group_by(metname.with.adduct) %>%
  dplyr::arrange(Average.Rt.min.) %>%
  dplyr::mutate(dupcount = dplyr::n(), # count total duplicates for a given metname.with.adduct
                dupnum = dplyr::row_number(), # assign numeric labels for each duplicate
                dupletter = toupper(letters[dupnum]), # convert numeric labels to capital letters
                # Rename metabolites by appending duplicate labels where appropriate
                # Metabolite.name = ifelse(dupcount > 1, paste0(Metabolite.name, "__", dupletter), Metabolite.name),
                metname.with.adduct = ifelse(dupcount > 1, paste0(metname.with.adduct, "__", dupletter), metname.with.adduct)) %>%
  dplyr::ungroup() %>%
  dplyr::relocate(metname.with.adduct, .after = `Adduct.type`) %>%
  dplyr::relocate(dupcount, dupnum, dupletter, .after = `MS.MS.spectrum`)

pos_dup_df <- lpos_dat %>% 
  dplyr::filter(dupcount > 1) %>% 
  dplyr::select(Average.Rt.min.:metname.with.adduct) %>% 
  dplyr::arrange(metname.with.adduct)

# -------------------------------------------------------------------------
### pmartR Component Creation --------------------------------------------------

# Neg ---------------------------------------------------------------------

# Create initial edata, emeta, and fdata objects
lneg_edat <- lneg_dat %>%
  dplyr::select(metname.with.adduct,
                dplyr::contains("_Neg_"))
lneg_emet <- lneg_dat %>%
  dplyr::select(-dplyr::contains("_Neg_"))

lneg_fdat <- data.frame(SampleID = colnames(lneg_edat)[-1]) %>%
  dplyr::mutate(Category = dplyr::case_when(
    grepl("BEAT_AML", SampleID) ~ "Sample",
    grepl("CPTAC4_AML_BM", SampleID) ~ "QC: CPTAC4_AML_BM",
    grepl("CPTAC4_AML_WB", SampleID) ~ "QC: CPTAC4_AML_WB",
    TRUE ~ ""
  )) %>%
  dplyr::mutate(SampleID_og = SampleID,
                SampleID = make.names(SampleID)) %>%
  dplyr::filter(!is.na(Category))

# Check consistency of edat and fdat names
all(lneg_fdat$SampleID %in% colnames(lneg_edat)[-1])

# Pos ---------------------------------------------------------------------

# Create initial edata, emeta, and fdata objects
lpos_edat <- lpos_dat %>%
  dplyr::select(metname.with.adduct,
                dplyr::contains("_Pos_"))
lpos_emet <- lpos_dat %>%
  dplyr::select(-dplyr::contains("_Pos_"))

lpos_fdat <- data.frame(SampleID = colnames(lpos_edat)[-1]) %>%
  dplyr::mutate(Category = dplyr::case_when(
    grepl("BEAT_AML", SampleID) ~ "Sample",
    grepl("CPTAC4_AML_BM", SampleID) ~ "QC: CPTAC4_AML_BM",
    grepl("CPTAC4_AML_WB", SampleID) ~ "QC: CPTAC4_AML_WB",
    TRUE ~ ""
  )) %>%
  dplyr::mutate(SampleID_og = SampleID,
                SampleID = make.names(SampleID)) %>%
  dplyr::filter(!is.na(Category))

# Check consistency of edat and fdat names
all(lpos_fdat$SampleID %in% colnames(lpos_edat)[-1])

# -------------------------------------------------------------------------
### Differing Adducts Among Duplicates --------------------------------------
## Isolation ---------------------------------------------------------------
# Neg ---------------------------------------------------------------------

# Isolate cases where there are duplicated metabolite names with different adducts.
lneg_tempdat <- lneg_dat %>%
  dplyr::group_by(Metabolite.name) %>%
  dplyr::mutate(metname_dupcount = dplyr::n()) %>%
  dplyr::relocate(metname_dupcount, dupcount) %>%
  dplyr::filter(metname_dupcount > dupcount) %>%
  dplyr::ungroup()

# temp <- lneg_tempdat %>%
#   dplyr::mutate(round_Average.Rt.min. = round(Average.Rt.min., 2),
#                 round_Average.Rt.min. = dplyr::case_when(
#                   metname.with.adduct == "neg_CL 72:7|CL 18:1_18:2_18:2_18:2_[M-2H]2-" ~ 27.35,
#                   metname.with.adduct == "neg_Cer 32:0;O2|Cer 18:0;O2/14:0_[M-H]-" ~ 22.33,
#                   metname.with.adduct == "neg_Cer 40:1;O2|Cer 18:0;O2/22:1_[M+CH3COO]-" ~ 26.64,
#                   TRUE ~ round_Average.Rt.min.
#                 ),
#                 metname.with.rt = paste0(Metabolite.name, "_", round_Average.Rt.min.)) %>%
#   dplyr::arrange(Metabolite.name, round_Average.Rt.min.) %>%
#   dplyr::relocate(Average.Rt.min., round_Average.Rt.min., .after = metname.with.adduct)

# Note: Only true duplicates if the retention times are similar enough.
# After the summation option below, treat the unpaired cases with unique retention times
# as non-duplicates. 
# The setpair IDs below are used to denote which duplicates should be paired.
lneg_tempdat <- lneg_tempdat %>%
  dplyr::mutate(round_Average.Rt.min. = round(Average.Rt.min., 2),
                round_Average.Rt.min. = dplyr::case_when(
                  metname.with.adduct == "neg_CL 72:7|CL 18:1_18:2_18:2_18:2_[M-2H]2-" ~ 27.35,
                  metname.with.adduct == "neg_Cer 32:0;O2|Cer 18:0;O2/14:0_[M-H]-" ~ 22.33,
                  metname.with.adduct == "neg_Cer 40:1;O2|Cer 18:0;O2/22:1_[M+CH3COO]-" ~ 26.64,
                  TRUE ~ round_Average.Rt.min.
                ),
                metname.with.rt = paste0(Metabolite.name, "_", round_Average.Rt.min.)) %>%
  dplyr::group_by(metname.with.rt) %>%
  dplyr::arrange(round_Average.Rt.min.) %>% 
  dplyr::mutate(rtdupcount = dplyr::n(), # count total duplicates for a given metname.with.rt
                setpair = ifelse(rtdupcount > 1, "AA", "")) %>%
  dplyr::ungroup() %>%
  dplyr::select(-rtdupcount, -round_Average.Rt.min., -metname.with.rt) %>%
  dplyr::mutate(metname.with.setpair = ifelse(setpair != "", paste0(Metabolite.name, "__", setpair), metname.with.adduct)) %>%
  dplyr::group_by(metname.with.setpair) %>%
  dplyr::arrange(Average.Rt.min.) %>%
  dplyr::mutate(dupcount2 = dplyr::n(), # count total duplicates for a given metname.with.setpair
                dupnum2 = ceiling(dplyr::row_number()/2), # assign numeric labels for each duplicate
                dupletter2 = toupper(letters[dupnum2]), # convert numeric labels to capital letters
                # Rename metabolites by appending duplicate labels where appropriate
                # Metabolite.name = ifelse(dupcount > 1, paste0(Metabolite.name, "__", dupletter), Metabolite.name),
                metname.with.setpair = ifelse(dupcount2 > 1, paste0(Metabolite.name, "__", dupletter2, dupletter2),
                                              metname.with.setpair)) %>%
  dplyr::ungroup() %>%
  # dplyr::relocate(dupcount2, dupnum2, dupletter2, metname.with.setpair, .after = `Average.Rt.min.`) %>%
  dplyr::select(-dupcount2, -dupnum2, -dupletter2) %>%
  dplyr::relocate(dupletter, setpair, Average.Rt.min., metname.with.setpair, .after = metname.with.adduct)

neg_dup_df2 <- lneg_tempdat %>%
  dplyr::filter(setpair != "") %>%
  dplyr::select(Average.Rt.min., Average.Mz, Metabolite.name, Adduct.type, metname.with.adduct) %>%
  dplyr::arrange(Metabolite.name, Average.Rt.min.)

# Pos ---------------------------------------------------------------------

# Isolate cases where there are duplicated metabolite names with different adducts.
lpos_tempdat <- lpos_dat %>%
  dplyr::group_by(Metabolite.name) %>%
  dplyr::mutate(metname_dupcount = dplyr::n()) %>%
  dplyr::relocate(metname_dupcount, dupcount) %>%
  dplyr::filter(metname_dupcount > dupcount) %>%
  dplyr::ungroup()

# temp <- lpos_tempdat %>%
#   dplyr::mutate(round_Average.Rt.min. = round(Average.Rt.min., 2),
#                 round_Average.Rt.min. = dplyr::case_when(
#                   metname.with.adduct == "pos_Cer 34:2;O2|Cer 18:2;O2/16:0_[M+H-H2O]+" ~ 21.89,
#                   metname.with.adduct == "pos_Cer 36:2;O2|Cer 18:2;O2/18:0_[M+H-H2O]+" ~ 23.86,
#                   metname.with.adduct == "pos_Cer 43:1;O2|Cer 18:1;O2/25:0_[M+H-H2O]+" ~ 28.43,
#                   metname.with.adduct == "pos_Cer 43:2;O2|Cer 19:1;O2/24:1_[M+H]+" ~ 27.45,
#                   metname.with.adduct == "pos_HexCer 36:1;O2|HexCer 18:1;O2/18:0_[M+H]+" ~ 23.80,
#                   TRUE ~ round_Average.Rt.min.
#                 ),
#                 metname.with.rt = paste0(Metabolite.name, "_", round_Average.Rt.min.)) %>%
#   dplyr::arrange(Metabolite.name, round_Average.Rt.min.) %>%
#   dplyr::relocate(Average.Rt.min., round_Average.Rt.min., .after = metname.with.adduct)

# Note: Only true duplicates if the retention times are similar enough.
# After the summation option below, treat the unpaired cases with unique retention times
# as non-duplicates. 
# The setpair IDs below are used to denote which duplicates should be paired.
lpos_tempdat <- lpos_tempdat %>%
  dplyr::mutate(round_Average.Rt.min. = round(Average.Rt.min., 2),
                round_Average.Rt.min. = dplyr::case_when(
                  metname.with.adduct == "pos_Cer 34:2;O2|Cer 18:2;O2/16:0_[M+H-H2O]+" ~ 21.89,
                  metname.with.adduct == "pos_Cer 36:2;O2|Cer 18:2;O2/18:0_[M+H-H2O]+" ~ 23.86,
                  metname.with.adduct == "pos_Cer 43:1;O2|Cer 18:1;O2/25:0_[M+H-H2O]+" ~ 28.43,
                  metname.with.adduct == "pos_Cer 43:2;O2|Cer 19:1;O2/24:1_[M+H]+" ~ 27.45,
                  metname.with.adduct == "pos_HexCer 36:1;O2|HexCer 18:1;O2/18:0_[M+H]+" ~ 23.80,
                  TRUE ~ round_Average.Rt.min.
                ),
                metname.with.rt = paste0(Metabolite.name, "_", round_Average.Rt.min.)) %>%
  dplyr::group_by(metname.with.rt) %>%
  dplyr::arrange(round_Average.Rt.min.) %>% 
  dplyr::mutate(rtdupcount = dplyr::n(), # count total duplicates for a given metname.with.rt
                setpair = ifelse(rtdupcount > 1, "AA", "")) %>%
  dplyr::ungroup() %>%
  dplyr::select(-rtdupcount, -round_Average.Rt.min., -metname.with.rt) %>%
  dplyr::mutate(metname.with.setpair = ifelse(setpair != "", paste0(Metabolite.name, "__", setpair), metname.with.adduct)) %>%
  dplyr::group_by(metname.with.setpair) %>%
  dplyr::arrange(Average.Rt.min.) %>%
  dplyr::mutate(dupcount2 = dplyr::n(), # count total duplicates for a given metname.with.setpair
                dupnum2 = ceiling(dplyr::row_number()/2), # assign numeric labels for each duplicate
                dupletter2 = toupper(letters[dupnum2]), # convert numeric labels to capital letters
                # Rename metabolites by appending duplicate labels where appropriate
                # Metabolite.name = ifelse(dupcount > 1, paste0(Metabolite.name, "__", dupletter), Metabolite.name),
                metname.with.setpair = ifelse(dupcount2 > 1, paste0(Metabolite.name, "__", dupletter2, dupletter2),
                                              metname.with.setpair)) %>%
  dplyr::ungroup() %>%
  # dplyr::relocate(dupcount2, dupnum2, dupletter2, metname.with.setpair, .after = `Average.Rt.min.`) %>%
  dplyr::select(-dupcount2, -dupnum2, -dupletter2) %>%
  dplyr::relocate(dupletter, setpair, Average.Rt.min., metname.with.setpair, .after = metname.with.adduct)

pos_dup_df2 <- lpos_tempdat %>%
  dplyr::filter(setpair != "") %>%
  dplyr::select(Average.Rt.min., Average.Mz, Metabolite.name, Adduct.type, metname.with.adduct) %>%
  dplyr::arrange(Metabolite.name, Average.Rt.min.)

# -------------------------------------------------------------------------
## Summing Across Adducts (Sum) --------------------------------------------------
# Neg ---------------------------------------------------------------------

lneg_tempdat_pairs <- lneg_tempdat %>%
  dplyr::filter(metname_dupcount == 2) %>%
  dplyr::select(metname.with.setpair, 
                dplyr::contains("_Neg_")) %>%
  dplyr::mutate(Metabolite.name = gsub("__AA", "", metname.with.setpair)) %>%
  dplyr::select(-metname.with.setpair) %>%
  dplyr::relocate(Metabolite.name) %>%
  purrrlyr::slice_rows("Metabolite.name") %>%
  purrrlyr::dmap(sum)

lneg_tempdat_pairs_meta <- lneg_tempdat %>%
  dplyr::filter(metname_dupcount == 2) %>%
  dplyr::select(-dplyr::contains("_Neg_")) %>%
  dplyr::mutate(Metabolite.name = gsub("__AA", "", metname.with.setpair)) %>%
  dplyr::select(-dplyr::contains("metname"), -dupnum, -dupletter, -dupcount,-setpair) %>%
  dplyr::relocate(Metabolite.name) %>%
  purrrlyr::slice_rows("Metabolite.name") %>%
  purrrlyr::dmap(~ifelse(length(unique(.x)) == 1, as.character(unique(.x)), paste(.x, collapse = ", ")))

lneg_tempdat_pairs <- lneg_tempdat_pairs %>%
  dplyr::left_join(lneg_tempdat_pairs_meta, by = "Metabolite.name") %>%
  dplyr::relocate(names(lneg_tempdat_pairs_meta)) %>%
  dplyr::relocate(Metabolite.name)
rm(lneg_tempdat_pairs_meta)

lneg_tempdat_mults <- lneg_tempdat %>%
  dplyr::filter(metname_dupcount > 2) 

lneg_tempdat_mults_sum1 <- lneg_tempdat_mults %>%
  dplyr::filter(setpair == "") %>%
  dplyr::select(-Metabolite.name) %>%
  dplyr::rename(Metabolite.name = metname.with.setpair) %>%
  dplyr::select(-dplyr::contains("metname"), -dupnum, -dupletter, -dupcount,-setpair)  

lneg_tempdat_mults_sum2 <- lneg_tempdat_mults %>%
  dplyr::filter(setpair != "") %>%
  dplyr::select(metname.with.setpair,
                dplyr::contains("_Neg_")) %>%
  purrrlyr::slice_rows("metname.with.setpair") %>%
  purrrlyr::dmap(sum) %>%
  dplyr::rename(Metabolite.name = metname.with.setpair)

lneg_tempdat_mults_sum2_meta <- lneg_tempdat_mults %>%
  dplyr::filter(setpair != "") %>%
  dplyr::select(-dplyr::contains("_Neg_")) %>%
  dplyr::select(-Metabolite.name) %>%
  dplyr::rename(Metabolite.name = metname.with.setpair) %>%
  dplyr::select(-dplyr::contains("metname"), -dupnum, -dupletter, -dupcount,-setpair) %>%
  dplyr::relocate(Metabolite.name) %>%
  purrrlyr::slice_rows("Metabolite.name") %>%
  purrrlyr::dmap(~ifelse(length(unique(.x)) == 1, as.character(unique(.x)), paste(.x, collapse = ", ")))

lneg_tempdat_mults_sum2 <- lneg_tempdat_mults_sum2 %>%
  dplyr::left_join(lneg_tempdat_mults_sum2_meta, by = "Metabolite.name") %>%
  dplyr::relocate(names(lneg_tempdat_mults_sum2_meta)) %>%
  dplyr::relocate(Metabolite.name)
rm(lneg_tempdat_mults_sum2_meta)

lneg_tempdat_mults_sum3 <- rbind.data.frame(lneg_tempdat_mults_sum1,
                                            lneg_tempdat_mults_sum2) %>%
  dplyr::mutate(Metabolite.name2 = gsub("__.*", "", Metabolite.name),
                dupid = gsub(".*__", "", Metabolite.name)) %>%
  dplyr::relocate(Metabolite.name2, dupid, .after = Metabolite.name)

# Check that no lipid was mistakenly dropped
all(unique(lneg_tempdat$Metabolite.name) %in% c(gsub("_\\[.*", "", unique(lneg_tempdat_mults_sum3$Metabolite.name2)),
                                                gsub("_\\[.*", "", unique(lneg_tempdat_pairs$Metabolite.name))))

all(c(gsub("_\\[.*", "", unique(lneg_tempdat_mults_sum3$Metabolite.name2)),
      gsub("_\\[.*", "", unique(lneg_tempdat_pairs$Metabolite.name))) %in% unique(lneg_tempdat$Metabolite.name))

lneg_tempdat_mults_sum3 <- lneg_tempdat_mults_sum3 %>%
  dplyr::select(-Metabolite.name2, dupid)

lneg_dat_sum <- lneg_dat %>%
  dplyr::filter(!(metname.with.adduct %in% 
                    c(lneg_tempdat$metname.with.adduct))) %>%
  dplyr::select(-dupcount, -dupnum, -dupletter) %>%
  dplyr::mutate(Alignment.ID        = as.character(Alignment.ID),
                Average.Rt.min.     = as.character(Average.Rt.min.),
                Average.Mz          = as.character(Average.Mz),
                Reference.m.z       = as.character(Reference.m.z)) %>%#,
  # Fill..              = as.character(Fill..)) %>%
  dplyr::bind_rows(lneg_tempdat_pairs %>%
                     dplyr::mutate(Alignment.ID        = as.character(Alignment.ID),
                                   Average.Rt.min.     = as.character(Average.Rt.min.),
                                   Average.Mz          = as.character(Average.Mz),
                                   Reference.m.z       = as.character(Reference.m.z))) %>% #,
  # Fill..              = as.character(Fill..))) %>%
  dplyr::bind_rows(lneg_tempdat_mults_sum3 %>%
                     dplyr::mutate(Alignment.ID        = as.character(Alignment.ID),
                                   Average.Rt.min.     = as.character(Average.Rt.min.),
                                   Average.Mz          = as.character(Average.Mz),
                                   Reference.m.z       = as.character(Reference.m.z))) %>% #,
  # Fill..              = as.character(Fill..))) %>%
  dplyr::mutate(metname.with.adduct = ifelse(is.na(metname.with.adduct), Metabolite.name, metname.with.adduct))

length(unique(lneg_dat_sum$metname.with.adduct)) == nrow(lneg_dat_sum)

# Create accompanying edata and emeta objects
lneg_edat_sum <- lneg_dat_sum %>%
  dplyr::select(metname.with.adduct, 
                dplyr::contains("_Neg_"))
lneg_emet_sum <- lneg_dat_sum %>%
  dplyr::select(-dplyr::contains("_Neg_"))

rm(lneg_tempdat_mults_sum1, lneg_tempdat_mults_sum2, lneg_tempdat_mults)

# Pos ---------------------------------------------------------------------

lpos_tempdat_pairs <- lpos_tempdat %>%
  dplyr::filter(metname_dupcount == 2) %>%
  dplyr::select(metname.with.setpair,
                dplyr::contains("_Pos_")) %>%
  dplyr::mutate(Metabolite.name = gsub("__AA", "", metname.with.setpair)) %>%
  dplyr::select(-metname.with.setpair) %>%
  dplyr::relocate(Metabolite.name) %>%
  purrrlyr::slice_rows("Metabolite.name") %>%
  purrrlyr::dmap(sum)

lpos_tempdat_pairs_meta <- lpos_tempdat %>%
  dplyr::filter(metname_dupcount == 2) %>%
  dplyr::select(-dplyr::contains("_Pos_")) %>%
  dplyr::mutate(Metabolite.name = gsub("__AA", "", metname.with.setpair)) %>%
  dplyr::select(-dplyr::contains("metname"), -dupnum, -dupletter, -dupcount,-setpair) %>%
  dplyr::relocate(Metabolite.name) %>%
  purrrlyr::slice_rows("Metabolite.name") %>%
  purrrlyr::dmap(~ifelse(length(unique(.x)) == 1, as.character(unique(.x)), paste(.x, collapse = ", ")))

lpos_tempdat_pairs <- lpos_tempdat_pairs %>%
  dplyr::left_join(lpos_tempdat_pairs_meta, by = "Metabolite.name") %>%
  dplyr::relocate(names(lpos_tempdat_pairs_meta)) %>%
  dplyr::relocate(Metabolite.name)
rm(lpos_tempdat_pairs_meta)


lpos_tempdat_mults <- lpos_tempdat %>%
  dplyr::filter(metname_dupcount > 2) 

lpos_tempdat_mults_sum1 <- lpos_tempdat_mults %>%
  dplyr::filter(setpair == "") %>%
  dplyr::select(-Metabolite.name) %>%
  dplyr::rename(Metabolite.name = metname.with.setpair) %>%
  dplyr::select(-dplyr::contains("metname"), -dupnum, -dupletter, -dupcount,-setpair) 

lpos_tempdat_mults_sum2 <- lpos_tempdat_mults %>%
  dplyr::filter(setpair != "") %>%
  dplyr::select(metname.with.setpair,
                dplyr::contains("_Pos_")) %>%
  purrrlyr::slice_rows("metname.with.setpair") %>%
  purrrlyr::dmap(sum) %>%
  dplyr::rename(Metabolite.name = metname.with.setpair)

lpos_tempdat_mults_sum2_meta <- lpos_tempdat_mults %>%
  dplyr::filter(setpair != "") %>%
  dplyr::select(-dplyr::contains("_Pos_")) %>%
  dplyr::select(-Metabolite.name) %>%
  dplyr::rename(Metabolite.name = metname.with.setpair) %>%
  dplyr::select(-dplyr::contains("metname"), -dupnum, -dupletter, -dupcount,-setpair) %>%
  dplyr::relocate(Metabolite.name) %>%
  purrrlyr::slice_rows("Metabolite.name") %>%
  purrrlyr::dmap(~ifelse(length(unique(.x)) == 1, as.character(unique(.x)), paste(.x, collapse = ", ")))

lpos_tempdat_mults_sum2 <- lpos_tempdat_mults_sum2 %>%
  dplyr::left_join(lpos_tempdat_mults_sum2_meta, by = "Metabolite.name") %>%
  dplyr::relocate(names(lpos_tempdat_mults_sum2_meta)) %>%
  dplyr::relocate(Metabolite.name)
rm(lpos_tempdat_mults_sum2_meta)

lpos_tempdat_mults_sum3 <- rbind.data.frame(lpos_tempdat_mults_sum1,
                                            lpos_tempdat_mults_sum2) %>%
  dplyr::mutate(Metabolite.name2 = gsub("__.*", "", Metabolite.name),
                dupid = gsub(".*__", "", Metabolite.name)) %>%
  dplyr::relocate(Metabolite.name2, dupid, .after = Metabolite.name)

# Check that no lipid was mistakenly dropped
all(unique(lpos_tempdat$Metabolite.name) %in% c(gsub("_\\[.*", "", unique(lpos_tempdat_mults_sum3$Metabolite.name2)),
                                                gsub("_\\[.*", "", unique(lpos_tempdat_pairs$Metabolite.name))))

all(c(gsub("_\\[.*", "", unique(lpos_tempdat_mults_sum3$Metabolite.name2)),
      gsub("_\\[.*", "", unique(lpos_tempdat_pairs$Metabolite.name))) %in% unique(lpos_tempdat$Metabolite.name))


lpos_tempdat_mults_sum3 <- lpos_tempdat_mults_sum3 %>%
  dplyr::select(-Metabolite.name2, dupid)

lpos_dat_sum <- lpos_dat %>%
  dplyr::filter(!(metname.with.adduct %in% 
                    c(lpos_tempdat$metname.with.adduct))) %>%
  dplyr::select(-dupcount, -dupnum, -dupletter) %>%
  dplyr::mutate(Alignment.ID        = as.character(Alignment.ID),
                Average.Rt.min.     = as.character(Average.Rt.min.),
                Average.Mz          = as.character(Average.Mz),
                Reference.m.z       = as.character(Reference.m.z)) %>% #,
  # Fill..              = as.character(Fill..)) %>%
  dplyr::bind_rows(lpos_tempdat_pairs %>%
                     dplyr::mutate(Alignment.ID        = as.character(Alignment.ID),
                                   Average.Rt.min.     = as.character(Average.Rt.min.),
                                   Average.Mz          = as.character(Average.Mz),
                                   Reference.m.z       = as.character(Reference.m.z))) %>% #,
  # Fill..              = as.character(Fill..))) %>%
  dplyr::bind_rows(lpos_tempdat_mults_sum3 %>%
                     dplyr::mutate(Alignment.ID        = as.character(Alignment.ID),
                                   Average.Rt.min.     = as.character(Average.Rt.min.),
                                   Average.Mz          = as.character(Average.Mz),
                                   Reference.m.z       = as.character(Reference.m.z))) %>% #,
  # Fill..              = as.character(Fill..))) %>%
  dplyr::mutate(metname.with.adduct = ifelse(is.na(metname.with.adduct), Metabolite.name, metname.with.adduct))

length(unique(lpos_dat_sum$metname.with.adduct)) == nrow(lpos_dat_sum)

# Create accompanying edata and emeta objects
lpos_edat_sum <- lpos_dat_sum %>%
  dplyr::select(metname.with.adduct, 
                dplyr::contains("_Pos_"))
lpos_emet_sum <- lpos_dat_sum %>%
  dplyr::select(-dplyr::contains("_Pos_"))

rm(lpos_tempdat_mults_sum1, lpos_tempdat_mults_sum2, lpos_tempdat_mults)

# -------------------------------------------------------------------------
### Common Formatting of Modes ---------------------------------------------------------
## Differing Adduct Among Duplicates Ignored -------------------------------

# Make sure sample names are identical in pos and neg modes for when they are later combined
colnames(lneg_edat) %in% colnames(lpos_edat)
colnames(lpos_edat) %in% colnames(lneg_edat)

colnames(lneg_edat) <- gsub("_L_Neg_22Feb23_Crater.WCSH315305", "", colnames(lneg_edat))
colnames(lneg_edat) <- gsub("_Lumos_Neg_22Feb23_Crater.WCSH315305", "", colnames(lneg_edat))

colnames(lpos_edat) <- gsub("_L_Pos_18Feb23_Crater.WCSH315305", "", colnames(lpos_edat))
colnames(lpos_edat) <- gsub("_Lumos_Pos_18Feb23_Crater.WCSH315305", "", colnames(lpos_edat))

# All pos mode samples are in neg mode
colnames(lpos_edat)[!(colnames(lpos_edat) %in% colnames(lneg_edat))] 

# All neg mode samples are in pos mode
colnames(lneg_edat)[!(colnames(lneg_edat) %in% colnames(lpos_edat))] 

# Apply above changes in colnames to SampleIDs in the respective fdatasets
lneg_fdat <- lneg_fdat %>%
  dplyr::mutate(SampleID = gsub("_L_Neg_22Feb23_Crater.WCSH315305", "", SampleID),
                SampleID = gsub("_Lumos_Neg_22Feb23_Crater.WCSH315305", "", SampleID))

all(lneg_fdat$SampleID %in% colnames(lneg_edat)[-1])

lpos_fdat <- lpos_fdat %>%
  dplyr::mutate(SampleID = gsub("_L_Pos_18Feb23_Crater.WCSH315305", "", SampleID),
                SampleID = gsub("_Lumos_Pos_18Feb23_Crater.WCSH315305", "", SampleID))

all(lpos_fdat$SampleID %in% colnames(lpos_edat)[-1])

## Sum ---------------------------------------------------------------------

# Make sure sample names are identical in pos and neg modes for when they are later combined
colnames(lneg_edat_sum) %in% colnames(lpos_edat_sum)
colnames(lpos_edat_sum) %in% colnames(lneg_edat_sum)

colnames(lneg_edat_sum) <- gsub("_L_Neg_22Feb23_Crater.WCSH315305", "", colnames(lneg_edat_sum))
colnames(lneg_edat_sum) <- gsub("_Lumos_Neg_22Feb23_Crater.WCSH315305", "", colnames(lneg_edat_sum))

colnames(lpos_edat_sum) <- gsub("_L_Pos_18Feb23_Crater.WCSH315305", "", colnames(lpos_edat_sum))
colnames(lpos_edat_sum) <- gsub("_Lumos_Pos_18Feb23_Crater.WCSH315305", "", colnames(lpos_edat_sum))

# All pos mode samples are in neg mode
colnames(lpos_edat_sum)[!(colnames(lpos_edat_sum) %in% colnames(lneg_edat_sum))] 

# All neg mode samples are in pos mode
colnames(lneg_edat_sum)[!(colnames(lneg_edat_sum) %in% colnames(lpos_edat_sum))] 


all(lneg_fdat$SampleID %in% colnames(lneg_edat_sum)[-1])

all(lpos_fdat$SampleID %in% colnames(lpos_edat_sum)[-1])

# -------------------------------------------------------------------------

### Process Blanks Assessment (insufficient data to do for these data) -----------------------------------------------

## Negative Mode -----------------------------------------------------------

# # We take the average values across the blanks and compare to the 
# # medacross the analytic samples 
# pb_means <- apply(lneg_edat_sum %>% 
#                     dplyr::select(dplyr::contains("PB")), 1, 
#                   function(x){mean(as.numeric(x), na.rm = TRUE)})
# samp_means <- apply(lneg_edat_sum %>% 
#                       dplyr::select(-dplyr::contains("PB"),
#                                     -dplyr::contains("QC"),
#                                     -dplyr::contains("Blank")), 1, 
#                     function(x){mean(as.numeric(x), na.rm = TRUE)})
# 
# # Graphic to visualize all pointwise differences (PB1)
# tempdf_pb1 <- lneg_edat_sum %>% 
#   dplyr::select(-dplyr::contains("PB"),
#                 -dplyr::contains("QC"),
#                 -dplyr::contains("Blank"))
# for(i in 2:ncol(tempdf_pb1)){
#   # Want to remove the cases where both are 0
#   temp <- mapply(function(x,y){
#     if(x == 0 & y == 0){
#       return(NA)
#     } else{
#       x-y
#     }
#   }, tempdf_pb1[,i, drop = TRUE], lneg_edat_sum[,"Lip_PB_01",drop = TRUE])
#   tempdf_pb1[,i] <- temp
# }
# tempdf_pb1 <- tempdf_pb1 %>%
#   tidyr::pivot_longer(!metname.with.adduct, names_to = "Sample", values_to = "Difference") %>%
#   as.data.frame() %>%
#   dplyr::filter(!is.na(Difference)) %>%
#   dplyr::group_by(metname.with.adduct) %>%
#   dplyr::mutate(idx = dplyr::row_number()) %>%
#   dplyr::ungroup() %>%
#   dplyr::mutate(Diffcat = cut(Difference, breaks = c(min(Difference), 0, 5000, max(Difference)), include.lowest = TRUE))
# 
# tempdf_pb1_counts <- tempdf_pb1 %>%
#   dplyr::group_by(metname.with.adduct, Diffcat) %>%
#   dplyr::summarise(count = dplyr::n())
# 
# # png(filename = here("reports", "final", "report_figures", "kftall_lipids_worm", "pb1_negmode_pointdiff.png"), units = "in", width = 16, height = 5, res = 300)
# neg_p1 <- ggplot(data = tempdf_pb1_counts, aes(fill = Diffcat, y = count, x = metname.with.adduct)) +
#   geom_bar(position = "stack", stat = "identity") + theme(axis.text.x = element_blank()) +
#   ylab("Samples") + ggtitle("Intensity Differences: PB 01 (Neg Mode)") + xlab("Metabolite")
# # dev.off()
# 
# # Export cases where there is a negative diffcat
# tempdf_pb1_export <- tempdf_pb1 %>%
#   dplyr::filter(Difference < 0) %>%
#   dplyr::select(-Diffcat)
# 
# # Graphic to visualize all pointwise differences (PB2)
# tempdf_pb2 <- lneg_edat_sum %>% 
#   dplyr::select(-dplyr::contains("PB"),
#                 -dplyr::contains("QC"),
#                 -dplyr::contains("Blank"))
# for(i in 2:ncol(tempdf_pb2)){
#   # Want to remove the cases where both are 0
#   temp <- mapply(function(x,y){
#     if(x == 0 & y == 0){
#       return(NA)
#     } else{
#       x-y
#     }
#   }, tempdf_pb2[,i, drop = TRUE], lneg_edat_sum[,"Lip_PB_02",drop = TRUE])
#   tempdf_pb2[,i] <- temp
# }
# tempdf_pb2 <- tempdf_pb2 %>%
#   tidyr::pivot_longer(!metname.with.adduct, names_to = "Sample", values_to = "Difference") %>%
#   as.data.frame() %>%
#   dplyr::filter(!is.na(Difference)) %>%
#   dplyr::group_by(metname.with.adduct) %>%
#   dplyr::mutate(idx = dplyr::row_number()) %>%
#   dplyr::ungroup() %>%
#   dplyr::mutate(Diffcat = cut(Difference, breaks = c(min(Difference), 0, 5000, max(Difference)), include.lowest = TRUE))
# 
# tempdf_pb2_counts <- tempdf_pb2 %>%
#   dplyr::group_by(metname.with.adduct, Diffcat) %>%
#   dplyr::summarise(count = dplyr::n())
# 
# # png(filename = here("reports", "final", "report_figures", "kftall_lipids_worm", "pb2_negmode_pointdiff.png"), units = "in", width = 16, height = 5, res = 300)
# neg_p2 <- ggplot(data = tempdf_pb2_counts, aes(fill = Diffcat, y = count, x = metname.with.adduct)) +
#   geom_bar(position = "stack", stat = "identity") + theme(axis.text.x = element_blank()) +
#   ylab("Samples") + ggtitle("Intensity Differences: PB 02 (Neg Mode)") + xlab("Metabolite")
# # dev.off()
# 
# # Export cases where there is a negative diffcat
# tempdf_pb2_export <- tempdf_pb2 %>%
#   dplyr::filter(Difference < 0) %>%
#   dplyr::select(-Diffcat)
# 
# # Graphic to visualize all pointwise differences (PB3)
# tempdf_pb3 <- lneg_edat_sum %>% 
#   dplyr::select(-dplyr::contains("PB"),
#                 -dplyr::contains("QC"),
#                 -dplyr::contains("Blank"))
# for(i in 2:ncol(tempdf_pb3)){
#   # Want to remove the cases where both are 0
#   temp <- mapply(function(x,y){
#     if(x == 0 & y == 0){
#       return(NA)
#     } else{
#       x-y
#     }
#   }, tempdf_pb3[,i, drop = TRUE], lneg_edat_sum[,"Lip_PB_03",drop = TRUE])
#   tempdf_pb3[,i] <- temp
# }
# tempdf_pb3 <- tempdf_pb3 %>%
#   tidyr::pivot_longer(!metname.with.adduct, names_to = "Sample", values_to = "Difference") %>%
#   as.data.frame() %>%
#   dplyr::filter(!is.na(Difference)) %>%
#   dplyr::group_by(metname.with.adduct) %>%
#   dplyr::mutate(idx = dplyr::row_number()) %>%
#   dplyr::ungroup() %>%
#   dplyr::mutate(Diffcat = cut(Difference, breaks = c(min(Difference), 0, 5000, max(Difference)), include.lowest = TRUE))
# 
# tempdf_pb3_counts <- tempdf_pb3 %>%
#   dplyr::group_by(metname.with.adduct, Diffcat) %>%
#   dplyr::summarise(count = dplyr::n())
# 
# # png(filename = here("reports", "final", "report_figures", "kftall_lipids_worm", "pb3_negmode_pointdiff.png"), units = "in", width = 16, height = 5, res = 300)
# neg_p3 <- ggplot(data = tempdf_pb3_counts, aes(fill = Diffcat, y = count, x = metname.with.adduct)) +
#   geom_bar(position = "stack", stat = "identity") + theme(axis.text.x = element_blank()) +
#   ylab("Samples") + ggtitle("Intensity Differences: PB 03 (Neg Mode)") + xlab("Metabolite")
# # dev.off()
# 
# # Export cases where there is a negative diffcat
# tempdf_pb3_export <- tempdf_pb3 %>%
#   dplyr::filter(Difference < 0) %>%
#   dplyr::select(-Diffcat)
# 
# tempdf_pb1_export_neg <- tempdf_pb1_export
# tempdf_pb2_export_neg <- tempdf_pb2_export
# tempdf_pb3_export_neg <- tempdf_pb3_export

## Positive Mode -----------------------------------------------------------

# # We take the average values across the blanks and compare to the 
# # medacross the analytic samples 
# pb_means <- apply(lpos_edat_sum %>% 
#                     dplyr::select(dplyr::contains("PB")), 1, 
#                   function(x){mean(as.numeric(x), na.rm = TRUE)})
# samp_means <- apply(lpos_edat_sum %>% 
#                       dplyr::select(-dplyr::contains("PB"),
#                                     -dplyr::contains("QC"),
#                                     -dplyr::contains("Blank")), 1, 
#                     function(x){mean(as.numeric(x), na.rm = TRUE)})
# 
# # Graphic to visualize all pointwise differences (PB1)
# tempdf_pb1 <- lpos_edat_sum %>% 
#   dplyr::select(-dplyr::contains("PB"),
#                 -dplyr::contains("QC"),
#                 -dplyr::contains("Blank"))
# for(i in 2:ncol(tempdf_pb1)){
#   # Want to remove the cases where both are 0
#   temp <- mapply(function(x,y){
#     if(x == 0 & y == 0){
#       return(NA)
#     } else{
#       x-y
#     }
#   }, tempdf_pb1[,i, drop = TRUE], lpos_edat_sum[,"Lip_PB_01",drop = TRUE])
#   tempdf_pb1[,i] <- temp
# }
# tempdf_pb1 <- tempdf_pb1 %>%
#   tidyr::pivot_longer(!metname.with.adduct, names_to = "Sample", values_to = "Difference") %>%
#   as.data.frame() %>%
#   dplyr::filter(!is.na(Difference)) %>%
#   dplyr::group_by(metname.with.adduct) %>%
#   dplyr::mutate(idx = dplyr::row_number()) %>%
#   dplyr::ungroup() %>%
#   dplyr::mutate(Diffcat = cut(Difference, breaks = c(min(Difference), 0, 5000, max(Difference)), include.lowest = TRUE))
# 
# tempdf_pb1_counts <- tempdf_pb1 %>%
#   dplyr::group_by(metname.with.adduct, Diffcat) %>%
#   dplyr::summarise(count = dplyr::n())
# 
# # png(filename = here("reports", "final", "report_figures", "kftall_lipids_worm", "pb1_posmode_pointdiff.png"), units = "in", width = 16, height = 5, res = 300)
# pos_p1 <- ggplot(data = tempdf_pb1_counts, aes(fill = Diffcat, y = count, x = metname.with.adduct)) +
#   geom_bar(position = "stack", stat = "identity") + theme(axis.text.x = element_blank()) +
#   ylab("Samples") + ggtitle("Intensity Differences: PB 01 (Pos Mode)") + xlab("Metabolite")
# # dev.off()
# 
# # Export cases where there is a positive diffcat
# tempdf_pb1_export <- tempdf_pb1 %>%
#   dplyr::filter(Difference < 0) %>%
#   dplyr::select(-Diffcat)
# 
# # Graphic to visualize all pointwise differences (PB2)
# tempdf_pb2 <- lpos_edat_sum %>% 
#   dplyr::select(-dplyr::contains("PB"),
#                 -dplyr::contains("QC"),
#                 -dplyr::contains("Blank"))
# for(i in 2:ncol(tempdf_pb2)){
#   # Want to remove the cases where both are 0
#   temp <- mapply(function(x,y){
#     if(x == 0 & y == 0){
#       return(NA)
#     } else{
#       x-y
#     }
#   }, tempdf_pb2[,i, drop = TRUE], lpos_edat_sum[,"Lip_PB_02",drop = TRUE])
#   tempdf_pb2[,i] <- temp
# }
# tempdf_pb2 <- tempdf_pb2 %>%
#   tidyr::pivot_longer(!metname.with.adduct, names_to = "Sample", values_to = "Difference") %>%
#   as.data.frame() %>%
#   dplyr::filter(!is.na(Difference)) %>%
#   dplyr::group_by(metname.with.adduct) %>%
#   dplyr::mutate(idx = dplyr::row_number()) %>%
#   dplyr::ungroup() %>%
#   dplyr::mutate(Diffcat = cut(Difference, breaks = c(min(Difference), 0, 5000, max(Difference)), include.lowest = TRUE))
# 
# tempdf_pb2_counts <- tempdf_pb2 %>%
#   dplyr::group_by(metname.with.adduct, Diffcat) %>%
#   dplyr::summarise(count = dplyr::n())
# 
# # png(filename = here("reports", "final", "report_figures", "kftall_lipids_worm", "pb2_posmode_pointdiff.png"), units = "in", width = 16, height = 5, res = 300)
# pos_p2 <- ggplot(data = tempdf_pb2_counts, aes(fill = Diffcat, y = count, x = metname.with.adduct)) +
#   geom_bar(position = "stack", stat = "identity") + theme(axis.text.x = element_blank()) +
#   ylab("Samples") + ggtitle("Intensity Differences: PB 02 (Pos Mode)") + xlab("Metabolite")
# # dev.off()
# 
# # Export cases where there is a positive diffcat
# tempdf_pb2_export <- tempdf_pb2 %>%
#   dplyr::filter(Difference < 0) %>%
#   dplyr::select(-Diffcat)
# 
# # Graphic to visualize all pointwise differences (PB3)
# tempdf_pb3 <- lpos_edat_sum %>% 
#   dplyr::select(-dplyr::contains("PB"),
#                 -dplyr::contains("QC"),
#                 -dplyr::contains("Blank"))
# for(i in 2:ncol(tempdf_pb3)){
#   # Want to remove the cases where both are 0
#   temp <- mapply(function(x,y){
#     if(x == 0 & y == 0){
#       return(NA)
#     } else{
#       x-y
#     }
#   }, tempdf_pb3[,i, drop = TRUE], lpos_edat_sum[,"Lip_PB_03",drop = TRUE])
#   tempdf_pb3[,i] <- temp
# }
# tempdf_pb3 <- tempdf_pb3 %>%
#   tidyr::pivot_longer(!metname.with.adduct, names_to = "Sample", values_to = "Difference") %>%
#   as.data.frame() %>%
#   dplyr::filter(!is.na(Difference)) %>%
#   dplyr::group_by(metname.with.adduct) %>%
#   dplyr::mutate(idx = dplyr::row_number()) %>%
#   dplyr::ungroup() %>%
#   dplyr::mutate(Diffcat = cut(Difference, breaks = c(min(Difference), 0, 5000, max(Difference)), include.lowest = TRUE))
# 
# tempdf_pb3_counts <- tempdf_pb3 %>%
#   dplyr::group_by(metname.with.adduct, Diffcat) %>%
#   dplyr::summarise(count = dplyr::n())
# 
# # png(filename = here("reports", "final", "report_figures", "kftall_lipids_worm", "pb3_posmode_pointdiff.png"), units = "in", width = 16, height = 5, res = 300)
# pos_p3 <- ggplot(data = tempdf_pb3_counts, aes(fill = Diffcat, y = count, x = metname.with.adduct)) +
#   geom_bar(position = "stack", stat = "identity") + theme(axis.text.x = element_blank()) +
#   ylab("Samples") + ggtitle("Intensity Differences: PB 03 (Pos Mode)") + xlab("Metabolite")
# # dev.off()
# 
# # Export cases where there is a positive diffcat
# tempdf_pb3_export <- tempdf_pb3 %>%
#   dplyr::filter(Difference < 0) %>%
#   dplyr::select(-Diffcat)
# 
# tempdf_pb1_export_pos <- tempdf_pb1_export
# tempdf_pb2_export_pos <- tempdf_pb2_export
# tempdf_pb3_export_pos <- tempdf_pb3_export


# -------------------------------------------------------------------------

# png(filename = here("figures", "ptrc_lipids", "pb_negmode_pointdiff.png"), units = "in", width = 16, height = 10, res = 300)
# ggpubr::ggarrange(neg_p1, neg_p2, neg_p3, nrow = 3, ncol = 1, common.legend = TRUE,
#                   legend = "right")
# dev.off()
# 
# png(filename = here("figures", "ptrc_lipids", "pb_posmode_pointdiff.png"), units = "in", width = 16, height = 10, res = 300)
# ggpubr::ggarrange(pos_p1, pos_p2, pos_p3, nrow = 3, ncol = 1, common.legend = TRUE,
#                   legend = "right")
# dev.off()
# 
# pb_pointdiffs_imgs <- c(here("figures", "ptrc_lipids", "pb_negmode_pointdiff.png"),
#                         here("figures", "ptrc_lipids", "pb_posmode_pointdiff.png"))
# 
# rm(pb_pointdiffs_imgs)

# -------------------------------------------------------------------------


### Flagging (not done for these data) ----------------------------------------------------------------
## Sum ---------------------------------------------------------------------

# # Sum: Check if any of the lipids are flagged for the QC sample and sample intensity check
# qcsamp_check <- apply(lneg_dat_sum %>%
#                         dplyr::select(dplyr::contains("QC_Pool")),
#                       1, function(x){sum(x<3000)/length(x)*100})
# lneg_emet_sum$QCpool_flag <- as.numeric(qcsamp_check >= 40) 
# 
# samp_check <- apply(lneg_dat_sum %>%
#                       dplyr::select(-dplyr::contains("QC_Pool"),
#                                     -dplyr::contains("PB"),
#                                     -dplyr::contains("Blank")) %>%
#                       dplyr::select(dplyr::contains("PTRC")),
#                     1, function(x){sum(x<3000)/length(x)*100})
# lneg_emet_sum$Sample_flag <- as.numeric(samp_check >= 40) 
# rm(qcsamp_check, samp_check)
# 
# qcsamp_check <- apply(lpos_dat_sum %>%
#                         dplyr::select(dplyr::contains("QC_Pool")),
#                       1, function(x){sum(x<5000)/length(x)*100})
# lpos_emet_sum$QCpool_flag <- as.numeric(qcsamp_check >= 40) 
# 
# samp_check <- apply(lpos_dat_sum %>%
#                       dplyr::select(-dplyr::contains("QC_Pool"),
#                                     -dplyr::contains("PB"),
#                                     -dplyr::contains("Blank")) %>%
#                       dplyr::select(dplyr::contains("PTRC")),
#                     1, function(x){sum(x<5000)/length(x)*100})
# lpos_emet_sum$Sample_flag <- as.numeric(samp_check >= 40) 
# rm(qcsamp_check, samp_check)


## Choice ------------------------------------------------------------------

# # Choice/Higher: Check if any of the lipids are flagged for the QC sample and sample intensity check
# qcsamp_check <- apply(lneg_dat_choice %>%
#                         dplyr::select(dplyr::contains("QC_Pool")),
#                       1, function(x){sum(x<3000)/length(x)*100})
# lneg_emet_choice$QCpool_flag <- as.numeric(qcsamp_check >= 40) 
# 
# samp_check <- apply(lneg_dat_choice %>%
#                       dplyr::select(-dplyr::contains("QC_Pool"),
#                                     -dplyr::contains("PB"),
#                                     -dplyr::contains("Blank")) %>%
#                       dplyr::select(dplyr::contains("PTRC")),
#                     1, function(x){sum(x<3000)/length(x)*100})
# lneg_emet_choice$Sample_flag <- as.numeric(samp_check >= 40) 
# rm(qcsamp_check, samp_check)
# 
# qcsamp_check <- apply(lpos_dat_choice %>%
#                         dplyr::select(dplyr::contains("QC_Pool")),
#                       1, function(x){sum(x<5000)/length(x)*100})
# lpos_emet_choice$QCpool_flag <- as.numeric(qcsamp_check >= 40) 
# 
# samp_check <- apply(lpos_dat_choice %>%
#                       dplyr::select(-dplyr::contains("QC_Pool"),
#                                     -dplyr::contains("PB"),
#                                     -dplyr::contains("Blank")) %>%
#                       dplyr::select(dplyr::contains("PTRC")),
#                     1, function(x){sum(x<5000)/length(x)*100})
# lpos_emet_choice$Sample_flag <- as.numeric(samp_check >= 40) 
# rm(qcsamp_check, samp_check)

# -------------------------------------------------------------------------

# neg_wid <- slick_div(
#   x = DT::datatable(data = lneg_emet_sum %>%
#                       dplyr::filter(QCpool_flag == 1 | Sample_flag == 1) %>%
#                       dplyr::select(Average.Rt.min., Average.Mz, Metabolite.name, Adduct.type, metname.with.adduct, QCpool_flag, Sample_flag) %>%
#                       dplyr::rename(Lipid.name = Metabolite.name,
#                                     lipname.with.adduct = metname.with.adduct), rownames = FALSE,
#                     options = list(scrollX = TRUE)),
#   css = htmltools::css(height = '490px', width = '800px', marginLeft='auto',marginRight='auto')
# )
# 
# pos_wid <- slick_div(
#   x = DT::datatable(data = lpos_emet_sum %>%
#                       dplyr::filter(QCpool_flag == 1 | Sample_flag == 1) %>%
#                       dplyr::select(Average.Rt.min., Average.Mz, Metabolite.name, Adduct.type, metname.with.adduct, QCpool_flag, Sample_flag) %>%
#                       dplyr::rename(Lipid.name = Metabolite.name,
#                                     lipname.with.adduct = metname.with.adduct), rownames = FALSE,
#                     options = list(scrollX = TRUE)),
#   css = htmltools::css(height = '490px', width = '800px', marginLeft='auto',marginRight='auto')
# )
# 
# doms <- slick_list(neg_wid, pos_wid)
# 
# rm(neg_wid, pos_wid, doms)

### CV Computation ----------------------------------------------------------------
## Sum ---------------------------------------------------------------------

# Sum: Check if any of the lipids are flagged for the QC sample and sample intensity check
qccv_check <- apply(lneg_dat_sum %>%
                      dplyr::select(dplyr::contains("_QC_")),
                    1, function(x){sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE)*100})
lneg_emet_sum$QC_CV <- as.numeric(qccv_check) 
rm(qccv_check)

qccv_check <- apply(lpos_dat_sum %>%
                      dplyr::select(dplyr::contains("_QC_")),
                    1, function(x){sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE)*100})
lpos_emet_sum$QC_CV <- as.numeric(qccv_check) 
rm(qccv_check)

# -------------------------------------------------------------------------


### Log2-transform / 0-to-NA ------------------------------------------------
## Sum ---------------------------------------------------------------------

# 318 instances of 0 have been replaced with NA
pm_lneg_sum <- as.lipidData(e_data = lneg_edat_sum, 
                            f_data = lneg_fdat, 
                            e_meta = lneg_emet_sum,
                            edata_cname = "metname.with.adduct", 
                            fdata_cname = "SampleID",
                            emeta_cname = "metname.with.adduct",
                            data_scale = "abundance",
                            data_types = "Negative")


# 1022 instances of 0 have been replaced with NA
pm_lpos_sum <- as.lipidData(e_data = lpos_edat_sum, 
                            f_data = lpos_fdat, 
                            e_meta = lpos_emet_sum,
                            edata_cname = "metname.with.adduct", 
                            fdata_cname = "SampleID",
                            emeta_cname = "metname.with.adduct",
                            data_scale = "abundance",
                            data_types = "Positive")


# Transform data to the log2 scale
pm_lneg_sum <- edata_transform(pm_lneg_sum, data_scale = "log2")
pm_lpos_sum <- edata_transform(pm_lpos_sum, data_scale = "log2")

# -------------------------------------------------------------------------
### Export Data (No Run Order Correction) -------------------------------------------------------------
# Sum ---------------------------------------------------------------------

lneg_sum_export_dat <- pm_lneg_sum$e_data %>%
  dplyr::left_join(pm_lneg_sum$e_meta) %>%
  dplyr::relocate(Alignment.ID:QC_CV, .after = metname.with.adduct) %>%
  dplyr::select(-dupid) %>%
  dplyr::mutate(Metabolite.name = gsub("neg_", "", Metabolite.name),
                metname.with.adduct = gsub("neg_", "", metname.with.adduct)) %>% #,
  # metname.with.adduct = gsub("_\\[M-H\\]-", "", metname.with.adduct),
  # metname.with.adduct = gsub("_\\[M\\+CH3COO\\]-", "", metname.with.adduct)) %>%
  dplyr::rename(original_name = Metabolite.name,
                formatted_name = metname.with.adduct)


lpos_sum_export_dat <- pm_lpos_sum$e_data %>%
  dplyr::left_join(pm_lpos_sum$e_meta) %>%
  dplyr::relocate(Alignment.ID:QC_CV, .after = metname.with.adduct) %>%
  dplyr::select(-dupid) %>%
  dplyr::mutate(Metabolite.name = gsub("pos_", "", Metabolite.name),
                metname.with.adduct = gsub("pos_", "", metname.with.adduct)) %>% #,
  # metname.with.adduct = gsub("_\\[M\\+H-H2O\\]\\+", "", metname.with.adduct),
  # metname.with.adduct = gsub("_\\[M\\+H\\]\\+", "", metname.with.adduct),
  # metname.with.adduct = gsub("_\\[M\\+Na\\]\\+", "", metname.with.adduct),
  # metname.with.adduct = gsub("_\\[M\\+NH4\\]\\+", "", metname.with.adduct)) %>%
  dplyr::rename(original_name = Metabolite.name,
                formatted_name = metname.with.adduct)


writexl::write_xlsx(x = list(`Negative` = lneg_sum_export_dat ,
                             `Positive` = lpos_sum_export_dat),
                    path = here("beataml_lipids_log2_sum.xlsx"))
saveRDS(lneg_sum_export_dat, here("beataml_lipids_neg_log2_sum.rds"))
saveRDS(lpos_sum_export_dat, here("beataml_lipids_pos_log2_sum.rds"))


# -------------------------------------------------------------------------

# Save fdat
saveRDS(lneg_fdat, here("beataml_lipids_neg_fdat.rds"))
saveRDS(lpos_fdat, here("beataml_lipids_pos_fdat.rds"))

