library(magrittr)
library(pmartR)
library(ggplot2)
library(here)

i_am("scripts/ptrc_lipids_proc.R")

# refer to ptrc_lipids_report for corresponding report.

##### Data Processing ---------------------------------------------------------
#### Import ------------------------------------------------------------------

data_filenames <- list.files(here("data"))

# Neg ---------------------------------------------------------------------

## Skip first row since contains fdata
lneg_dat <- readxl::read_excel(path = here("data", "PTRC_lipids_NEG_for_stats.xlsx"),
                               sheet = "Height_0_2025630936_NEG",
                               skip = 1)
colnames(lneg_dat) <- make.names(colnames(lneg_dat))
## Extract first three rows only to obtain precursor to fdata
lneg_dat_samplabels <- readxl::read_excel(path = here("data", "PTRC_lipids_NEG_for_stats.xlsx"),
                                          sheet = "Height_0_2025630936_NEG",
                                          n_max = 2,
                                          col_names = FALSE)

# Pos ---------------------------------------------------------------------

## Skip first row since contains fdata
lpos_dat <- readxl::read_excel(path = here("data", "PTRC_lipids_POS_for_stats.xlsx"),
                               sheet = "Height_0_20256261733_POS",
                               skip = 1)
colnames(lpos_dat) <- make.names(colnames(lpos_dat))
## Extract first two rows only to obtain precursor to fdata
lpos_dat_samplabels <- readxl::read_excel(path = here("data", "PTRC_lipids_POS_for_stats.xlsx"),
                                          sheet = "Height_0_20256261733_POS",
                                          n_max = 2,
                                          col_names = FALSE)

# 11/13/2025: Add in alignment columns for later integration with BeatAML data.
# Same is not done for neg mode because neg mode alignment was done after the 
# processing in this script, whereas the same was not done for the positive mode
# data for whatever reason.
lpos_dat_align <- readxl::read_excel(path = here("data", "PTRC_lipids_POS_for_stats_aligned_with_BEAT.xlsx"),
                                     skip = 1) 
colnames(lpos_dat_align) <- make.names(colnames(lpos_dat_align))
lpos_dat_align <- lpos_dat_align %>%
  dplyr::select(Alignment.ID, BEAT.AML.Alignment.ID, BEAT.AML.RT)
lpos_dat <- lpos_dat %>%
  dplyr::left_join(lpos_dat_align) %>%
  dplyr::relocate(BEAT.AML.Alignment.ID, BEAT.AML.RT, .after = Metabolite.name)

# -------------------------------------------------------------------------

#### Filtering ---------------------------------------------------------------
### S/N Filter ---------------------------------------------------------------------
# Neg ---------------------------------------------------------------------

filt_count_neg <- nrow(lneg_dat %>%
                         dplyr::filter(S.N.average < 10))

lneg_dat <- lneg_dat %>%
  dplyr::filter(S.N.average >= 10)

# Pos ---------------------------------------------------------------------

filt_count_pos <- nrow(lpos_dat %>%
                         dplyr::filter(S.N.average < 10))

lpos_dat <- lpos_dat %>%
  dplyr::filter(S.N.average >= 10)

# -------------------------------------------------------------------------

### Process Blank / Instrument Blank Filter ---------------------------------
# Neg ---------------------------------------------------------------------

# Filter out if average process blank value for lipid is greater than 10% or more of the
# values in the samples, remove it.
# Compare between using the average of the process blanks to the average of the instrument
# blanks
pb_means <- apply(lneg_dat %>% 
                    dplyr::select(dplyr::contains("PB")), 1, function(x){mean(as.numeric(x), na.rm = TRUE)})
ib_means <- apply(lneg_dat %>% 
                    dplyr::select(dplyr::contains("Blank_Co")), 1, function(x){mean(as.numeric(x), na.rm = TRUE)})
blank_comparison <- vector("list", length = length(pb_means))
for(i in 1:length(pb_means)){
  tempdat <- as.numeric(lneg_dat %>%
                          dplyr::select(dplyr::contains("PTRC_exp26_Lip"),
                                        -dplyr::contains("PB"),
                                        -dplyr::contains("Blank_Co")) %>%
                          .[i,])
  
  pb_prop <- sum(pb_means[i] > tempdat)/length(tempdat)*100
  ib_prop <- sum(ib_means[i] > tempdat)/length(tempdat)*100
  blank_comparison[[i]] <- data.frame(Metabolite.name = lneg_dat$Metabolite.name[i],
                                      Average.Rt.min. = lneg_dat$Average.Rt.min.[i],
                                      Adduct.type = lneg_dat$Adduct.type[i],
                                      pb_prop = pb_prop,
                                      ib_prop = ib_prop,
                                      pb_flag = ifelse(pb_prop >= 10, 1, 0),
                                      ib_flag = ifelse(ib_prop >= 10, 1, 0))
}
blank_comparison <- Reduce("rbind", blank_comparison) %>%
  dplyr::mutate(join_flag = ifelse(pb_flag == 1 & ib_flag == 1,1,0),
                different_flag = ifelse(pb_flag == ib_flag, 0, 1))

blank_comparison_neg <- blank_comparison %>%
  dplyr::mutate(blank_remove = ifelse(pb_flag == 1 | ib_flag == 1, 1, 0)) %>%
  dplyr::select(-pb_prop, -ib_prop)

blank_comparison <- blank_comparison %>%
  dplyr::mutate(blank_remove = ifelse(pb_flag == 1 | ib_flag == 1, 1, 0)) %>%
  dplyr::select(Metabolite.name, Average.Rt.min., Adduct.type, blank_remove)

# As per JK 4/25/24...remove based on both pb and ib
lneg_dat <- lneg_dat %>%
  dplyr::left_join(blank_comparison) %>%
  dplyr::filter(blank_remove == 0) %>%
  dplyr::select(-blank_remove)
nrow(blank_comparison_neg) - nrow(lneg_dat)

rm(i, ib_means, ib_prop, pb_means, pb_prop, tempdat, blank_comparison)

# Pos ---------------------------------------------------------------------

# Filter out if average process blank value for lipid is greater than 10% or more of the
# values in the samples, remove it.
# Compare between using the average of the process blanks to the average of the instrument
# blanks
pb_means <- apply(lpos_dat %>% 
                    dplyr::select(dplyr::contains("PB")), 1, function(x){mean(as.numeric(x), na.rm = TRUE)})
ib_means <- apply(lpos_dat %>% 
                    dplyr::select(dplyr::contains("Blank_Co")), 1, function(x){mean(as.numeric(x), na.rm = TRUE)})
blank_comparison <- vector("list", length = length(pb_means))
for(i in 1:length(pb_means)){
  tempdat <- as.numeric(lpos_dat %>%
                          dplyr::select(dplyr::contains("PTRC_exp26_Lip"),
                                        -dplyr::contains("PB"),
                                        -dplyr::contains("Blank_Co")) %>%
                          .[i,])
  pb_prop <- sum(pb_means[i] > tempdat)/length(tempdat)*100
  ib_prop <- sum(ib_means[i] > tempdat)/length(tempdat)*100
  blank_comparison[[i]] <- data.frame(Metabolite.name = lpos_dat$Metabolite.name[i],
                                      Average.Rt.min. = lpos_dat$Average.Rt.min.[i],
                                      Adduct.type = lpos_dat$Adduct.type[i],
                                      pb_prop = pb_prop,
                                      ib_prop = ib_prop,
                                      pb_flag = ifelse(pb_prop >= 10, 1, 0),
                                      ib_flag = ifelse(ib_prop >= 10, 1, 0))
}
blank_comparison <- Reduce("rbind", blank_comparison) %>%
  dplyr::mutate(join_flag = ifelse(pb_flag == 1 & ib_flag == 1,1,0),
                different_flag = ifelse(pb_flag == ib_flag, 0, 1))

blank_comparison_pos <- blank_comparison %>%
  dplyr::mutate(blank_remove = ifelse(pb_flag == 1 | ib_flag == 1, 1, 0)) %>%
  dplyr::select(-pb_prop, -ib_prop)

blank_comparison <- blank_comparison %>%
  dplyr::mutate(blank_remove = ifelse(pb_flag == 1 | ib_flag == 1, 1, 0)) %>%
  dplyr::select(Metabolite.name, Average.Rt.min., Adduct.type, blank_remove)

# As per JK 4/25/24...remove based on both pb and ib
lpos_dat <- lpos_dat %>%
  dplyr::left_join(blank_comparison) %>%
  dplyr::filter(blank_remove == 0) %>%
  dplyr::select(-blank_remove)
nrow(blank_comparison_pos) - nrow(lpos_dat)

rm(i, ib_means, ib_prop, pb_means, pb_prop, tempdat, blank_comparison)

writexl::write_xlsx(x = list(`Neg` = blank_comparison_neg,
                             `Pos` = blank_comparison_pos),
                    path = here("data", "processed", "ptrc_lipids_blank_flags.xlsx"))

# -------------------------------------------------------------------------

rm(blank_comparison_neg, blank_comparison_pos)

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

run_order_pos <- readxl::read_xlsx(here("data", data_filenames[grepl("run", data_filenames)]),
                                   sheet = "POS")

run_order_neg <- readxl::read_xlsx(here("data", data_filenames[grepl("run", data_filenames)]),
                                   sheet = "NEG")

# Neg ---------------------------------------------------------------------

# Create initial edata, emeta, and fdata objects
lneg_edat <- lneg_dat %>%
  dplyr::select(metname.with.adduct,
                dplyr::contains("PTRC"))
lneg_emet <- lneg_dat %>%
  dplyr::select(-dplyr::contains("PTRC"))

lneg_fdat <- lneg_dat_samplabels %>%
  dplyr::select(-(`...1`:`...12`)) %>%
  t(.) %>%
  data.frame(.) %>%
  setNames(c("Category", "SampleID")) %>%
  dplyr::mutate(SampleID_og = SampleID,
                SampleID = make.names(SampleID)) %>%
  dplyr::filter(!is.na(Category))

lneg_fdat <- lneg_fdat %>%
  dplyr::left_join(run_order_neg %>%
                     dplyr::select(DATA, `Run Order`) %>%
                     dplyr::rename(SampleID = DATA,
                                   run_order = `Run Order`))

# Check consistency of edat and fdat names
all(lneg_fdat$SampleID %in% colnames(lneg_edat)[-1])

# Pos ---------------------------------------------------------------------

# Create initial edata, emeta, and fdata objects
lpos_edat <- lpos_dat %>%
  dplyr::select(metname.with.adduct,
                dplyr::contains("PTRC"))
lpos_emet <- lpos_dat %>%
  dplyr::select(-dplyr::contains("PTRC"))

lpos_fdat <- lpos_dat_samplabels %>%
  dplyr::select(-(`...1`:`...12`)) %>%
  t(.) %>%
  data.frame(.) %>%
  setNames(c("Category", "SampleID")) %>%
  dplyr::mutate(SampleID_og = SampleID,
                SampleID = make.names(SampleID)) %>%
  dplyr::filter(!is.na(Category))

lpos_fdat <- lpos_fdat %>%
  dplyr::left_join(run_order_pos %>%
                     dplyr::select(DATA, `Run Order`) %>%
                     dplyr::rename(SampleID = DATA,
                                   run_order = `Run Order`))

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
#                   metname.with.adduct == "neg_Cer 34:0;O2|Cer 18:0;O2/16:0_[M-H]-__B" ~ 23.59,
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
                  metname.with.adduct == "neg_Cer 34:0;O2|Cer 18:0;O2/16:0_[M-H]-__B" ~ 23.59,
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
#                   metname.with.adduct == "pos_Cer 40:1;O2|Cer 18:1;O2/22:0_[M+H-H2O]+__C" ~ 26.19,
#                   metname.with.adduct == "pos_Cer 43:2;O2|Cer 19:1;O2/24:1_[M+H]+__A" ~ 26.36,
#                   metname.with.adduct == "pos_Cer 44:2;O2|Cer 18:1;O2/26:1_[M+H]+" ~ 26.90,
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
                  metname.with.adduct == "pos_Cer 40:1;O2|Cer 18:1;O2/22:0_[M+H-H2O]+__C" ~ 26.19,
                  metname.with.adduct == "pos_Cer 43:2;O2|Cer 19:1;O2/24:1_[M+H]+__A" ~ 26.36,
                  metname.with.adduct == "pos_Cer 44:2;O2|Cer 18:1;O2/26:1_[M+H]+" ~ 26.90,
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
                dplyr::contains("PTRC")) %>%
  dplyr::mutate(Metabolite.name = gsub("__AA", "", metname.with.setpair)) %>%
  dplyr::select(-metname.with.setpair) %>%
  dplyr::relocate(Metabolite.name) %>%
  purrrlyr::slice_rows("Metabolite.name") %>%
  purrrlyr::dmap(sum)

lneg_tempdat_pairs_meta <- lneg_tempdat %>%
  dplyr::filter(metname_dupcount == 2) %>%
  dplyr::select(-dplyr::contains("PTRC")) %>%
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
                dplyr::contains("PTRC")) %>%
  purrrlyr::slice_rows("metname.with.setpair") %>%
  purrrlyr::dmap(sum) %>%
  dplyr::rename(Metabolite.name = metname.with.setpair)

lneg_tempdat_mults_sum2_meta <- lneg_tempdat_mults %>%
  dplyr::filter(setpair != "") %>%
  dplyr::select(-dplyr::contains("PTRC")) %>%
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
                Reference.m.z       = as.character(Reference.m.z),
                ppm.error           = as.character(ppm.error),
                # Total.score         = as.character(Total.score),
                # Dot.product         = as.character(Dot.product),
                # Reverse.dot.product = as.character(Reverse.dot.product),
                S.N.average         = as.character(S.N.average),
                Std.dev = as.character(Std.dev),
                median = as.character(median),
                X.CV = as.character(X.CV)) %>%#,
  # Fill..              = as.character(Fill..)) %>%
  dplyr::bind_rows(lneg_tempdat_pairs %>%
                     dplyr::mutate(Alignment.ID        = as.character(Alignment.ID),
                                   Average.Rt.min.     = as.character(Average.Rt.min.),
                                   Average.Mz          = as.character(Average.Mz),
                                   Reference.m.z       = as.character(Reference.m.z),
                                   ppm.error           = as.character(ppm.error),
                                   # Total.score         = as.character(Total.score),
                                   # Dot.product         = as.character(Dot.product),
                                   # Reverse.dot.product = as.character(Reverse.dot.product),
                                   S.N.average         = as.character(S.N.average),
                                   Std.dev = as.character(Std.dev),
                                   median = as.character(median),
                                   X.CV = as.character(X.CV))) %>% #,
  # Fill..              = as.character(Fill..))) %>%
  dplyr::bind_rows(lneg_tempdat_mults_sum3 %>%
                     dplyr::mutate(Alignment.ID        = as.character(Alignment.ID),
                                   Average.Rt.min.     = as.character(Average.Rt.min.),
                                   Average.Mz          = as.character(Average.Mz),
                                   Reference.m.z       = as.character(Reference.m.z),
                                   ppm.error           = as.character(ppm.error),
                                   # Total.score         = as.character(Total.score),
                                   # Dot.product         = as.character(Dot.product),
                                   # Reverse.dot.product = as.character(Reverse.dot.product),
                                   S.N.average         = as.character(S.N.average),
                                   Std.dev = as.character(Std.dev),
                                   median = as.character(median),
                                   X.CV = as.character(X.CV))) %>% #,
  # Fill..              = as.character(Fill..))) %>%
  dplyr::mutate(metname.with.adduct = ifelse(is.na(metname.with.adduct), Metabolite.name, metname.with.adduct))

length(unique(lneg_dat_sum$metname.with.adduct)) == nrow(lneg_dat_sum)

# Create accompanying edata and emeta objects
lneg_edat_sum <- lneg_dat_sum %>%
  dplyr::select(metname.with.adduct, 
                dplyr::contains("PTRC"))
lneg_emet_sum <- lneg_dat_sum %>%
  dplyr::select(-dplyr::contains("PTRC"))

rm(lneg_tempdat_mults_sum1, lneg_tempdat_mults_sum2, lneg_tempdat_mults)

# Pos ---------------------------------------------------------------------

lpos_tempdat_pairs <- lpos_tempdat %>%
  dplyr::filter(metname_dupcount == 2) %>%
  dplyr::select(metname.with.setpair,
                dplyr::contains("PTRC")) %>%
  dplyr::mutate(Metabolite.name = gsub("__AA", "", metname.with.setpair)) %>%
  dplyr::select(-metname.with.setpair) %>%
  dplyr::relocate(Metabolite.name) %>%
  purrrlyr::slice_rows("Metabolite.name") %>%
  purrrlyr::dmap(sum)

lpos_tempdat_pairs_meta <- lpos_tempdat %>%
  dplyr::filter(metname_dupcount == 2) %>%
  dplyr::select(-dplyr::contains("PTRC")) %>%
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
                dplyr::contains("PTRC")) %>%
  purrrlyr::slice_rows("metname.with.setpair") %>%
  purrrlyr::dmap(sum) %>%
  dplyr::rename(Metabolite.name = metname.with.setpair)

lpos_tempdat_mults_sum2_meta <- lpos_tempdat_mults %>%
  dplyr::filter(setpair != "") %>%
  dplyr::select(-dplyr::contains("PTRC")) %>%
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
                Reference.m.z       = as.character(Reference.m.z),
                ppm.error           = as.character(ppm.error),
                BEAT.AML.RT         = as.character(BEAT.AML.RT),
                # Total.score         = as.character(Total.score),
                # Dot.product         = as.character(Dot.product),
                # Reverse.dot.product = as.character(Reverse.dot.product),
                S.N.average         = as.character(S.N.average),
                Std.dev = as.character(Std.dev),
                median = as.character(median),
                X.CV = as.character(X.CV)) %>% #,
  # Fill..              = as.character(Fill..)) %>%
  dplyr::bind_rows(lpos_tempdat_pairs %>%
                     dplyr::mutate(Alignment.ID        = as.character(Alignment.ID),
                                   Average.Rt.min.     = as.character(Average.Rt.min.),
                                   Average.Mz          = as.character(Average.Mz),
                                   Reference.m.z       = as.character(Reference.m.z),
                                   ppm.error           = as.character(ppm.error),
                                   BEAT.AML.RT         = as.character(BEAT.AML.RT),
                                   # Total.score         = as.character(Total.score),
                                   # Dot.product         = as.character(Dot.product),
                                   # Reverse.dot.product = as.character(Reverse.dot.product),
                                   S.N.average         = as.character(S.N.average),
                                   Std.dev = as.character(Std.dev),
                                   median = as.character(median),
                                   X.CV = as.character(X.CV))) %>% #,
  # Fill..              = as.character(Fill..))) %>%
  dplyr::bind_rows(lpos_tempdat_mults_sum3 %>%
                     dplyr::mutate(Alignment.ID        = as.character(Alignment.ID),
                                   Average.Rt.min.     = as.character(Average.Rt.min.),
                                   Average.Mz          = as.character(Average.Mz),
                                   Reference.m.z       = as.character(Reference.m.z),
                                   ppm.error           = as.character(ppm.error),
                                   BEAT.AML.RT         = as.character(BEAT.AML.RT),
                                   # Total.score         = as.character(Total.score),
                                   # Dot.product         = as.character(Dot.product),
                                   # Reverse.dot.product = as.character(Reverse.dot.product),
                                   S.N.average         = as.character(S.N.average),
                                   Std.dev = as.character(Std.dev),
                                   median = as.character(median),
                                   X.CV = as.character(X.CV))) %>% #,
  # Fill..              = as.character(Fill..))) %>%
  dplyr::mutate(metname.with.adduct = ifelse(is.na(metname.with.adduct), Metabolite.name, metname.with.adduct))

length(unique(lpos_dat_sum$metname.with.adduct)) == nrow(lpos_dat_sum)

# Create accompanying edata and emeta objects
lpos_edat_sum <- lpos_dat_sum %>%
  dplyr::select(metname.with.adduct, 
                dplyr::contains("PTRC"))
lpos_emet_sum <- lpos_dat_sum %>%
  dplyr::select(-dplyr::contains("PTRC"))

rm(lpos_tempdat_mults_sum1, lpos_tempdat_mults_sum2, lpos_tempdat_mults)

# -------------------------------------------------------------------------
## Choosing Higher Intensity Adduct (Choice) ----------------------------------------
# Neg ---------------------------------------------------------------------

# Filter out cases that were "unmatched", meaning that their retention times were
# different. These are the cases in lneg_tempdat_mults_sum3 with singular letters
# appended. 
lneg_tempdat <- lneg_tempdat %>%
  dplyr::filter(!(metname.with.adduct %in% lneg_tempdat_mults_sum3$Metabolite.name)) %>%
  dplyr::filter(!(metname.with.adduct %in% lneg_tempdat_pairs$Metabolite.name)) %>%
  dplyr::mutate(metname.with.setpair = gsub("_\\[M\\+CH3COO\\]-", "", metname.with.setpair),
                metname.with.setpair = gsub("_\\[M-H\\]-", "", metname.with.setpair))

# For each case, determine whether one adduct is uniformly higher in its intensities than another
# Note that we average in cases where there are duplicated name/adduct combos. 
temp_metnames <- unique(lneg_tempdat$metname.with.setpair)
tempres_adduct_neg <- vector("list", length = length(temp_metnames))
for(i in 1:length(temp_metnames)){
  tempdat <- lneg_tempdat %>%
    dplyr::filter(metname.with.setpair %in% temp_metnames[i]) 
  preres <- tempdat %>% 
    dplyr::select(metname.with.adduct, metname.with.setpair, Metabolite.name, Adduct.type)
  tempdat <- tempdat %>%
    dplyr::select(Adduct.type, 
                  dplyr::contains("PTRC")) %>%
    purrrlyr::slice_rows("Adduct.type") %>%
    purrrlyr::dmap(mean)
  
  higher_idx <- apply(tempdat[,-(1)], 2, function(x){ifelse(diff(x) == 0, NA, which.max(x))})
  # intens_diff <- apply(tempdat[,-(1:21)], 2, function(x){max(x)-max(x[-which.max(x)])})
  higher_idx <- higher_idx[!is.na(higher_idx)]
  higher_idx_adduct <- tempdat$Adduct.type[higher_idx]
  
  tempres_adduct_neg[[i]] <- data.frame(table(higher_idx_adduct)) %>%
    dplyr::mutate(Perc = round(Freq/length(higher_idx)*100,3))
  if(!("[M-H]-" %in% tempres_adduct_neg[[i]]$higher_idx_adduct)){
    tempres_adduct_neg[[i]] <- rbind.data.frame(tempres_adduct_neg[[i]],
                                                data.frame(higher_idx_adduct = "[M-H]-",
                                                           Freq = 0, Perc = 0))
  } else if(!("[M+CH3COO]-" %in% tempres_adduct_neg[[i]]$higher_idx_adduct)){
    tempres_adduct_neg[[i]] <- rbind.data.frame(tempres_adduct_neg[[i]],
                                                data.frame(higher_idx_adduct = "[M+CH3COO]-",
                                                           Freq = 0, Perc = 0))
  }
  tempres_adduct_neg[[i]] <- tempres_adduct_neg[[i]] %>%
    dplyr::rename(Adduct.type = higher_idx_adduct) %>%
    dplyr::right_join(preres, by = "Adduct.type")
}
rm(i, tempdat, preres, higher_idx, higher_idx_adduct, temp_metnames)

temp2 <- Reduce("rbind", tempres_adduct_neg) %>%
  dplyr::group_by(Adduct.type) %>%
  dplyr::arrange(Adduct.type, desc(Perc)) %>%
  dplyr::mutate(Metabolite.name = factor(metname.with.setpair, levels = unique(.$metname.with.setpair))) %>%
  dplyr::ungroup() 

png(filename = here("figures", "ptrc_lipids", "negmode_adduct_intensity_compare.png"), units = "in", width = 16, height = 7, res = 300)
set.seed(1)
ggplot(data = temp2, aes(x = Metabolite.name, y = Perc, color = Adduct.type)) + 
  geom_jitter(height = 0, width = 0.2, shape = 1) + theme_bw() + 
  ylab("Percentage of Samples Higher") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

lneg_tempdat <- lneg_tempdat %>%
  dplyr::left_join(temp2 %>% dplyr::select(-Adduct.type, -Metabolite.name, -metname.with.setpair), by = "metname.with.adduct") %>%
  dplyr::rename(Freq_higher = Freq,
                Perc_higher = Perc)

# Since there appears to be a biological trend with intensities among adduct duplicates,
# we decide which to retain based on this
lneg_tempdat_choice <- lneg_tempdat %>% 
  dplyr::group_by(metname.with.setpair) %>%
  dplyr::arrange(desc(Perc_higher)) %>%
  dplyr::mutate(duporder = dplyr::row_number()) %>%
  dplyr::relocate(duporder, Perc_higher, .after = metname.with.setpair) %>%
  dplyr::ungroup() %>%
  dplyr::filter(duporder == 1)

lneg_dat_choice <- lneg_dat %>%
  dplyr::filter(!(metname.with.adduct %in% 
                    c(lneg_tempdat$metname.with.adduct))) %>%
  dplyr::bind_rows(lneg_tempdat_choice)

# Create accompanying edata and emeta objects
lneg_edat_choice <- lneg_dat_choice %>%
  dplyr::select(metname.with.adduct, 
                dplyr::contains("PTRC"))
lneg_emet_choice <- lneg_dat_choice %>%
  dplyr::select(-dplyr::contains("PTRC"))

rm(lneg_tempdat, lneg_tempdat_choice, lneg_tempdat_mults_sum3, lneg_tempdat_pairs,
   temp2, lneg_dat_samplabels, tempres_adduct_neg)

# Pos ---------------------------------------------------------------------

# Filter out cases that were "unmatched", meaning that their retention times were
# different. 
lpos_tempdat <- lpos_tempdat %>%
  dplyr::filter(!(metname.with.adduct %in% lpos_tempdat_mults_sum3$Metabolite.name)) %>%
  dplyr::filter(!(metname.with.adduct %in% lpos_tempdat_pairs$Metabolite.name[grepl("\\[M", lpos_tempdat_pairs$Metabolite.name)])) %>%
  dplyr::mutate(metname.with.setpair = gsub("_\\[M\\+H-H2O\\]\\+", "", metname.with.setpair),
                metname.with.setpair = gsub("_\\[M\\+H\\]\\+", "", metname.with.setpair))

# For each case, determine whether one adduct is uniformly higher in its intensities than another
# Note that we average in cases where there are duplicated name/adduct combos. 
temp_metnames <- unique(lpos_tempdat$metname.with.setpair)
tempres_adduct_pos <- vector("list", length = length(temp_metnames))
for(i in 1:length(temp_metnames)){
  tempdat <- lpos_tempdat %>%
    dplyr::filter(metname.with.setpair %in% temp_metnames[i]) 
  preres <- tempdat %>% dplyr::select(metname.with.adduct, metname.with.setpair, Metabolite.name, Adduct.type)
  tempdat <- tempdat %>%
    dplyr::select(Adduct.type, 
                  dplyr::contains("PTRC")) %>%
    purrrlyr::slice_rows("Adduct.type") %>%
    purrrlyr::dmap(mean)
  
  higher_idx <- apply(tempdat[,-(1)], 2, function(x){ifelse(diff(x) == 0, NA, which.max(x))})
  # intens_diff <- apply(tempdat[,-(1:21)], 2, function(x){max(x)-max(x[-which.max(x)])})
  higher_idx <- higher_idx[!is.na(higher_idx)]
  higher_idx_adduct <- tempdat$Adduct.type[higher_idx]
  
  tempres_adduct_pos[[i]] <- data.frame(table(higher_idx_adduct)) %>%
    dplyr::mutate(Perc = round(Freq/length(higher_idx)*100,3))
  if(!("[M+H]+" %in% tempres_adduct_pos[[i]]$higher_idx_adduct)){
    tempres_adduct_pos[[i]] <- rbind.data.frame(tempres_adduct_pos[[i]],
                                                data.frame(higher_idx_adduct = "[M+H]+",
                                                           Freq = 0, Perc = 0))
  } else if(!("[M+H-H2O]+" %in% tempres_adduct_pos[[i]]$higher_idx_adduct)){
    tempres_adduct_pos[[i]] <- rbind.data.frame(tempres_adduct_pos[[i]],
                                                data.frame(higher_idx_adduct = "[M+H-H2O]+",
                                                           Freq = 0, Perc = 0))
  }
  tempres_adduct_pos[[i]] <- tempres_adduct_pos[[i]] %>%
    dplyr::rename(Adduct.type = higher_idx_adduct) %>%
    dplyr::right_join(preres, by = "Adduct.type")
}
rm(i, tempdat, preres, higher_idx, higher_idx_adduct, temp_metnames)

temp2 <- Reduce("rbind", tempres_adduct_pos) %>%
  dplyr::group_by(Adduct.type) %>%
  dplyr::arrange(Adduct.type, desc(Perc)) %>%
  dplyr::mutate(Metabolite.name = factor(metname.with.setpair, levels = unique(.$metname.with.setpair))) %>%
  dplyr::ungroup() 

png(filename = here("figures", "ptrc_lipids", "posmode_adduct_intensity_compare.png"), units = "in", width = 16, height = 7, res = 300)
set.seed(1)
ggplot(data = temp2, aes(x = Metabolite.name, y = Perc, color = Adduct.type)) + 
  geom_jitter(height = 0, width = 0.2, shape = 1) + theme_bw() + 
  ylab("Percentage of Samples Higher") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

lpos_tempdat <- lpos_tempdat %>%
  dplyr::left_join(temp2 %>% dplyr::select(-Adduct.type, -Metabolite.name, -metname.with.setpair), by = "metname.with.adduct") %>%
  dplyr::rename(Freq_higher = Freq,
                Perc_higher = Perc)

# Since there appears to be a biological trend with intensities among adduct duplicates,
# we decide which to retain based on this
lpos_tempdat_choice <- lpos_tempdat %>% 
  dplyr::group_by(metname.with.setpair) %>%
  dplyr::arrange(desc(Perc_higher)) %>%
  dplyr::mutate(duporder = dplyr::row_number()) %>%
  dplyr::relocate(duporder, Perc_higher, .after = metname.with.setpair) %>%
  dplyr::ungroup() %>%
  dplyr::filter(duporder == 1)

lpos_dat_choice <- lpos_dat %>%
  dplyr::filter(!(metname.with.adduct %in% 
                    c(lpos_tempdat$metname.with.adduct))) %>%
  dplyr::bind_rows(lpos_tempdat_choice)

# Create accompanying edata and emeta objects
lpos_edat_choice <- lpos_dat_choice %>%
  dplyr::select(metname.with.adduct, 
                dplyr::contains("PTRC"))
lpos_emet_choice <- lpos_dat_choice %>%
  dplyr::select(-dplyr::contains("PTRC"))


rm(lpos_tempdat, lpos_tempdat_choice, lpos_tempdat_mults_sum3, lpos_tempdat_pairs,
   temp2, lpos_dat_samplabels, tempres_adduct_pos)


# -------------------------------------------------------------------------
### Common Formatting of Modes ---------------------------------------------------------
## Differing Adduct Among Duplicates Ignored -------------------------------

# Make sure sample names are identical in pos and neg modes for when they are later combined
colnames(lneg_edat) %in% colnames(lpos_edat)
colnames(lpos_edat) %in% colnames(lneg_edat)

colnames(lneg_edat) <- gsub("_Neg_20Jun25_Glacier_WCSH315302", "", colnames(lneg_edat))
colnames(lneg_edat) <- gsub("PTRC_exp26_", "", colnames(lneg_edat))

colnames(lpos_edat) <- gsub("_Pos_14Jun25_Glacier_WCSH315302", "", colnames(lpos_edat))
colnames(lpos_edat) <- gsub("_POS_14Jun25_Glacier_WCSH315302", "", colnames(lpos_edat))
colnames(lpos_edat) <- gsub("PTRC_exp26_", "", colnames(lpos_edat))

# All pos mode samples are in neg mode
colnames(lpos_edat)[!(colnames(lpos_edat) %in% colnames(lneg_edat))] 

# All neg mode samples are in pos mode
colnames(lneg_edat)[!(colnames(lneg_edat) %in% colnames(lpos_edat))] 

# Apply above changes in colnames to SampleIDs in the respective fdatasets
lneg_fdat <- lneg_fdat %>%
  dplyr::mutate(SampleID = gsub("_Neg_20Jun25_Glacier_WCSH315302", "", SampleID),
                SampleID = gsub("PTRC_exp26_", "", SampleID))

all(lneg_fdat$SampleID %in% colnames(lneg_edat)[-1])

lpos_fdat <- lpos_fdat %>%
  dplyr::mutate(SampleID = gsub("_Pos_14Jun25_Glacier_WCSH315302", "", SampleID),
                SampleID = gsub("_POS_14Jun25_Glacier_WCSH315302", "", SampleID),
                SampleID = gsub("PTRC_exp26_", "", SampleID))

all(lpos_fdat$SampleID %in% colnames(lpos_edat)[-1])

## Sum ---------------------------------------------------------------------

# Make sure sample names are identical in pos and neg modes for when they are later combined
colnames(lneg_edat_sum) %in% colnames(lpos_edat_sum)
colnames(lpos_edat_sum) %in% colnames(lneg_edat_sum)

colnames(lneg_edat_sum) <- gsub("_Neg_20Jun25_Glacier_WCSH315302", "", colnames(lneg_edat_sum))
colnames(lneg_edat_sum) <- gsub("PTRC_exp26_", "", colnames(lneg_edat_sum))

colnames(lpos_edat_sum) <- gsub("_Pos_14Jun25_Glacier_WCSH315302", "", colnames(lpos_edat_sum))
colnames(lpos_edat_sum) <- gsub("_POS_14Jun25_Glacier_WCSH315302", "", colnames(lpos_edat_sum))
colnames(lpos_edat_sum) <- gsub("PTRC_exp26_", "", colnames(lpos_edat_sum))

# All pos mode samples are in neg mode
colnames(lpos_edat_sum)[!(colnames(lpos_edat_sum) %in% colnames(lneg_edat_sum))] 

# All neg mode samples are in pos mode
colnames(lneg_edat_sum)[!(colnames(lneg_edat_sum) %in% colnames(lpos_edat_sum))] 


all(lneg_fdat$SampleID %in% colnames(lneg_edat_sum)[-1])

all(lpos_fdat$SampleID %in% colnames(lpos_edat_sum)[-1])

## Choice ------------------------------------------------------------------

# Make sure sample names are identical in pos and neg modes for when they are later combined
colnames(lneg_edat_choice) %in% colnames(lpos_edat_choice)
colnames(lpos_edat_choice) %in% colnames(lneg_edat_choice)

colnames(lneg_edat_choice) <- gsub("_Neg_20Jun25_Glacier_WCSH315302", "", colnames(lneg_edat_choice))
colnames(lneg_edat_choice) <- gsub("PTRC_exp26_", "", colnames(lneg_edat_choice))

colnames(lpos_edat_choice) <- gsub("_Pos_14Jun25_Glacier_WCSH315302", "", colnames(lpos_edat_choice))
colnames(lpos_edat_choice) <- gsub("_POS_14Jun25_Glacier_WCSH315302", "", colnames(lpos_edat_choice))
colnames(lpos_edat_choice) <- gsub("PTRC_exp26_", "", colnames(lpos_edat_choice))

# All pos mode samples are in neg mode
colnames(lpos_edat_choice)[!(colnames(lpos_edat_choice) %in% colnames(lneg_edat_choice))] 

# All neg mode samples are in pos mode
colnames(lneg_edat_choice)[!(colnames(lneg_edat_choice) %in% colnames(lpos_edat_choice))] 


all(lneg_fdat$SampleID %in% colnames(lneg_edat_choice)[-1])

all(lpos_fdat$SampleID %in% colnames(lpos_edat_choice)[-1])

# -------------------------------------------------------------------------

### Process Blanks Assessment -----------------------------------------------

## Negative Mode -----------------------------------------------------------

# We take the average values across the blanks and compare to the 
# medacross the analytic samples 
pb_means <- apply(lneg_edat_sum %>% 
                    dplyr::select(dplyr::contains("PB")), 1, 
                  function(x){mean(as.numeric(x), na.rm = TRUE)})
samp_means <- apply(lneg_edat_sum %>% 
                      dplyr::select(-dplyr::contains("PB"),
                                    -dplyr::contains("QC"),
                                    -dplyr::contains("Blank"),
                                    -metname.with.adduct), 1, 
                    function(x){mean(as.numeric(x), na.rm = TRUE)})

# Graphic to visualize all pointwise differences (PB1)
tempdf_pb1 <- lneg_edat_sum %>% 
  dplyr::select(-dplyr::contains("PB"),
                -dplyr::contains("QC"),
                -dplyr::contains("Blank"))
for(i in 2:ncol(tempdf_pb1)){
  # Want to remove the cases where both are 0
  temp <- mapply(function(x,y){
    if(x == 0 & y == 0){
      return(NA)
    } else{
      x-y
    }
  }, tempdf_pb1[,i, drop = TRUE], lneg_edat_sum[,"Lip_PB_01",drop = TRUE])
  tempdf_pb1[,i] <- temp
}
tempdf_pb1 <- tempdf_pb1 %>%
  tidyr::pivot_longer(!metname.with.adduct, names_to = "Sample", values_to = "Difference") %>%
  as.data.frame() %>%
  dplyr::filter(!is.na(Difference)) %>%
  dplyr::group_by(metname.with.adduct) %>%
  dplyr::mutate(idx = dplyr::row_number()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Diffcat = cut(Difference, breaks = c(min(Difference), 0, 5000, max(Difference)), include.lowest = TRUE))

tempdf_pb1_counts <- tempdf_pb1 %>%
  dplyr::group_by(metname.with.adduct, Diffcat) %>%
  dplyr::summarise(count = dplyr::n())

# png(filename = here("reports", "final", "report_figures", "kftall_lipids_worm", "pb1_negmode_pointdiff.png"), units = "in", width = 16, height = 5, res = 300)
neg_p1 <- ggplot(data = tempdf_pb1_counts, aes(fill = Diffcat, y = count, x = metname.with.adduct)) +
  geom_bar(position = "stack", stat = "identity") + theme(axis.text.x = element_blank()) +
  ylab("Samples") + ggtitle("Intensity Differences: PB 01 (Neg Mode)") + xlab("Metabolite")
# dev.off()

# Export cases where there is a negative diffcat
tempdf_pb1_export <- tempdf_pb1 %>%
  dplyr::filter(Difference < 0) %>%
  dplyr::select(-Diffcat)

# Graphic to visualize all pointwise differences (PB2)
tempdf_pb2 <- lneg_edat_sum %>% 
  dplyr::select(-dplyr::contains("PB"),
                -dplyr::contains("QC"),
                -dplyr::contains("Blank"))
for(i in 2:ncol(tempdf_pb2)){
  # Want to remove the cases where both are 0
  temp <- mapply(function(x,y){
    if(x == 0 & y == 0){
      return(NA)
    } else{
      x-y
    }
  }, tempdf_pb2[,i, drop = TRUE], lneg_edat_sum[,"Lip_PB_02",drop = TRUE])
  tempdf_pb2[,i] <- temp
}
tempdf_pb2 <- tempdf_pb2 %>%
  tidyr::pivot_longer(!metname.with.adduct, names_to = "Sample", values_to = "Difference") %>%
  as.data.frame() %>%
  dplyr::filter(!is.na(Difference)) %>%
  dplyr::group_by(metname.with.adduct) %>%
  dplyr::mutate(idx = dplyr::row_number()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Diffcat = cut(Difference, breaks = c(min(Difference), 0, 5000, max(Difference)), include.lowest = TRUE))

tempdf_pb2_counts <- tempdf_pb2 %>%
  dplyr::group_by(metname.with.adduct, Diffcat) %>%
  dplyr::summarise(count = dplyr::n())

# png(filename = here("reports", "final", "report_figures", "kftall_lipids_worm", "pb2_negmode_pointdiff.png"), units = "in", width = 16, height = 5, res = 300)
neg_p2 <- ggplot(data = tempdf_pb2_counts, aes(fill = Diffcat, y = count, x = metname.with.adduct)) +
  geom_bar(position = "stack", stat = "identity") + theme(axis.text.x = element_blank()) +
  ylab("Samples") + ggtitle("Intensity Differences: PB 02 (Neg Mode)") + xlab("Metabolite")
# dev.off()

# Export cases where there is a negative diffcat
tempdf_pb2_export <- tempdf_pb2 %>%
  dplyr::filter(Difference < 0) %>%
  dplyr::select(-Diffcat)

# Graphic to visualize all pointwise differences (PB3)
tempdf_pb3 <- lneg_edat_sum %>% 
  dplyr::select(-dplyr::contains("PB"),
                -dplyr::contains("QC"),
                -dplyr::contains("Blank"))
for(i in 2:ncol(tempdf_pb3)){
  # Want to remove the cases where both are 0
  temp <- mapply(function(x,y){
    if(x == 0 & y == 0){
      return(NA)
    } else{
      x-y
    }
  }, tempdf_pb3[,i, drop = TRUE], lneg_edat_sum[,"Lip_PB_03",drop = TRUE])
  tempdf_pb3[,i] <- temp
}
tempdf_pb3 <- tempdf_pb3 %>%
  tidyr::pivot_longer(!metname.with.adduct, names_to = "Sample", values_to = "Difference") %>%
  as.data.frame() %>%
  dplyr::filter(!is.na(Difference)) %>%
  dplyr::group_by(metname.with.adduct) %>%
  dplyr::mutate(idx = dplyr::row_number()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Diffcat = cut(Difference, breaks = c(min(Difference), 0, 5000, max(Difference)), include.lowest = TRUE))

tempdf_pb3_counts <- tempdf_pb3 %>%
  dplyr::group_by(metname.with.adduct, Diffcat) %>%
  dplyr::summarise(count = dplyr::n())

# png(filename = here("reports", "final", "report_figures", "kftall_lipids_worm", "pb3_negmode_pointdiff.png"), units = "in", width = 16, height = 5, res = 300)
neg_p3 <- ggplot(data = tempdf_pb3_counts, aes(fill = Diffcat, y = count, x = metname.with.adduct)) +
  geom_bar(position = "stack", stat = "identity") + theme(axis.text.x = element_blank()) +
  ylab("Samples") + ggtitle("Intensity Differences: PB 03 (Neg Mode)") + xlab("Metabolite")
# dev.off()

# Export cases where there is a negative diffcat
tempdf_pb3_export <- tempdf_pb3 %>%
  dplyr::filter(Difference < 0) %>%
  dplyr::select(-Diffcat)

tempdf_pb1_export_neg <- tempdf_pb1_export
tempdf_pb2_export_neg <- tempdf_pb2_export
tempdf_pb3_export_neg <- tempdf_pb3_export

## Positive Mode -----------------------------------------------------------

# We take the average values across the blanks and compare to the 
# medacross the analytic samples 
pb_means <- apply(lpos_edat_sum %>% 
                    dplyr::select(dplyr::contains("PB")), 1, 
                  function(x){mean(as.numeric(x), na.rm = TRUE)})
samp_means <- apply(lpos_edat_sum %>% 
                      dplyr::select(-dplyr::contains("PB"),
                                    -dplyr::contains("QC"),
                                    -dplyr::contains("Blank"),
                                    -metname.with.adduct), 1, 
                    function(x){mean(as.numeric(x), na.rm = TRUE)})

# Graphic to visualize all pointwise differences (PB1)
tempdf_pb1 <- lpos_edat_sum %>% 
  dplyr::select(-dplyr::contains("PB"),
                -dplyr::contains("QC"),
                -dplyr::contains("Blank"))
for(i in 2:ncol(tempdf_pb1)){
  # Want to remove the cases where both are 0
  temp <- mapply(function(x,y){
    if(x == 0 & y == 0){
      return(NA)
    } else{
      x-y
    }
  }, tempdf_pb1[,i, drop = TRUE], lpos_edat_sum[,"Lip_PB_01",drop = TRUE])
  tempdf_pb1[,i] <- temp
}
tempdf_pb1 <- tempdf_pb1 %>%
  tidyr::pivot_longer(!metname.with.adduct, names_to = "Sample", values_to = "Difference") %>%
  as.data.frame() %>%
  dplyr::filter(!is.na(Difference)) %>%
  dplyr::group_by(metname.with.adduct) %>%
  dplyr::mutate(idx = dplyr::row_number()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Diffcat = cut(Difference, breaks = c(min(Difference), 0, 5000, max(Difference)), include.lowest = TRUE))

tempdf_pb1_counts <- tempdf_pb1 %>%
  dplyr::group_by(metname.with.adduct, Diffcat) %>%
  dplyr::summarise(count = dplyr::n())

# png(filename = here("reports", "final", "report_figures", "kftall_lipids_worm", "pb1_posmode_pointdiff.png"), units = "in", width = 16, height = 5, res = 300)
pos_p1 <- ggplot(data = tempdf_pb1_counts, aes(fill = Diffcat, y = count, x = metname.with.adduct)) +
  geom_bar(position = "stack", stat = "identity") + theme(axis.text.x = element_blank()) +
  ylab("Samples") + ggtitle("Intensity Differences: PB 01 (Pos Mode)") + xlab("Metabolite")
# dev.off()

# Export cases where there is a positive diffcat
tempdf_pb1_export <- tempdf_pb1 %>%
  dplyr::filter(Difference < 0) %>%
  dplyr::select(-Diffcat)

# Graphic to visualize all pointwise differences (PB2)
tempdf_pb2 <- lpos_edat_sum %>% 
  dplyr::select(-dplyr::contains("PB"),
                -dplyr::contains("QC"),
                -dplyr::contains("Blank"))
for(i in 2:ncol(tempdf_pb2)){
  # Want to remove the cases where both are 0
  temp <- mapply(function(x,y){
    if(x == 0 & y == 0){
      return(NA)
    } else{
      x-y
    }
  }, tempdf_pb2[,i, drop = TRUE], lpos_edat_sum[,"Lip_PB_02",drop = TRUE])
  tempdf_pb2[,i] <- temp
}
tempdf_pb2 <- tempdf_pb2 %>%
  tidyr::pivot_longer(!metname.with.adduct, names_to = "Sample", values_to = "Difference") %>%
  as.data.frame() %>%
  dplyr::filter(!is.na(Difference)) %>%
  dplyr::group_by(metname.with.adduct) %>%
  dplyr::mutate(idx = dplyr::row_number()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Diffcat = cut(Difference, breaks = c(min(Difference), 0, 5000, max(Difference)), include.lowest = TRUE))

tempdf_pb2_counts <- tempdf_pb2 %>%
  dplyr::group_by(metname.with.adduct, Diffcat) %>%
  dplyr::summarise(count = dplyr::n())

# png(filename = here("reports", "final", "report_figures", "kftall_lipids_worm", "pb2_posmode_pointdiff.png"), units = "in", width = 16, height = 5, res = 300)
pos_p2 <- ggplot(data = tempdf_pb2_counts, aes(fill = Diffcat, y = count, x = metname.with.adduct)) +
  geom_bar(position = "stack", stat = "identity") + theme(axis.text.x = element_blank()) +
  ylab("Samples") + ggtitle("Intensity Differences: PB 02 (Pos Mode)") + xlab("Metabolite")
# dev.off()

# Export cases where there is a positive diffcat
tempdf_pb2_export <- tempdf_pb2 %>%
  dplyr::filter(Difference < 0) %>%
  dplyr::select(-Diffcat)

# Graphic to visualize all pointwise differences (PB3)
tempdf_pb3 <- lpos_edat_sum %>% 
  dplyr::select(-dplyr::contains("PB"),
                -dplyr::contains("QC"),
                -dplyr::contains("Blank"))
for(i in 2:ncol(tempdf_pb3)){
  # Want to remove the cases where both are 0
  temp <- mapply(function(x,y){
    if(x == 0 & y == 0){
      return(NA)
    } else{
      x-y
    }
  }, tempdf_pb3[,i, drop = TRUE], lpos_edat_sum[,"Lip_PB_03",drop = TRUE])
  tempdf_pb3[,i] <- temp
}
tempdf_pb3 <- tempdf_pb3 %>%
  tidyr::pivot_longer(!metname.with.adduct, names_to = "Sample", values_to = "Difference") %>%
  as.data.frame() %>%
  dplyr::filter(!is.na(Difference)) %>%
  dplyr::group_by(metname.with.adduct) %>%
  dplyr::mutate(idx = dplyr::row_number()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Diffcat = cut(Difference, breaks = c(min(Difference), 0, 5000, max(Difference)), include.lowest = TRUE))

tempdf_pb3_counts <- tempdf_pb3 %>%
  dplyr::group_by(metname.with.adduct, Diffcat) %>%
  dplyr::summarise(count = dplyr::n())

# png(filename = here("reports", "final", "report_figures", "kftall_lipids_worm", "pb3_posmode_pointdiff.png"), units = "in", width = 16, height = 5, res = 300)
pos_p3 <- ggplot(data = tempdf_pb3_counts, aes(fill = Diffcat, y = count, x = metname.with.adduct)) +
  geom_bar(position = "stack", stat = "identity") + theme(axis.text.x = element_blank()) +
  ylab("Samples") + ggtitle("Intensity Differences: PB 03 (Pos Mode)") + xlab("Metabolite")
# dev.off()

# Export cases where there is a positive diffcat
tempdf_pb3_export <- tempdf_pb3 %>%
  dplyr::filter(Difference < 0) %>%
  dplyr::select(-Diffcat)

tempdf_pb1_export_pos <- tempdf_pb1_export
tempdf_pb2_export_pos <- tempdf_pb2_export
tempdf_pb3_export_pos <- tempdf_pb3_export


# -------------------------------------------------------------------------

png(filename = here("figures", "ptrc_lipids", "pb_negmode_pointdiff.png"), units = "in", width = 16, height = 10, res = 300)
ggpubr::ggarrange(neg_p1, neg_p2, neg_p3, nrow = 3, ncol = 1, common.legend = TRUE,
                  legend = "right")
dev.off()

png(filename = here("figures", "ptrc_lipids", "pb_posmode_pointdiff.png"), units = "in", width = 16, height = 10, res = 300)
ggpubr::ggarrange(pos_p1, pos_p2, pos_p3, nrow = 3, ncol = 1, common.legend = TRUE,
                  legend = "right")
dev.off()

pb_pointdiffs_imgs <- c(here("figures", "ptrc_lipids", "pb_negmode_pointdiff.png"),
                        here("figures", "ptrc_lipids", "pb_posmode_pointdiff.png"))

rm(pb_pointdiffs_imgs)

# -------------------------------------------------------------------------


### Flagging ----------------------------------------------------------------
## Sum ---------------------------------------------------------------------

# Sum: Check if any of the lipids are flagged for the QC sample and sample intensity check
qcsamp_check <- apply(lneg_dat_sum %>%
                        dplyr::select(dplyr::contains("QC_Pool")),
                      1, function(x){sum(x<3000)/length(x)*100})
lneg_emet_sum$QCpool_flag <- as.numeric(qcsamp_check >= 40) 

samp_check <- apply(lneg_dat_sum %>%
                      dplyr::select(-dplyr::contains("QC_Pool"),
                                    -dplyr::contains("PB"),
                                    -dplyr::contains("Blank")) %>%
                      dplyr::select(dplyr::contains("PTRC")),
                    1, function(x){sum(x<3000)/length(x)*100})
lneg_emet_sum$Sample_flag <- as.numeric(samp_check >= 40) 
rm(qcsamp_check, samp_check)

qcsamp_check <- apply(lpos_dat_sum %>%
                        dplyr::select(dplyr::contains("QC_Pool")),
                      1, function(x){sum(x<5000)/length(x)*100})
lpos_emet_sum$QCpool_flag <- as.numeric(qcsamp_check >= 40) 

samp_check <- apply(lpos_dat_sum %>%
                      dplyr::select(-dplyr::contains("QC_Pool"),
                                    -dplyr::contains("PB"),
                                    -dplyr::contains("Blank")) %>%
                      dplyr::select(dplyr::contains("PTRC")),
                    1, function(x){sum(x<5000)/length(x)*100})
lpos_emet_sum$Sample_flag <- as.numeric(samp_check >= 40) 
rm(qcsamp_check, samp_check)


## Choice ------------------------------------------------------------------

# Choice/Higher: Check if any of the lipids are flagged for the QC sample and sample intensity check
qcsamp_check <- apply(lneg_dat_choice %>%
                        dplyr::select(dplyr::contains("QC_Pool")),
                      1, function(x){sum(x<3000)/length(x)*100})
lneg_emet_choice$QCpool_flag <- as.numeric(qcsamp_check >= 40) 

samp_check <- apply(lneg_dat_choice %>%
                      dplyr::select(-dplyr::contains("QC_Pool"),
                                    -dplyr::contains("PB"),
                                    -dplyr::contains("Blank")) %>%
                      dplyr::select(dplyr::contains("PTRC")),
                    1, function(x){sum(x<3000)/length(x)*100})
lneg_emet_choice$Sample_flag <- as.numeric(samp_check >= 40) 
rm(qcsamp_check, samp_check)

qcsamp_check <- apply(lpos_dat_choice %>%
                        dplyr::select(dplyr::contains("QC_Pool")),
                      1, function(x){sum(x<5000)/length(x)*100})
lpos_emet_choice$QCpool_flag <- as.numeric(qcsamp_check >= 40) 

samp_check <- apply(lpos_dat_choice %>%
                      dplyr::select(-dplyr::contains("QC_Pool"),
                                    -dplyr::contains("PB"),
                                    -dplyr::contains("Blank")) %>%
                      dplyr::select(dplyr::contains("PTRC")),
                    1, function(x){sum(x<5000)/length(x)*100})
lpos_emet_choice$Sample_flag <- as.numeric(samp_check >= 40) 
rm(qcsamp_check, samp_check)



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
                      dplyr::select(dplyr::contains("QC_Pool")),
                    1, function(x){sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE)*100})
lneg_emet_sum$QC_CV <- as.numeric(qccv_check) 
rm(qccv_check)

qccv_check <- apply(lpos_dat_sum %>%
                      dplyr::select(dplyr::contains("QC_Pool")),
                    1, function(x){sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE)*100})
lpos_emet_sum$QC_CV <- as.numeric(qccv_check) 
rm(qccv_check)


## Choice ------------------------------------------------------------------

# Choice/Higher: Check if any of the lipids are flagged for the QC sample and sample intensity check
qccv_check <- apply(lneg_dat_choice %>%
                      dplyr::select(dplyr::contains("QC_Pool")),
                    1, function(x){sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE)*100})
lneg_emet_choice$QC_CV <- as.numeric(qccv_check) 
rm(qccv_check)

qccv_check <- apply(lpos_dat_choice %>%
                      dplyr::select(dplyr::contains("QC_Pool")),
                    1, function(x){sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE)*100})
lpos_emet_choice$QC_CV <- as.numeric(qccv_check) 
rm(qccv_check)


# Visualizations ----------------------------------------------------------

png(filename = here("figures", "ptrc_lipids", "cv_distribution_neg.png"), units = "in", width = 16, height = 5, res = 300)
ggplot(data = lneg_emet_sum, aes(x = QC_CV)) + 
  geom_histogram(bins = 100, color = "black", fill = "gray") + 
  theme_bw() + 
  geom_vline(xintercept = median(lneg_emet_sum$QC_CV, na.rm = TRUE), 
             linetype = "dashed", color = "red", linewidth = 1.2) + 
  annotate("text", x = Inf, y = Inf, label = paste("Median %CV: ", 
                                                   round(median(lneg_emet_sum$QC_CV, na.rm = TRUE), 3)), vjust = 2, hjust = 1.2)
dev.off()

png(filename = here("figures", "ptrc_lipids", "cv_distribution_pos.png"), units = "in", width = 16, height = 5, res = 300)
ggplot(data = lpos_emet_sum, aes(x = QC_CV)) + 
  geom_histogram(bins = 100, color = "black", fill = "gray") + 
  theme_bw() + 
  geom_vline(xintercept = median(lpos_emet_sum$QC_CV, na.rm = TRUE), 
             linetype = "dashed", color = "red", linewidth = 1.2) + 
  annotate("text", x = Inf, y = Inf, label = paste("Median %CV: ", 
                                                   round(median(lpos_emet_sum$QC_CV, na.rm = TRUE), 3)), vjust = 2, hjust = 1.2)
dev.off()

# -------------------------------------------------------------------------


### Log2-transform / 0-to-NA ------------------------------------------------
## Sum ---------------------------------------------------------------------

# 10898 instances of 0 have been replaced with NA
pm_lneg_sum <- as.lipidData(e_data = lneg_edat_sum, 
                            f_data = lneg_fdat, 
                            e_meta = lneg_emet_sum,
                            edata_cname = "metname.with.adduct", 
                            fdata_cname = "SampleID",
                            emeta_cname = "metname.with.adduct",
                            data_scale = "abundance",
                            data_types = "Negative")


# 12756 instances of 0 have been replaced with NA
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

## Choice ------------------------------------------------------------------

# 10965 instances of 0 have been replaced with NA
pm_lneg_choice <- as.lipidData(e_data = lneg_edat_choice, 
                               f_data = lneg_fdat, 
                               e_meta = lneg_emet_choice,
                               edata_cname = "metname.with.adduct", 
                               fdata_cname = "SampleID",
                               emeta_cname = "metname.with.adduct",
                               data_scale = "abundance",
                               data_types = "Negative")


# 12903 instances of 0 have been replaced with NA
pm_lpos_choice <- as.lipidData(e_data = lpos_edat_choice, 
                               f_data = lpos_fdat, 
                               e_meta = lpos_emet_choice,
                               edata_cname = "metname.with.adduct", 
                               fdata_cname = "SampleID",
                               emeta_cname = "metname.with.adduct",
                               data_scale = "abundance",
                               data_types = "Positive")


# Transform data to the log2 scale
pm_lneg_choice <- edata_transform(pm_lneg_choice, data_scale = "log2")
pm_lpos_choice <- edata_transform(pm_lpos_choice, data_scale = "log2")

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
                    path = here("data", "processed", "ptrc_lipids_log2_sum.xlsx"))
saveRDS(lneg_sum_export_dat, here("data", "processed", "ptrc_lipids_neg_log2_sum.rds"))
saveRDS(lpos_sum_export_dat, here("data", "processed", "ptrc_lipids_pos_log2_sum.rds"))

# Choice ------------------------------------------------------------------

lneg_choice_export_dat <- pm_lneg_choice$e_data %>%
  dplyr::left_join(pm_lneg_choice$e_meta) %>%
  dplyr::relocate(Alignment.ID:QC_CV, .after = metname.with.adduct) %>%
  dplyr::select(-(dupcount:dupletter), -(metname_dupcount:Freq_higher)) %>%
  dplyr::mutate(Metabolite.name = gsub("neg_", "", Metabolite.name),
                metname.with.adduct = gsub("neg_", "", metname.with.adduct)) %>% #,
  # metname.with.adduct = gsub("_\\[M-H\\]-", "", metname.with.adduct),
  # metname.with.adduct = gsub("_\\[M\\+CH3COO\\]-", "", metname.with.adduct)) %>%
  dplyr::rename(original_name = Metabolite.name,
                formatted_name = metname.with.adduct)

lpos_choice_export_dat <- pm_lpos_choice$e_data %>%
  dplyr::left_join(pm_lpos_choice$e_meta) %>%
  dplyr::relocate(Alignment.ID:QC_CV, .after = metname.with.adduct) %>%
  dplyr::select(-(dupcount:dupletter), -(metname_dupcount:Freq_higher)) %>%
  dplyr::mutate(Metabolite.name = gsub("pos_", "", Metabolite.name),
                metname.with.adduct = gsub("pos_", "", metname.with.adduct)) %>% #,
  # metname.with.adduct = gsub("_\\[M\\+H-H2O\\]\\+", "", metname.with.adduct),
  # metname.with.adduct = gsub("_\\[M\\+H\\]\\+", "", metname.with.adduct),
  # metname.with.adduct = gsub("_\\[M\\+Na\\]\\+", "", metname.with.adduct),
  # metname.with.adduct = gsub("_\\[M\\+NH4\\]\\+", "", metname.with.adduct)) %>%
  dplyr::rename(original_name = Metabolite.name,
                formatted_name = metname.with.adduct)

writexl::write_xlsx(x = list(`Negative` = lneg_choice_export_dat ,
                             `Positive` = lpos_choice_export_dat),
                    path = here("data", "processed", "ptrc_lipids_log2_choice.xlsx"))
saveRDS(lneg_choice_export_dat, here("data", "processed", "ptrc_lipids_neg_log2_choice.rds"))
saveRDS(lpos_choice_export_dat, here("data", "processed", "ptrc_lipids_pos_log2_choice.rds"))

# -------------------------------------------------------------------------

### Run Order Correction ----------------------------------------------------

run_order_neg <- lneg_fdat %>%
  dplyr::select(SampleID, run_order)
run_order_pos <- lpos_fdat %>%
  dplyr::select(SampleID, run_order)

## Sum ---------------------------------------------------------------------
# Processing -------------------------------------------------------------

run_order_neg$SampleID[!(run_order_neg$SampleID %in% names(lneg_sum_export_dat))]

names(lneg_sum_export_dat)[!(names(lneg_sum_export_dat) %in% make.names(run_order_neg$SampleID))]

# Extract QC samples to be used for loess curve fitting (i.e. every other QC sample)

lneg_qc_train <- lneg_sum_export_dat %>%
  dplyr::select(formatted_name:QC_CV, QC_Pool_Lip_001,
                QC_Pool_Lip_003, QC_Pool_Lip_01, QC_Pool_Lip_03,
                QC_Pool_Lip_05, QC_Pool_Lip_07)
lneg_qc_val <- lneg_sum_export_dat %>%
  dplyr::select(formatted_name:QC_CV, QC_Pool_Lip_002,
                QC_Pool_Lip_004, QC_Pool_Lip_02, QC_Pool_Lip_04,
                QC_Pool_Lip_06, QC_Pool_Lip_08)


run_order_pos$SampleID[!(run_order_pos$SampleID %in% names(lpos_sum_export_dat))]

names(lpos_sum_export_dat)[!(names(lpos_sum_export_dat) %in% make.names(run_order_pos$SampleID))]

# Extract QC samples to be used for loess curve fitting (i.e. every other QC sample)

lpos_qc_train <- lpos_sum_export_dat %>%
  dplyr::select(formatted_name:QC_CV, QC_Pool_Lip_001,
                QC_Pool_Lip_003, QC_Pool_Lip_01, QC_Pool_Lip_03,
                QC_Pool_Lip_05, QC_Pool_Lip_07)
lpos_qc_val <- lpos_sum_export_dat %>%
  dplyr::select(formatted_name:QC_CV, QC_Pool_Lip_002,
                QC_Pool_Lip_004, QC_Pool_Lip_02, QC_Pool_Lip_04,
                QC_Pool_Lip_06, QC_Pool_Lip_08)

# Apply correction/normalization (to qc and samples) ------------------------------------------

lneg_runorder_tc <- lneg_sum_export_dat
for(i in 1:nrow(lneg_runorder_tc)){
  
  tempdat <- lneg_runorder_tc[i,] %>%
    dplyr::select(dplyr::contains("Lip_"),
                  -dplyr::contains("PB"),
                  -dplyr::contains("Blank")) %>%
    t(.) %>%
    as.data.frame(.) %>%
    setNames("Abundance") %>%
    dplyr::mutate(SampleID = row.names(.)) %>%
    dplyr::left_join(run_order_neg)
  
  tempfitdat <- lneg_qc_train[i,] %>%
    dplyr::select(-dplyr::contains("PB"),
                  -dplyr::contains("Blank"),
                  -(formatted_name:QC_CV)) %>%
    t(.) %>%
    as.data.frame(.) %>%
    setNames("Abundance") %>%
    dplyr::mutate(SampleID = row.names(.)) %>%
    dplyr::left_join(run_order_neg)
  
  tempmod <- lm(Abundance ~ run_order, data = tempfitdat)
  predvals <- predict(tempmod, tempdat)
  
  # all(names(lneg_runorder_tc[i, 47:166]) == tempdat$SampleID)
  lneg_runorder_tc[i, 47:166] <- lneg_runorder_tc[i, 47:166] - predvals + median(tempdat %>%
                                                                                   dplyr::filter(!grepl("PB", SampleID)) %>%
                                                                                   dplyr::filter(!grepl("Blankk", SampleID)) %>%
                                                                                   .$Abundance, na.rm = TRUE)
}

lpos_runorder_tc <- lpos_sum_export_dat
for(i in 1:nrow(lpos_runorder_tc)){
  
  tempdat <- lpos_runorder_tc[i,] %>%
    dplyr::select(dplyr::contains("Lip_"),
                  -dplyr::contains("PB"),
                  -dplyr::contains("Blank")) %>%
    t(.) %>%
    as.data.frame(.) %>%
    setNames("Abundance") %>%
    dplyr::mutate(SampleID = row.names(.)) %>%
    dplyr::left_join(run_order_pos)
  
  tempfitdat <- lpos_qc_train[i,] %>%
    dplyr::select(-dplyr::contains("PB"),
                  -dplyr::contains("Blank"),
                  -(formatted_name:QC_CV)) %>%
    t(.) %>%
    as.data.frame(.) %>%
    setNames("Abundance") %>%
    dplyr::mutate(SampleID = row.names(.)) %>%
    dplyr::left_join(run_order_pos)
  
  tempmod <- lm(Abundance ~ run_order, data = tempfitdat)
  predvals <- predict(tempmod, tempdat)
  
  # all(names(lpos_runorder_tc[i, 49:168]) == tempdat$SampleID)
  lpos_runorder_tc[i, 49:168] <- lpos_runorder_tc[i, 49:168] - predvals + median(tempdat %>%
                                                                                   dplyr::filter(!grepl("Blank", SampleID)) %>%
                                                                                   dplyr::filter(!grepl("PB", SampleID)) %>%
                                                                                   .$Abundance, na.rm = TRUE)
}

writexl::write_xlsx(x = list(`Negative` = lneg_runorder_tc ,
                             `Positive` = lpos_runorder_tc),
                    path = here("data", "processed", "ptrc_lipids_log2_sum_runordertc.xlsx"))
saveRDS(lneg_runorder_tc, here("data", "processed", "ptrc_lipids_neg_log2_sum_runordertc.rds"))
saveRDS(lpos_runorder_tc, here("data", "processed", "ptrc_lipids_pos_log2_sum_runordertc.rds"))

lneg_runorder_tc_sum <- lneg_runorder_tc
lpos_runorder_tc_sum <- lpos_runorder_tc

## Choice ------------------------------------------------------------------
# Processing -------------------------------------------------------------

run_order_neg$SampleID[!(run_order_neg$SampleID %in% names(lneg_choice_export_dat))]

names(lneg_choice_export_dat)[!(names(lneg_choice_export_dat) %in% make.names(run_order_neg$SampleID))]

# Extract QC samples to be used for loess curve fitting (i.e. every other QC sample)

lneg_qc_train <- lneg_choice_export_dat %>%
  dplyr::select(formatted_name:QC_CV, QC_Pool_Lip_001,
                QC_Pool_Lip_003, QC_Pool_Lip_01, QC_Pool_Lip_03,
                QC_Pool_Lip_05, QC_Pool_Lip_07)
lneg_qc_val <- lneg_choice_export_dat %>%
  dplyr::select(formatted_name:QC_CV, QC_Pool_Lip_002,
                QC_Pool_Lip_004, QC_Pool_Lip_02, QC_Pool_Lip_04,
                QC_Pool_Lip_06, QC_Pool_Lip_08)


run_order_pos$SampleID[!(run_order_pos$SampleID %in% names(lpos_choice_export_dat))]

names(lpos_choice_export_dat)[!(names(lpos_choice_export_dat) %in% make.names(run_order_pos$SampleID))]

# Extract QC samples to be used for loess curve fitting (i.e. every other QC sample)

lpos_qc_train <- lpos_choice_export_dat %>%
  dplyr::select(formatted_name:QC_CV, QC_Pool_Lip_001,
                QC_Pool_Lip_003, QC_Pool_Lip_01, QC_Pool_Lip_03,
                QC_Pool_Lip_05, QC_Pool_Lip_07)
lpos_qc_val <- lpos_choice_export_dat %>%
  dplyr::select(formatted_name:QC_CV, QC_Pool_Lip_002,
                QC_Pool_Lip_004, QC_Pool_Lip_02, QC_Pool_Lip_04,
                QC_Pool_Lip_06, QC_Pool_Lip_08)

# Apply correction/normalization (to qc and samples) ------------------------------------------

lneg_runorder_tc <- lneg_choice_export_dat
for(i in 1:nrow(lneg_runorder_tc)){
  
  tempdat <- lneg_runorder_tc[i,] %>%
    dplyr::select(dplyr::contains("Lip_"),
                  -dplyr::contains("PB"),
                  -dplyr::contains("Blank")) %>%
    t(.) %>%
    as.data.frame(.) %>%
    setNames("Abundance") %>%
    dplyr::mutate(SampleID = row.names(.)) %>%
    dplyr::left_join(run_order_neg)
  
  tempfitdat <- lneg_qc_train[i,] %>%
    dplyr::select(-dplyr::contains("Blk"),
                  -dplyr::contains("Blank"),
                  -(formatted_name:QC_CV)) %>%
    t(.) %>%
    as.data.frame(.) %>%
    setNames("Abundance") %>%
    dplyr::mutate(SampleID = row.names(.)) %>%
    dplyr::left_join(run_order_neg)
  
  tempmod <- lm(Abundance ~ run_order, data = tempfitdat)
  predvals <- predict(tempmod, tempdat)
  
  # all(names(lneg_runorder_tc[i, 47:166]) == tempdat$SampleID)
  lneg_runorder_tc[i, 47:166] <- lneg_runorder_tc[i, 47:166] - predvals + median(tempdat %>%
                                                                                   dplyr::filter(!grepl("Blank", SampleID)) %>%
                                                                                   dplyr::filter(!grepl("PB", SampleID)) %>%
                                                                                   .$Abundance, na.rm = TRUE)
}

lpos_runorder_tc <- lpos_choice_export_dat
for(i in 1:nrow(lpos_runorder_tc)){
  
  tempdat <- lpos_runorder_tc[i,] %>%
    dplyr::select(dplyr::contains("Lip_"),
                  -dplyr::contains("PB"),
                  -dplyr::contains("Blank")) %>%
    t(.) %>%
    as.data.frame(.) %>%
    setNames("Abundance") %>%
    dplyr::mutate(SampleID = row.names(.)) %>%
    dplyr::left_join(run_order_pos)
  
  tempfitdat <- lpos_qc_train[i,] %>%
    dplyr::select(-dplyr::contains("Blk"),
                  -dplyr::contains("Blank"),
                  -(formatted_name:QC_CV)) %>%
    t(.) %>%
    as.data.frame(.) %>%
    setNames("Abundance") %>%
    dplyr::mutate(SampleID = row.names(.)) %>%
    dplyr::left_join(run_order_pos)
  
  tempmod <- lm(Abundance ~ run_order, data = tempfitdat)
  predvals <- predict(tempmod, tempdat)
  
  # all(names(lpos_runorder_tc[i, 49:168]) == tempdat$SampleID)
  lpos_runorder_tc[i, 49:168] <- lpos_runorder_tc[i, 49:168] - predvals + median(tempdat %>%
                                                                                   dplyr::filter(!grepl("Blank", SampleID)) %>%
                                                                                   dplyr::filter(!grepl("PB", SampleID)) %>%
                                                                                   .$Abundance, na.rm = TRUE)
}

writexl::write_xlsx(x = list(`Negative` = lneg_runorder_tc ,
                             `Positive` = lpos_runorder_tc),
                    path = here("data", "processed", "ptrc_lipids_log2_choice_runordertc.xlsx"))
saveRDS(lneg_runorder_tc, here("data", "processed", "ptrc_lipids_neg_log2_choice_runordertc.rds"))
saveRDS(lpos_runorder_tc, here("data", "processed", "ptrc_lipids_pos_log2_choice_runordertc.rds"))

lneg_runorder_tc_choice <- lneg_runorder_tc
lpos_runorder_tc_choice <- lpos_runorder_tc




# -------------------------------------------------------------------------

# ### Prep Order Correction ----------------------------------------------------
# 
# liver_metadata <- readxl::read_excel("data/CPTAC_Liver_MPLex_Nov2024-BG_JK.xlsx", 
#                                      sheet = "PrepOrder") %>%
#   dplyr::select(`PNL ID`, `Mass (mg)`, `Sample Type`, `Overall MPLex Order`) %>%
#   setNames(c("SampleID", "mass", "samp_type", "prep_order")) %>%
#   dplyr::mutate(SampleID = paste0(SampleID, "_L"))
# 
# prep_order_neg <- liver_metadata %>%
#   dplyr::select(SampleID, prep_order)
# prep_order_pos <- liver_metadata %>%
#   dplyr::select(SampleID, prep_order)
# 
# ## Sum ---------------------------------------------------------------------
# # Processing -------------------------------------------------------------
# 
# prep_order_neg$SampleID[!(prep_order_neg$SampleID %in% names(lneg_runorder_tc_sum))]
# 
# names(lneg_runorder_tc_sum)[!(names(lneg_runorder_tc_sum) %in% make.names(prep_order_neg$SampleID))]
# 
# # Extract QC samples to be used for loess curve fitting (i.e. every other QC sample)
# # Need to find a way of extracting every other study sample efficiently.
# samp_names <- lneg_fdat %>%
#   dplyr::filter(Category == "Sample") %>%
#   .$SampleID
# oddseq <- seq(from = 1, to = length(samp_names), by = 2)
# evenseq <- seq(from = 2, to = length(samp_names)-1, by = 2)
# samp_names_train <- samp_names[oddseq]
# samp_names_test <- samp_names[evenseq]
# 
# lneg_samp_train <- lneg_runorder_tc_sum %>%
#   dplyr::select(formatted_name:QC_CV, dplyr::all_of(samp_names_train))
# lneg_samp_val <- lneg_runorder_tc_sum %>%
#   dplyr::select(formatted_name:QC_CV, dplyr::all_of(samp_names_test))
# 
# 
# prep_order_pos$SampleID[!(prep_order_pos$SampleID %in% names(lpos_runorder_tc_sum))]
# 
# names(lpos_runorder_tc_sum)[!(names(lpos_runorder_tc_sum) %in% make.names(prep_order_pos$SampleID))]
# 
# # Extract QC samples to be used for loess curve fitting (i.e. every other QC sample)
# 
# lpos_samp_train <- lpos_runorder_tc_sum %>%
#   dplyr::select(formatted_name:QC_CV, dplyr::all_of(samp_names_train))
# lpos_samp_val <- lpos_runorder_tc_sum %>%
#   dplyr::select(formatted_name:QC_CV, dplyr::all_of(samp_names_test))
# 
# # Apply correction/normalization (to qc and samples) ------------------------------------------
# 
# lneg_runpreporder_tc <- lneg_runorder_tc_sum
# for(i in 1:nrow(lneg_runpreporder_tc)){
#   
#   tempdat <- lneg_runpreporder_tc[i,] %>%
#     dplyr::select(dplyr::contains("CPTAC4_"),
#                   -dplyr::contains("Blk"),
#                   -dplyr::contains("Blank"),
#                   -dplyr::contains("QC_Pool")) %>%
#     t(.) %>%
#     as.data.frame(.) %>%
#     setNames("Abundance") %>%
#     dplyr::mutate(SampleID = row.names(.)) %>%
#     dplyr::left_join(prep_order_neg)
#   
#   tempfitdat <- lneg_samp_train[i,] %>%
#     dplyr::select(-dplyr::contains("Blk"),
#                   -dplyr::contains("Blank"),
#                   -dplyr::contains("QC_Pool"),
#                   -(formatted_name:QC_CV)) %>%
#     t(.) %>%
#     as.data.frame(.) %>%
#     setNames("Abundance") %>%
#     dplyr::mutate(SampleID = row.names(.)) %>%
#     dplyr::left_join(prep_order_neg)
#   
#   tempmod <- lm(Abundance ~ prep_order, data = tempfitdat)
#   predvals <- predict(tempmod, tempdat)
#   
#   # all(names(lneg_runpreporder_tc[i, 62:248]) == tempdat$SampleID)
#   lneg_runpreporder_tc[i, 62:248] <- lneg_runpreporder_tc[i, 62:248] - predvals + median(tempdat %>%
#                                                                                            dplyr::filter(!grepl("Blank", SampleID)) %>%
#                                                                                            dplyr::filter(!grepl("Blk", SampleID)) %>%
#                                                                                            dplyr::filter(!grepl("QC", SampleID)) %>%
#                                                                                            .$Abundance, na.rm = TRUE)
# }
# 
# lpos_runpreporder_tc <- lpos_runorder_tc_sum
# for(i in 1:nrow(lpos_runpreporder_tc)){
#   
#   tempdat <- lpos_runpreporder_tc[i,] %>%
#     dplyr::select(dplyr::contains("CPTAC4_"),
#                   -dplyr::contains("Blk"),
#                   -dplyr::contains("Blank"),
#                   -dplyr::contains("QC_Pool")) %>%
#     t(.) %>%
#     as.data.frame(.) %>%
#     setNames("Abundance") %>%
#     dplyr::mutate(SampleID = row.names(.)) %>%
#     dplyr::left_join(prep_order_pos)
#   
#   tempfitdat <- lpos_samp_train[i,] %>%
#     dplyr::select(-dplyr::contains("Blk"),
#                   -dplyr::contains("Blank"),
#                   -dplyr::contains("QC_Pool"),
#                   -(formatted_name:QC_CV)) %>%
#     t(.) %>%
#     as.data.frame(.) %>%
#     setNames("Abundance") %>%
#     dplyr::mutate(SampleID = row.names(.)) %>%
#     dplyr::left_join(prep_order_pos)
#   
#   tempmod <- lm(Abundance ~ prep_order, data = tempfitdat)
#   predvals <- predict(tempmod, tempdat)
#   
#   # all(names(lpos_runpreporder_tc[i, 62:248]) == tempdat$SampleID)
#   lpos_runpreporder_tc[i, 62:248] <- lpos_runpreporder_tc[i, 62:248] - predvals + median(tempdat %>%
#                                                                                            dplyr::filter(!grepl("Blank", SampleID)) %>%
#                                                                                            dplyr::filter(!grepl("Blk", SampleID)) %>%
#                                                                                            dplyr::filter(!grepl("QC", SampleID)) %>%
#                                                                                            .$Abundance, na.rm = TRUE)
# }
# 
# writexl::write_xlsx(x = list(`Negative` = lneg_runpreporder_tc ,
#                              `Positive` = lpos_runpreporder_tc),
#                     path = here("data", "processed", "ptrc_lipids_log2_sum_runprepordertc.xlsx"))
# saveRDS(lneg_runpreporder_tc, here("data", "processed", "ptrc_lipids_neg_log2_sum_runprepordertc.rds"))
# saveRDS(lpos_runpreporder_tc, here("data", "processed", "ptrc_lipids_pos_log2_sum_runprepordertc.rds"))
# 
# ## Choice ------------------------------------------------------------------
# # Processing -------------------------------------------------------------
# 
# prep_order_neg$SampleID[!(prep_order_neg$SampleID %in% names(lneg_choice_export_dat))]
# 
# names(lneg_choice_export_dat)[!(names(lneg_choice_export_dat) %in% make.names(prep_order_neg$SampleID))]
# 
# # Extract QC samples to be used for loess curve fitting (i.e. every other QC sample)
# 
# lneg_samp_train <- lneg_runorder_tc_choice %>%
#   dplyr::select(formatted_name:QC_CV, dplyr::all_of(samp_names_train))
# lneg_samp_val <- lneg_runorder_tc_choice %>%
#   dplyr::select(formatted_name:QC_CV, dplyr::all_of(samp_names_train))
# 
# 
# prep_order_pos$SampleID[!(prep_order_pos$SampleID %in% names(lpos_choice_export_dat))]
# 
# names(lpos_choice_export_dat)[!(names(lpos_choice_export_dat) %in% make.names(prep_order_pos$SampleID))]
# 
# # Extract QC samples to be used for loess curve fitting (i.e. every other QC sample)
# 
# lpos_samp_train <- lpos_runorder_tc_choice %>%
#   dplyr::select(formatted_name:QC_CV, dplyr::all_of(samp_names_train))
# lpos_samp_val <- lpos_runorder_tc_choice %>%
#   dplyr::select(formatted_name:QC_CV, dplyr::all_of(samp_names_train))
# 
# # Apply correction/normalization (to qc and samples) ------------------------------------------
# 
# lneg_runpreporder_tc <- lneg_runorder_tc_choice
# for(i in 1:nrow(lneg_runpreporder_tc)){
#   
#   tempdat <- lneg_runpreporder_tc[i,] %>%
#     dplyr::select(dplyr::contains("CPTAC4_"),
#                   -dplyr::contains("Blk"),
#                   -dplyr::contains("Blank"),
#                   -dplyr::contains("QC_Pool")) %>%
#     t(.) %>%
#     as.data.frame(.) %>%
#     setNames("Abundance") %>%
#     dplyr::mutate(SampleID = row.names(.)) %>%
#     dplyr::left_join(prep_order_neg)
#   
#   tempfitdat <- lneg_samp_train[i,] %>%
#     dplyr::select(-dplyr::contains("Blk"),
#                   -dplyr::contains("Blank"),
#                   -dplyr::contains("QC_Pool"),
#                   -(formatted_name:QC_CV)) %>%
#     t(.) %>%
#     as.data.frame(.) %>%
#     setNames("Abundance") %>%
#     dplyr::mutate(SampleID = row.names(.)) %>%
#     dplyr::left_join(prep_order_neg)
#   
#   tempmod <- lm(Abundance ~ prep_order, data = tempfitdat)
#   predvals <- predict(tempmod, tempdat)
#   
#   # all(names(lneg_runpreporder_tc[i, 62:248]) == tempdat$SampleID)
#   lneg_runpreporder_tc[i, 62:248] <- lneg_runpreporder_tc[i, 62:248] - predvals + median(tempdat %>%
#                                                                                            dplyr::filter(!grepl("Blank", SampleID)) %>%
#                                                                                            dplyr::filter(!grepl("Blk", SampleID)) %>%
#                                                                                            dplyr::filter(!grepl("QC", SampleID)) %>%
#                                                                                            .$Abundance, na.rm = TRUE)
# }
# 
# lpos_runpreporder_tc <- lpos_runorder_tc_choice
# for(i in 1:nrow(lpos_runpreporder_tc)){
#   
#   tempdat <- lpos_runpreporder_tc[i,] %>%
#     dplyr::select(dplyr::contains("CPTAC4_"),
#                   -dplyr::contains("Blk"),
#                   -dplyr::contains("Blank"),
#                   -dplyr::contains("QC_Pool")) %>%
#     t(.) %>%
#     as.data.frame(.) %>%
#     setNames("Abundance") %>%
#     dplyr::mutate(SampleID = row.names(.)) %>%
#     dplyr::left_join(prep_order_pos)
#   
#   tempfitdat <- lpos_samp_train[i,] %>%
#     dplyr::select(-dplyr::contains("Blk"),
#                   -dplyr::contains("Blank"),
#                   -dplyr::contains("QC_Pool"),
#                   -(formatted_name:QC_CV)) %>%
#     t(.) %>%
#     as.data.frame(.) %>%
#     setNames("Abundance") %>%
#     dplyr::mutate(SampleID = row.names(.)) %>%
#     dplyr::left_join(prep_order_pos)
#   
#   tempmod <- lm(Abundance ~ prep_order, data = tempfitdat)
#   predvals <- predict(tempmod, tempdat)
#   
#   # all(names(lneg_runpreporder_tc[i, 62:248]) == tempdat$SampleID)
#   lpos_runpreporder_tc[i, 62:248] <- lpos_runpreporder_tc[i, 62:248] - predvals + median(tempdat %>%
#                                                                                            dplyr::filter(!grepl("Blank", SampleID)) %>%
#                                                                                            dplyr::filter(!grepl("Blk", SampleID)) %>%
#                                                                                            dplyr::filter(!grepl("QC", SampleID)) %>%
#                                                                                            .$Abundance, na.rm = TRUE)
# }
# 
# writexl::write_xlsx(x = list(`Negative` = lneg_runpreporder_tc ,
#                              `Positive` = lpos_runpreporder_tc),
#                     path = here("data", "processed", "ptrc_lipids_log2_choice_runprepordertc.xlsx"))
# saveRDS(lneg_runpreporder_tc, here("data", "processed", "ptrc_lipids_neg_log2_choice_runprepordertc.rds"))
# saveRDS(lpos_runpreporder_tc, here("data", "processed", "ptrc_lipids_pos_log2_choice_runprepordertc.rds"))
# 
# 
# 
# 
# 
# 
# # -------------------------------------------------------------------------

# Save fdat
saveRDS(lneg_fdat, here("data", "processed", "ptrc_lipids_neg_fdat.rds"))
saveRDS(lpos_fdat, here("data", "processed", "ptrc_lipids_pos_fdat.rds"))

