## load data from diverse sources


library(tidyverse)
library(synapser)

synLogin()


##### pilot data metabolomics
mfile <- synGet('syn68710369')$path

hp <- readxl::read_xlsx(mfile,sheet = 'HILIC Positive') |>
  tidyr::pivot_longer(starts_with('PTRC'),names_to='Sample',values_to = 'value')|>
  mutate(Sample=stringr::str_replace(Sample,'QC_Pool_','QCPool'))|>  
  mutate(Sample=stringr::str_replace(Sample,'Blank_Pr_','BlankPr'))|>
  tidyr::separate(Sample,into=c('PTRC','exp','datatype','SampleID','runType','posneg'),sep='_')

hn <- readxl::read_xlsx(mfile,sheet = 'HILIC Negative')|>
  tidyr::pivot_longer(starts_with('PTRC'),names_to='Sample',values_to = 'value')|>
  mutate(Sample=stringr::str_replace(Sample,'QC_Pool_','QCPool'))|>  
  mutate(Sample=stringr::str_replace(Sample,'Blank_Pr_','BlankPr'))|>
  tidyr::separate(Sample,into=c('PTRC','exp','datatype','SampleID','runType','posneg'),sep='_')

rp <- readxl::read_xlsx(mfile,sheet = 'RP Positive')|>
  tidyr::pivot_longer(starts_with('PTRC'),names_to='Sample',values_to = 'value')|>
  mutate(Sample=stringr::str_replace(Sample,'QC_Pool_','QCPool'))|>  
  mutate(Sample=stringr::str_replace(Sample,'Blank_Pr_','BlankPr'))|>
  tidyr::separate(Sample,into=c('PTRC','exp','datatype','SampleID','runType','posneg'),sep='_')

rn <- readxl::read_xlsx(mfile,sheet = 'RP Negative')|>
  tidyr::pivot_longer(starts_with('PTRC'),names_to='Sample',values_to = 'value')|>
  mutate(Sample=stringr::str_replace(Sample,'QC_Pool_','QCPool'))|>  
  mutate(Sample=stringr::str_replace(Sample,'Blank_Pr_','BlankPr'))|>
  tidyr::separate(Sample,into=c('PTRC','exp','datatype','SampleID','runType','posneg'),sep='_')

full_metab <- rbind(hp,hn,rp,rn)