library(dplyr)
library(stringr)
library(openxlsx)

## load phosphosite level data 
df <- read.xlsx("PTRC_EXP28_InSilico_DiaNN_phosphosites_90_removesamples.xlsx", check.names=FALSE)

## make sure columns containing intensity values are in numeric
df <- df%>% mutate(across('PTRC_Exp28_Phos_01':'PTRC_Exp28_Phos_41', as.numeric))

## different uniprot protein accession IDs could be mapped to the same gene name. Create a new SITE2 column that is formatted as GeneName-Residue#, e.g.TADA2A-S6
df <- df %>% mutate(SITE = str_c(Gene.Names, "-", Residue, Site))

## subset site and sample columns
df <- df[,9:45]
## sum rows with same SITE ID
df <- df %>% group_by(SITE) %>%
  summarise(SITE=dplyr::first(SITE),
            across(everything(), sum, na.rm=TRUE))

## replace 0 with NA, followed by log2 transformation
df[df == 0] <- NA
df[,2:37] <- log(df[, 2:37], 2)
# check data distribution prior to median centering
boxplot(df[,2:37], cex.axis=1, las=2) 

# median centering
Zero_Center_Norm <- function(df) {
  med_norm <- function (df) 
  {
    norm.coeff <- apply(df, 2, median, na.rm = TRUE)# collect median of each sample from specified dataframe
    df1 <- sweep(df, 2, norm.coeff, "-") #subtract the median from each respective column in dataframe
    avg_of_median <- mean(norm.coeff) # calculate average of medians of each sample in group
    df1 <- df1 + avg_of_median #add average of averages back to each subtracted sample value
    return(df1)
  }
  df <- med_norm(df)
  return(df)
}

df[,2:37] <- Zero_Center_Norm(df[,2:37])

# check data distribution after median centering
boxplot(df[,2:37], cex.axis=1, las=2) 
# remove sites with NAs across all samples.
df <- df[rowSums(!is.na(df[ , 2:37])) > 0, ] 

write.xlsx(df, "PTRC_EXP28_InSilico_Cleaned_PhosphositeData.xlsx", rowNames=FALSE)
