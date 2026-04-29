library(dplyr)
library(MSnSet.utils)
library(openxlsx)

## load log2 transformed, median centered phosphosite level data 
df <- read.xlsx("PTRC_EXP28_InSilico_Cleaned_PhosphositeData.xlsx", check.names=FALSE)
meta <- read.xlsx("PTRC_metadata_Exp28_removesamples.xlsx",check.names=FALSE)


## build MSnSet
exprs <- df %>% arrange(., SITE) %>% tibble::column_to_rownames(var="SITE") %>% as.matrix()
fData <- data.frame(rownames(exprs)) %>% tibble::column_to_rownames(var="rownames.exprs.")
exprs <- exprs[ , paste(meta$Sample, sep = "")]
pData <- meta %>% mutate(SampleID = Sample) %>% tibble::column_to_rownames(var="Sample")

dfSet <- MSnSet(exprs = exprs,
                pData = pData, fData = fData)

# set up contrasts to compare 1) responders VS nonresponders 2) no-relapse VS refractory, relapse VS refractory, relapse VS no-relapse
contrasts1 <- c("ResponseGroupResponders-ResponseGroupNonResponders")
contrasts2 <- c("subcohortNorelapse-subcohortRefractory", "subcohortRelapse-subcohortRefractory", "subcohortRelapse-subcohortNorelapse")

tests1 <- limma_contrasts(eset = dfSet, model.str = "~ 0 + ResponseGroup", coef.str = "ResponseGroup",
                          contrasts = contrasts1, trend = TRUE, robust = TRUE, plot = FALSE)

tests2 <- limma_contrasts(eset = dfSet, model.str = "~ 0 + subcohort", coef.str = "subcohort",
                          contrasts = contrasts2, trend = TRUE, robust = TRUE, plot = FALSE)

write.xlsx(tests1, "PTRC_EXP28_Responder VS NonResponder.xlsx", rowNames=FALSE)
write.xlsx(tests2, "PTRC_EXP28_Refractory VS relapse and no-relapse.xlsx", rowNames=FALSE)
