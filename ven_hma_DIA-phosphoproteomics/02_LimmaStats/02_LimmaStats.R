library(dplyr)
library(MSnSet.utils)
library(openxlsx)
library(synapser)

## synapse login with .Renviron
synLogin()

## define folder ID and the target files
folder_id <-"syn68733653"

## load log2 transformed, median centered phosphosite level data 
file1_id <- synFindEntityId("PTRC_EXP28_InSilico_Cleaned_PhosphositeData.xlsx", parent = folder_id)
file1_entity <- synGet(file1_id)
df <- read.xlsx(file1_entity$path)

file2_id <- synFindEntityId("PTRC_metadata_Exp28_removesamples.xlsx", parent = folder_id)
file2_entity <- synGet(file2_id)
meta <- read.xlsx(file2_entity$path)

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


##now we can do volcano plots for figure 2
##Fig 5A
nr_vs_ref <- subset(tests2, contrast=='Norelapse-Refractory') |>
  drop_na()|>
  mutate(color = ifelse(logFC>0,'red','blue'))|>
  mutate(mLogp=-1*log10(P.Value))

nr_vs_ref$color[which(nr_vs_ref$P.Value >0.05)] <-'grey'

counts <- nr_vs_ref|>group_by(color)|> summarize(n=n_distinct(feature))
up = counts$n[counts$color=='red']
down = counts$n[counts$color=='blue']

ggplot(nr_vs_ref,aes(x=logFC,y=mLogp,col=color))+geom_point() + geom_hline(yintercept = 1.3,linetype='dashed')+
  scale_color_manual(values=list(red='darkred',blue='darkblue',grey='lightgrey'))+  
  theme_minimal(base_size = 12) +
  annotate('text',label=paste('n =',down),x = mean(nr_vs_ref$logFC[nr_vs_ref$color=='blue']), y = 6,col='darkblue')+
  annotate('text',label=paste('n =',up),x = mean(nr_vs_ref$logFC[nr_vs_ref$color=='red']), y = 6,col='darkred')+
  labs(x = "log2 fold change", y = '-log10(p value)', title = "")

ggsave('fig5a_volcano_nr_vs_ref.pdf')

##now do the same for relapse/refrac
##Fig 5B
r_vs_ref <- subset(tests2, contrast=='Relapse-Refractory') |>
  drop_na()|>
  mutate(color = ifelse(logFC>0,'red','blue'))|>
  mutate(mLogp=-1*log10(P.Value))

r_vs_ref$color[which(r_vs_ref$P.Value >0.05)] <-'grey'

counts <- r_vs_ref|>group_by(color)|> summarize(n=n_distinct(feature))
up = counts$n[counts$color=='red']
down = counts$n[counts$color=='blue']

ggplot(r_vs_ref,aes(x=logFC,y=mLogp,col=color))+geom_point() + geom_hline(yintercept = 1.3,linetype='dashed')+
  scale_color_manual(values=list(red='darkred',blue='darkblue',grey='lightgrey'))+  
  theme_minimal(base_size = 12) +
  annotate('text',label=paste('n =',down),x = mean(r_vs_ref$logFC[r_vs_ref$color=='blue']), y = 6,col='darkblue')+
  annotate('text',label=paste('n =',up),x = mean(r_vs_ref$logFC[r_vs_ref$color=='red']), y = 6,col='darkred')+
  labs(x = "log2 fold change", y = '-log10(p value)', title = "")

ggsave('fig5b_volcano_r_vs_ref.pdf')

write.xlsx(tests1, "PTRC_EXP28_Responder VS NonResponder.xlsx", rowNames=FALSE)
write.xlsx(tests2, "PTRC_EXP28_Refractory VS relapse and no-relapse.xlsx", rowNames=FALSE)

combined <- rbind(tests1,tests2) |>
  arrange(P.Value)

write.xlsx(combined,'Table_S3_LIMMA_results.xlsx',rowNames=F)
synStore(File('Table_S3_LIMMA_results.xlsx', parent='syn76730214'))