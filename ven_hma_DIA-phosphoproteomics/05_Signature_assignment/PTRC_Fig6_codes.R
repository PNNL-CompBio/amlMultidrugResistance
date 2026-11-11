## Function to compute signal value using two sets of markers. The signal will be mean(no relapse_markers) - mean(refractory_markers)
## NOTE: The markers are intersected with the set of COMPLETE (NO NA) features from the new data.
predict_helper_ <- function(new_msnset, case_markers, control_markers){
  ## intersect with complete features from exp28
  complete_features = featureNames(new_msnset)[(rowSums(is.na(exprs(new_msnset))) == 0)]
  case_markers = intersect(case_markers, complete_features)
  control_markers = intersect(control_markers, complete_features)
  case_mat = t(exprs(new_msnset)[case_markers, , drop = FALSE]) %>% as.data.frame()
  control_mat = t(exprs(new_msnset)[control_markers, , drop = FALSE]) %>% as.data.frame()
  if (ncol(control_mat) > 0){
    control_mat$control_avg = rowMeans(control_mat, na.rm = T)
    control_mat <- control_mat %>% select(control_avg)
  }
  if (ncol(case_mat) > 0){
    case_mat$case_avg = rowMeans(case_mat, na.rm = T)
    case_mat <- case_mat %>% select(case_avg)
  }
  logistic_df <- data.frame(sample = sampleNames(new_msnset)) %>% 
    cbind(case_mat, control_mat) %>%
    mutate(signal = case_avg - control_avg)
  return(logistic_df)
}



diffexp = filter(tests2, contrast=='Norelapse-Refractory')

refractory_markers = diffexp %>% filter(logFC < 0, P.Value < 0.05) %>% pull(feature)
longresponse_markers = diffexp %>% filter(logFC > 0, P.Value  < 0.05) %>% pull(feature)
xx = predict_helper_(dfSet, refractory_markers, longresponse_markers) %>%
  left_join(., pData(dfSet), by = c("sample" = "SampleID")) 

xx$subcohort <- factor(xx$subcohort, levels=c("Norelapse", "Relapse", "Refractory" ))
ggplot(xx, aes(x = subcohort, y = signal, fill = subcohort)) + 
  geom_boxplot()+
  scale_x_discrete(labels=c("Norelapse"="Long response", "Relapse"="Short response", "Refractory"="Refractory"))+
  theme(panel.background = element_blank(),
        axis.line=element_line(color="black"),
        axis.text.x = element_text(size=30, angle=0, face="bold", vjust=0.5, color="black"),
        axis.ticks.x =  element_line(colour = "black", size=1),
        axis.line.x = element_line(colour = "black", size=2),
        axis.line.y = element_line(colour = "black", size=2),
        axis.title.x = element_blank(),
        axis.text.y = element_text(face="bold", color="black", size=30),
        axis.title.y = element_text(face= "bold", size= 35),
        legend.position ="right",
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size=25),
        legend.title = element_text(size=25, face = "bold"))


## apply prediction trained from EXP28 to BeatAML phospho data
BeatAML <- read.xlsx("BeatAML_phosphosite.xlsx", check.names=FALSE)
BeatAMLmeta <- read.xlsx("BeatAML_metadata.xlsx", check.names=FALSE)

BeatAML <- BeatAML%>% mutate(across('01':'99', as.numeric))

## the site column containing multiple sites per row. e.g. AAK1-S618sT620t. Clean them up to two duplicate rows as AAK1-S618, AAK1-T620
BeatAML_expanded <- BeatAML %>%
  mutate(
    gene = str_extract(Site, "^[^-]+"),              
    sites_raw = str_extract(Site, "(?<=-).*"),      
    site_list = str_extract_all(sites_raw, "[A-Z]\\d+")  
  ) %>%
  unnest(site_list) %>%                                    
  mutate(gene_site = paste0(gene, "-", site_list)) %>%    
  select(-gene, -sites_raw, -site_list) %>%              
  relocate(gene_site) 

write.xlsx(BeatAML_expanded, "BeatAML_duplicateSite.xlsx",rowNames=FALSE)

BeatAML <- read.xlsx("BeatAML_duplicateSite.xlsx", check.names=FALSE)

BeatAML <- BeatAML%>% select(,-2) %>% group_by(gene_site) %>%
  summarise(gene_site=dplyr::first(gene_site),
            across(everything(), sum, na.rm=TRUE))

## check how many phosphosites overlap between EXP28 and BeatAML
comSite <- merge(BeatAML, df, by.x="gene_site", by.y="SITE2")

## build a MSnSet using BeatAML data
exprs <- BeatAML %>% arrange(., gene_site) %>% tibble::column_to_rownames(var="gene_site") %>% as.matrix()
fData <- data.frame(rownames(exprs)) %>% tibble::column_to_rownames(var="rownames.exprs.")
exprs <- exprs[ , paste(BeatAMLmeta$SampleID.abbrev, sep = "")]
pData <- BeatAMLmeta %>% mutate(SampleID = SampleID.abbrev) %>% tibble::column_to_rownames(var="SampleID.abbrev")

BeatAML_Set <- MSnSet(exprs = exprs,
                      pData = pData, fData = fData)

beataml_210_predictions = predict_helper_(BeatAML_Set, refractory_markers, longresponse_markers) %>%
  left_join(., pData(BeatAML_Set), by = c("sample" = "SampleID"))

##beataml_210_predictions %>%
##ggplot(aes(x = auc, y=signal)) + geom_density(alpha = 0.2)

## filter out BeatAML samples that are initial diagnosis or initial diagnosis;post-chemo
ll = beataml_210_predictions %>% filter(Specimen.access.group.concatenated %in% c("Initial Acute Leukemia Diagnosis",
                                                                                  "Initial Acute Leukemia Diagnosis; Post-Chemotherapy")) %>%
  filter(specimen.location == "OHSU")
ll = beataml_210_predictions %>% filter(specimen.location == "OHSU")
ll <- ll %>% mutate(auc = as.numeric(auc))

a = cor.test(ll$auc, ll$signal)
a
line_slope = coef(lm(ll$signal ~ ll$auc))[[2]]
line_intercept = coef(lm(ll$signal ~ ll$auc))[[1]]
plot_title = paste0("clinical signature vs BeatAML Ven+Aza drug AUC \n", "correlation = ", formatC(a$estimate, 2),
                    ", p_value = ", formatC(a$p.value, 2))
ggplot(ll, aes(x = auc, y = signal, label = sample)) + geom_point(size = 4, color="blue") +
  ggtitle(plot_title) + ylab("signal") + xlab("Ven+Aza drug AUC") + 
  geom_abline(intercept = line_intercept, slope = line_slope, size = 2, alpha = 0.3) +
  theme_bw() + theme(text = element_text(size = 25), axis.text=element_text(color="black"), plot.margin=margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"))