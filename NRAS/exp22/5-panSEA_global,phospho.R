# Differential expression & enrichment analyses: global & phospho
# PTRC2: exp 22
# Author: Belinda B. Garana
# Created: 2023-12-06
# Last edit: 2023-12-14

library(readxl); library(panSEA); library(synapser)
library(stringr); library(tidyr)

#### 1. Import metadata & crosstabs ####
setwd(
"~/OneDrive - PNNL/Documents/GitHub/Exp22_CRISPR-Cas9_NRAS_KO_cell_lines/proteomics/data/"
)
meta.df <- readxl::read_excel("Expt 22 Design.xlsx", 
                              sheet = 2)
global.df <- read.table(
  "global_data/ptrc_ex22_crosstab_global_gene_corrected.txt", 
  sep = "\t")
phospho.df <- read.table(
  "phospho_data/ptrc_ex22_crosstab_phospho_SiteID_corrected.txt", 
  sep = "\t")

# add column for feature names and later make it the first column
global.df$Gene <- rownames(global.df)
phospho.df$SUB_SITE <- rownames(phospho.df)

# add more descriptive treatment column
meta.df$Treatment <- meta.df$`Sample Name`
meta.df[grep("NRAS", meta.df$`Sample Name`), ]$Treatment <- substr(
  meta.df[grep("NRAS", meta.df$`Sample Name`), ]$'Sample Name', 1, 
  nchar(meta.df[grep("NRAS", meta.df$`Sample Name`), ]$'Sample Name') - 3) # remove last 3 chars

#### 2. Import BeatAML data ####
# import drug MOA annotations
moa.BeatAML <- utils::read.csv(
  "~/OneDrive - PNNL/Documents/PTRC2/BeatAML_single_drug_moa.csv",
  stringsAsFactors = FALSE, fileEncoding = "latin1")

# login to Synapse
synLogin()

# load data from Synapse
BeatAML.path <- "BeatAML_DMEA_inputs"
BeatAML_synapse_id <- list("drug_response.csv" = "syn51674470", 
                        "Ex10_metadata.txt" = "syn25807733",
                        "ptrc_ex10_crosstab_global_gene_corrected.txt" = "syn25714248",
                        "ptrc_ex10_crosstab_phospho_siteID_corrected(1).txt" = "syn25714921")
lapply(BeatAML_synapse_id, synapser::synGet, downloadLocation = BeatAML.path)

setwd(BeatAML.path)
drug.BeatAML <- read.csv(names(BeatAML_synapse_id)[1])
meta.BeatAML <- read.table(names(BeatAML_synapse_id)[2], sep = "\t", 
                           header = TRUE)
global.BeatAML <- read.table(names(BeatAML_synapse_id)[3], sep = "\t", 
                             header = TRUE)
phospho.BeatAML <- read.table(names(BeatAML_synapse_id)[4], sep = "\t", 
                              header = TRUE)

#### Step 4. Format BeatAML data for DMEA ####
sample.names <- "Barcode.ID"

## format drug sensitivity data frame
# format drug.BeatAML wide (samples in first column, drug names for rest of columns)
drug.BeatAML <- reshape2::dcast(drug.BeatAML, sample_id ~ inhibitor, 
                                value.var = "auc", fill = NA)

# remove drugs without moa annotations and drug combos
valid.drugs <- 
  names(drug.BeatAML)[names(drug.BeatAML) %in% 
                        moa.BeatAML[!is.na(moa.BeatAML),]$Drug] # 167 drugs
drug.BeatAML <- drug.BeatAML[ , c("sample_id", valid.drugs)] # 167 drugs
moa.BeatAML <- 
  moa.BeatAML[moa.BeatAML$Drug %in% names(drug.BeatAML)[2:ncol(drug.BeatAML)], ]

# change sample column name to match expression data
names(drug.BeatAML)[1] <- sample.names

## format global proteomics data frame
# change global.BeatAML column names from SampleID.abbrev to 
# Barcode.ID to match drug.BeatAML
global.ids <- names(global.BeatAML)

# remove X and any 0's from start of each column name and then
# replace SampleID.abbrev with Barcode.ID to match drug.BeatAML
for(i in seq_len(length(global.ids))){
  global.ids[i] <- substr(global.ids[i], 2, nchar(global.ids[i]))
  
  if(substring(global.ids[i], 1, 1) == 0){
    global.ids[i] <- substr(global.ids[i], 2, nchar(global.ids[i]))
  }
  
  if(global.ids[i] %in% meta.BeatAML$SampleID.abbrev){
    global.ids[i] <- meta.BeatAML[meta.BeatAML$SampleID.abbrev == global.ids[i], ]$Barcode.ID
  }
}

# replace global.BeatAML column names 
names(global.BeatAML) <- global.ids

# transpose global.BeatAML so that first column is Barcode.ID and 
# rest of columns are gene symbols
global.BeatAML <- as.data.frame(t(global.BeatAML))

# make first column Barcode.ID
global.BeatAML[, sample.names] <- rownames(global.BeatAML)
global.BeatAML <- 
  global.BeatAML[ , c(sample.names, 
                      names(global.BeatAML[ , 1:(ncol(global.BeatAML)-1)]))]

## format phospho-proteomics data frame
# change global.BeatAML column names from SampleID.abbrev to Barcode.ID to match drug.BeatAML
phospho.ids <- names(phospho.BeatAML)

# remove X and any 0's from start of each column name and then
# replace SampleID.abbrev with Barcode.ID to match drug.BeatAML
for(i in seq_len(length(phospho.ids))){
  phospho.ids[i] <- substr(phospho.ids[i], 2, nchar(phospho.ids[i]))

  if(substring(phospho.ids[i], 1, 1) == 0){
    phospho.ids[i] <- substr(phospho.ids[i], 2, nchar(phospho.ids[i]))
  }

  if(phospho.ids[i] %in% meta.BeatAML$SampleID.abbrev){
    phospho.ids[i] <- meta.BeatAML[
      meta.BeatAML$SampleID.abbrev == phospho.ids[i], ]$Barcode.ID
  }
}

# replace phospho.BeatAML column names
names(phospho.BeatAML) <- phospho.ids

# transpose phospho.BeatAML so that first column is Barcode.ID and rest of columns are gene symbols
phospho.BeatAML <- as.data.frame(t(phospho.BeatAML))

# make first column Barcode.ID
phospho.BeatAML[, sample.names] <- rownames(phospho.BeatAML)
phospho.BeatAML <- phospho.BeatAML[ , c(sample.names, names(phospho.BeatAML[ , 1:(ncol(phospho.BeatAML)-1)]))]

#### 2. Run panSEA across contrasts for each exp, omics type ####
# synapse IDs must match order of omics list
# make new folders in Analysis folder on Synapse: https://www.synapse.org/#!Synapse:syn39679397
synapse_id_map <- c("syn53212346" = "global_data/",
                    "syn53212349" = "phospho_data/")

## prepare other input parameters
# identify contrasts
treatments <- na.omit(unique(meta.df$Treatment))

# set short-hand names for treatments for better plots
# IMPORTANT: these must have the same order as treatments above
vols <- c("7L", "8L", "12L")
#treatment.names <- c("Gilt LR", "NT", "NRAS")
treatment.names <- c("NRAS KO", "Ctl")

# define types based on short-hand contrasts
types <- c(
  # compare to control
  "NRAS KO vs. Ctl"
)

types <- paste(vols, types, sep = ": ")

# get contrasts with full-length names to extract data
# contrasts <- c()
# for (i in 1:length(types)) {
#   # get treatment names (short-hand)
#   contrast.treatment.names <- stringr::str_split(types[i], " vs. ")[[1]]
#   
#   # get full treatment names
#   index1 <- treatments[treatment.names == contrast.treatment.names[1]]
#   index2 <- treatments[treatment.names == contrast.treatment.names[2]]
#   
#   # define contrast with full treatment name
#   contrasts[i] <- paste(index1, "vs.", index2)
# }

# in this case, the contrast is always NRAS vs. Ctl for each volume type

# identify column names containing data for each treatment
# add X in front of tube names to match column names 
# after import (e.g., X5 for tube 5)
#sample.groups <- list()
# for (i in 1:length(treatments)) {
#   sample.groups[[treatments[i]]] <- 
#    paste0("X", na.omit(meta.df[meta.df$Treatment == treatments[i], ]$Tube))
# }

#sample.groups[["NRAS KO"]] <- paste0("X", na.omit(meta.df[grep("NRAS", meta.df$Treatment), ]$Tube))
#sample.groups[["Ctl"]] <- paste0("X", na.omit(meta.df[grep("guide", meta.df$Treatment), ]$Tube))

#sample.groups <- list(1:3, 4)
#all.files2 <- c(DEG.files, GSEA.files, DMEA.files)
#subsets <- c("Differential expression", "GSEA", "DMEA")
#setwd("~/OneDrive - PNNL/Documents/GitHub/Exp22_CRISPR-Cas9_NRAS_KO_cell_lines/proteomics/data/")
base.path <- "~/OneDrive - PNNL/Documents/GitHub/Exp22_CRISPR-Cas9_NRAS_KO_cell_lines/proteomics/data/"
omics <- c("global", "phospho")
for (k in 1:length(omics)) {
  setwd(paste0(base.path, synapse_id_map[k]))
  
  ## prepare set annotations
  # generate gmt.features beforehand to save time
  if (omics[k] == "global") {
    msigdb.info <- msigdbr::msigdbr("Homo sapiens", "C2", "CP:KEGG")
    
    # extract necessary info into data frame
    msigdb.info <- as.data.frame(msigdb.info[, c(
      "gene_symbol",
      "gs_name",
      "gs_description"
    )])
    
    gmt <- DMEA::as_gmt(
      msigdb.info, "gene_symbol", "gs_name", min.per.set = 6,
      descriptions = "gs_description"
    ) 
  } else if (omics[k] == "phospho") {
    # # only create gmt the first time
    # use PNNL phospho feature IDs instead of ksdb
    SUB_SITE <- phospho.df$SUB_SITE
    phospho.ref <- data.frame(SUB_SITE)
    phospho.ref <- phospho.ref %>% tidyr::extract(SUB_SITE, "KINASE",
                                                  remove = FALSE)
    SUB_SITE <- NULL
    gmt <- DMEA::as_gmt(phospho.ref, "SUB_SITE", "KINASE", min.per.set = 6)

    # save gmt for future analyses
    saveRDS(gmt, "gmt_PNNL_kinase-substrate_PTRC2_exp22.rds")
    
    gmt <- readRDS("gmt_PNNL_kinase-substrate_PTRC2_exp22.rds")
  }
  
  # run panSEA for each omics type across all contrasts
  # CAUTION: this only works because the # of samples for each treatment type 
  # is equal; otherwise would have to run panSEA for each contrast separately 
  # and perhaps set the group.samples input parameter for panSEA
  
  # assemble inputs
  data.list <- list()
  expr.BeatAML <- list()
  gmt.features <- list()
  #sample.groups <- list()
  for (i in 1:length(vols)) {
    # identify treatments in this contrast
    vol.treatments <- treatments[grep(vols[i], treatments)]
    vol.treatment1 <- vol.treatments[grep("NRAS", vol.treatments)]
    vol.treatment2 <- vol.treatments[grep("guide", vol.treatments)]
    vol.treatments12 <- c(vol.treatment1, vol.treatment2)
    
    #vol.tubes <- paste0("X", meta.df[meta.df$Treatment %in% vol.treatments12, ]$Tube)
    vol.tubes1 <- paste0("X", meta.df[meta.df$Treatment %in% vol.treatment1, ]$Tube)
    vol.tubes2 <- paste0("X", meta.df[meta.df$Treatment %in% vol.treatment2, ]$Tube)
    vol.tubes12 <- c(vol.tubes1, vol.tubes2)
    #sample.groups[[types[i]]] <- list(vol.tubes1, vol.tubes2)
    
    # assemble each input data set for each contrast
    if (omics[k] == "global"){
      data.list[[types[i]]] <- 
        global.df[ , c("Gene", vol.tubes12)]
      
      colnames(data.list[[types[i]]]) <- c("Gene", "NRAS KO 1", "NRAS KO 2",
                                           "NRAS KO 3", "NT Guide")
      
      expr.BeatAML[[types[i]]] <- global.BeatAML 
      
      feature.names <- rep("Gene", length(types))
    } else if (omics[k] == "phospho") {
      data.list[[types[i]]] <- 
        phospho.df[ , c("SUB_SITE", vol.tubes12)]
      
      colnames(data.list[[types[i]]]) <- c("SUB_SITE", "NRAS KO 1", "NRAS KO 2",
                                           "NRAS KO 3", "NT Guide")
      
      expr.BeatAML[[types[i]]] <- phospho.BeatAML
      
      feature.names <- rep("SUB_SITE", length(types))
    }
    
    gmt.features[[types[i]]] <- gmt
  }
  
  # sample.groups <- list(c("NRAS KO 1", "NRAS KO 2", "NRAS KO 3"),
  #                       "NT Guide")
  sample.groups <- list(2:4, 5)
  
  # run panSEA by querying BeatAML data
  #Sys.setenv('R_MAX_VSIZE'=32000000000)
  #library(usethis) 
  #usethis::edit_r_environ()
  min.per.set <- 6
  panSEA.BeatAML <- panSEA::panSEA(data.list, types, 
                                   feature.names = feature.names,
                                   group.samples = sample.groups,
                                   gmt.features = gmt.features,
                                   gmt.drugs = DMEA::as_gmt(moa.BeatAML, 
                                                            min.per.set = 
                                                              min.per.set),
                                   drug.sensitivity = drug.BeatAML,
                                   expression = expr.BeatAML, 
                                   min.per.set = min.per.set)
  
  ## save results & upload to Synapse
  # store all results locally
  dir.create("analysis")
  setwd("analysis")
  
  # save phospho DMEA results because accidentally forgot
  #DEGs <- read.csv("Differential expression/Differential_expression_results.csv")
  #DEG.list <- list()
  # for (i in 1:length(types)) {
  #   temp.degs <- DEGs[DEGs$type == types[i], c("feature", "Log2FC")]
  #   colnames(temp.degs)[1] <- "SUB_SITE"
  #   DEG.list[[types[i]]] <- temp.degs
  # }
  # panSEA.BeatAML <- list()
  # panSEA.BeatAML[["mDMEA.results"]] <- panSEA::mDMEA(drug.sensitivity = drug.BeatAML, 
  #                                 gmt = DMEA::as_gmt(moa.BeatAML, 
  #                                                    min.per.set = min.per.set),
  #                                 expression = list(phospho.BeatAML, phospho.BeatAML, phospho.BeatAML),
  #                                 weights = DEG.list, types,
  #                                 feature.names = rep("SUB_SITE", length(types)),
  #                                 min.per.set = min.per.set)
  # # compile inputs & outputs for network graph
  # inputs <- list()
  # for (i in 1:length(types)) {
  #   inputs[[types[i]]] <- panSEA.BeatAML[["mDMEA.results"]]$all.results[[types[i]]]$corr.result
  # }
  # 
  # outputs <- list()
  # for (i in 1:length(types)) {
  #   outputs[[types[i]]] <- panSEA.BeatAML[["mDMEA.results"]]$all.results[[types[i]]]$result
  # }
  # 
  # panSEA.BeatAML[["mDMEA.network"]] <- panSEA::netSEA(
  #   inputs, outputs, rep("Drug", length(inputs)),
  #   rank.var = rep("Pearson.est", length(types)), p = 0.05, FDR = 0.25, 
  #   n.network.sets = 2*length(types), scale = 5
  # )
  # 
  # 
  # set file names
  DEG.files <- list("Differential_expression_results.csv" = 
                      panSEA.BeatAML$mDEG.results$compiled.results$results,
                    "Differential_expression_mean_results.csv" =
                      panSEA.BeatAML$mDEG.results$compiled.results$mean.results,
                    "Differential_expression_correlation_matrix.pdf" =
                      panSEA.BeatAML$mDEG.results$compiled.results$corr.matrix,
                    "Differential_expression_venn_diagram.pdf" =
                      panSEA.BeatAML$mDEG.results$compiled.results$venn.diagram,
                    "Differential_expression_dot_plot.pdf" =
                      panSEA.BeatAML$mDEG.results$compiled.results$dot.plot)
  DMEA.files <- list("DMEA_results.csv" =
                       panSEA.BeatAML$mDMEA.results$compiled.results$results,
                     "DMEA_mean_results.csv" =
                       panSEA.BeatAML$mDMEA.results$compiled.results$mean.results,
                     "DMEA_correlation_matrix.pdf" =
                       panSEA.BeatAML$mDMEA.results$compiled.results$corr.matrix,
                     "DMEA_venn_diagram.pdf" =
                       panSEA.BeatAML$mDMEA.results$compiled.results$venn.diagram,
                     "DMEA_dot_plot.pdf" =
                       panSEA.BeatAML$mDMEA.results$compiled.results$dot.plot,
                     "DMEA_interactive_network.graph.html" =
                       panSEA.BeatAML$mDMEA.network$interactive)
  if (omics[k] == "phospho") {
    KSEA.files <- list("KSEA_results.csv" =
                         panSEA.BeatAML$mGSEA.results$compiled.results$results,
                       "KSEA_mean_results.csv" =
                         panSEA.BeatAML$mGSEA.results$compiled.results$mean.results,
                       "KSEA_correlation_matrix.pdf" =
                         panSEA.BeatAML$mGSEA.results$compiled.results$corr.matrix,
                       "KSEA_venn_diagram.pdf" =
                         panSEA.BeatAML$mGSEA.results$compiled.results$venn.diagram,
                       "KSEA_dot_plot.pdf" =
                         panSEA.BeatAML$mGSEA.results$compiled.results$dot.plot,
                       "KSEA_interactive_network.graph.html" =
                         panSEA.BeatAML$mGSEA.network$interactive)
    all.files <- list('Differential expression' = DEG.files,
                      'KSEA' = KSEA.files,
                      'DMEA' = DMEA.files)
  } else {
    GSEA.files <- list("GSEA_results.csv" =
                         panSEA.BeatAML$mGSEA.results$compiled.results$results,
                       "GSEA_mean_results.csv" =
                         panSEA.BeatAML$mGSEA.results$compiled.results$mean.results,
                       "GSEA_correlation_matrix.pdf" =
                         panSEA.BeatAML$mGSEA.results$compiled.results$corr.matrix,
                       "GSEA_venn_diagram.pdf" =
                         panSEA.BeatAML$mGSEA.results$compiled.results$venn.diagram,
                       "GSEA_dot_plot.pdf" =
                         panSEA.BeatAML$mGSEA.results$compiled.results$dot.plot,
                       "GSEA_interactive_network.graph.html" =
                         panSEA.BeatAML$mGSEA.network$interactive)
    all.files <- list('Differential expression' = DEG.files,
                      'GSEA' = GSEA.files,
                      'DMEA' = DMEA.files)
  }
  
  for (i in 1:length(all.files)) {
    # create local folder for subset of results
    setwd(paste0(base.path, synapse_id_map[k], "analysis"))
    dir.create(names(all.files)[i])
    setwd(names(all.files)[i])
    
    
    # save results locally
    temp.files <- all.files[[i]]
    
    CSV.files <- names(temp.files)[grepl(".csv", names(temp.files))]
    for (j in 1:length(CSV.files)) {
      write.csv(temp.files[[CSV.files[j]]], CSV.files[j], row.names = FALSE)
    }
    
    PDF.files <- names(temp.files)[grepl(".pdf", names(temp.files))]
    for (j in 1:length(PDF.files)) {
      if (grepl("dot", PDF.files[j])) { # wide plot for dot plots
        ggplot2::ggsave(PDF.files[j], temp.files[[PDF.files[j]]],
                        device = "pdf", width = 11)
      } else {
        ggplot2::ggsave(PDF.files[j], temp.files[[PDF.files[j]]],
                        device = "pdf")
      }
      # ggplot2::ggsave(PDF.files[j], temp.files[[PDF.files[j]]],
      #                 device = "pdf")
    }
    
    HTML.files <- names(temp.files)[grepl(".html", names(temp.files))]
    if (length(HTML.files) > 0) {
      for (j in 1:length(HTML.files)) {
        visNetwork::visSave(temp.files[[HTML.files[j]]], HTML.files[j])
      }
    }
    
    # create folder on Synpase for subset of results
    dataFolder <- 
        synapser::synStore(synapser::Folder(names(all.files)[i],
                                            parent = names(synapse_id_map)[k]))
    
    # upload results to Synapse
    CSVs <- lapply(as.list(CSV.files), synapser::File,
                    parent = dataFolder)
    lapply(CSVs, synapser::synStore)
    
    PDFs <- lapply(as.list(PDF.files), synapser::File,
                   parent = dataFolder)
    lapply(PDFs, synapser::synStore)
    
    if (length(HTML.files) > 0) {
      HTMLs <- lapply(HTML.files, synapser::File,
                                    parent = dataFolder)
      lapply(HTMLs, synapser::synStore) 
    }
  }
  if (omics[k] != "phospho") {
    setwd(paste0(base.path, synapse_id_map[k], "analysis"))
    saveRDS(panSEA.BeatAML, file=paste0("exp22_", omics[k], "_panSEA_BeatAML.rds")) # 1.09 GB for 3 contrasts
    #panSEA.BeatAML <- readRDS(paste0("exp22_", omics[k], "_panSEA_BeatAML.rds")) 
  }
  panSEA.BeatAML <- NULL # make space to process next omics type
}
