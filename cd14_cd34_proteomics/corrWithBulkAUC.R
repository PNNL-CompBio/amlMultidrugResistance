# correlations between sorted expression and bulk AUC
library(synapser);library(ggplot2);library(DMEA)
library(plyr);library(dplyr)
setwd("~/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/analysis")
synapser::synLogin()
base.path <- getwd()

#### prep data ####
# get sorted data
dia.wo.out <- readRDS(synapser::synGet("syn70198439")$path) # DIA_2batches_noOutliers.rds
# if you want to use the TMT data instead, that RDS file is here: syn70198440

# get bulk drug sensitivity data
drug.BeatAML <- read.csv(synapser::synGet("syn65677650")$path) # other drug sensitivity data: syn51674470
ven <- drug.BeatAML[drug.BeatAML$Specimen..Lab.ID %in% dia.wo.out$meta$patient & 
                      drug.BeatAML$Inhibitor.Panel.Definition..Drug=="Venetoclax",
                            c("Specimen..Lab.ID","Probit.Interpretation..Area.Under.Curve")] # 20 / 23
azaVen <- drug.BeatAML[drug.BeatAML$Specimen..Lab.ID %in% dia.wo.out$meta$patient & 
                      drug.BeatAML$Inhibitor.Panel.Definition..Drug=="Azacitidine - Venetoclax",
                    c("Specimen..Lab.ID","Probit.Interpretation..Area.Under.Curve")] # 20 / 23
# any samples lost?
lost.samples.ven <- unique(dia.wo.out$meta$patient[!(dia.wo.out$meta$patient %in% ven$Specimen..Lab.ID)])
# 3 / 23 patients lost: "16-01184" "21-00839" "22-00251"
lost.samples.av <- unique(dia.wo.out$meta$patient[!(dia.wo.out$meta$patient %in% azaVen$Specimen..Lab.ID)])
# 3 / 23 patients lost: "16-01184" "21-00839" "22-00251"

# add data from another file to get info for these 3 patients
other.pts <- read.csv(synapser::synGet("syn53627410")$path)
colnames(ven) <- c("sample_id","auc")
colnames(azaVen) <- c("sample_id","auc")
other.ven <- other.pts[other.pts$patient %in% lost.samples.ven,c("patient","Ven_AUC")]
other.av <- other.pts[other.pts$patient %in% lost.samples.av,c("patient","Aza.Ven_AUC")]
colnames(other.ven) <- colnames(ven)
colnames(other.av) <- colnames(ven)
ven <- rbind(ven, other.ven)
azaVen <- rbind(azaVen, other.av)

# format proteomics for correlations
prot <- as.data.frame(t(dia.wo.out$global))
prot$sample_id <- rownames(prot)
prot <- prot[,c("sample_id",colnames(prot)[1:(ncol(prot)-1)])]
prot$sample_id <- sub("PTRC_","",prot$sample_id)
prot$sample_id <- sub("_.*","",prot$sample_id)
prot$sample_id <- sub("[.]","-",prot$sample_id)
id.map <- unique(prot$sample_id[grepl("X",prot$sample_id)])
names(id.map) <- c("16-01184","17-01060","18-00103","18-00105","19-00074",
                   "21-00432","21-00839","22-00117","22-00251","22-00571")
for (i in 1:length(id.map)) {
  prot[prot$sample_id==id.map[i],]$sample_id <- names(id.map)[i]
}

# filter proteomics for sorting method, cell type
prot.bead <- prot[!grepl("_f_",rownames(prot)) & 
                    !grepl("Flow",rownames(prot), ignore.case=TRUE),]
prot.flow <- prot[grepl("_f_",rownames(prot)) | 
                    grepl("Flow",rownames(prot), ignore.case=TRUE),]
prot.cd14 <- prot[grepl("cd14",rownames(prot),ignore.case=TRUE),]
prot.cd14.bead <- prot.bead[grepl("cd14",rownames(prot.bead),ignore.case=TRUE),]
prot.cd14.flow <- prot.flow[grepl("cd14",rownames(prot.flow),ignore.case=TRUE),]
prot.cd34 <- prot[grepl("cd34",rownames(prot),ignore.case=TRUE),]
prot.cd34.bead <- prot.bead[grepl("cd34",rownames(prot.bead),ignore.case=TRUE),]
prot.cd34.flow <- prot.flow[grepl("cd34",rownames(prot.flow),ignore.case=TRUE),]
prot.msc.flow <- prot.flow[grepl("msc",row.names(prot.flow),ignore.case=TRUE),]

#### run correlations ####
dir.create("correlations")
setwd("correlations")
inputs <- list("Overall" = prot, "Bead" = prot.bead, "Flow" = prot.flow,
               "CD14" = prot.cd14, "CD14_Bead" = prot.cd14.bead,
               "CD14_Flow" = prot.cd14.flow, "CD34" = prot.cd34,
               "CD34_Bead" = prot.cd34.bead, "CD34_Flow" = prot.cd34.flow,
               "MSC_Flow" = prot.msc.flow)
for (i in names(inputs)) {
  # prep data
  ven.prot <- merge(ven, inputs[[i]], by="sample_id")
  ven.prot$auc <- as.numeric(ven.prot$auc)
  ven.prot <- ven.prot[!is.na(ven.prot$auc),]
  av.prot <- merge(azaVen, inputs[[i]], by="sample_id")
  av.prot$auc <- as.numeric(av.prot$auc)
  av.prot <- av.prot[!is.na(av.prot$auc),]
  
  # run correlation
  ven.prot.corr <- DMEA::rank_corr(ven.prot, variable="Protein",value="normAbudance") # no q<0.05
  av.prot.corr <- DMEA::rank_corr(av.prot, variable="Protein",value="normAbudance") # no q<0.05
  
  # save results
  write.csv(ven.prot.corr$result,
            paste0("venAUC_correlationsWithSortedProteomicsDIA_",i,".csv"),row.names=FALSE)
  write.csv(av.prot.corr$result,
            paste0("azaVenAUC_correlationsWithSortedProteomicsDIA_",i,".csv"),row.names=FALSE)
  if (is.list(ven.prot.corr$scatter.plots)) {
    ggsave(paste0("venAUC_correlationsWithSortedProteomicsDIA_",i,".pdf"),
           ven.prot.corr$scatter.plots, width=5,height=5)
  }
  if (is.list(av.prot.corr$scatter.plots)) {
    ggsave(paste0("azaVenAUC_correlationsWithSortedProteomicsDIA_",i,".pdf"),
           av.prot.corr$scatter.plots, width=5,height=5)
  }
}

#### histograms ####
setwd(file.path(base.path,"correlations"))
dir.create("histograms")
setwd("histograms")
for (i in names(inputs)) {
  # prep data
  setwd(file.path(base.path,"correlations"))
  ven.prot.corr <- read.csv(paste0("venAUC_correlationsWithSortedProteomicsDIA_",i,".csv"))
  av.prot.corr <- read.csv(paste0("azaVenAUC_correlationsWithSortedProteomicsDIA_",i,".csv"))
  gsea.inputs <- list("Ven" = ven.prot.corr[,c("Protein","Spearman.est","Spearman.q")],
                      "AzaVen" = av.prot.corr[,c("Protein","Spearman.est","Spearman.q")])
  
  for (j in names(gsea.inputs)) {
    setwd(file.path(base.path,"correlations","histograms"))
    dir.create(j)
    setwd(j)
    temp.df <- na.omit(gsea.inputs[[j]])
    temp.df$Significance <- "Adjusted p > 0.05"
    if (any(temp.df$Spearman.q <= 0.05)) {
      temp.df[which(temp.df$Spearman.q <= 0.05),]$Significance <- "Adjusted p <= 0.05"
    }
    n.distinct <- length(unique(temp.df$Spearman.est))
    n.ties <- nrow(temp.df) - n.distinct
    perc.ties <- round(n.ties * 100 / nrow(temp.df),0)
    temp.plot <- ggplot2::ggplot(temp.df, aes(x=Spearman.est, 
                                              group=Significance, 
                                              fill=Significance)) +
      geom_histogram(position="identity", alpha=0.5) + theme_classic() +
      scale_fill_manual(values = c("#00BFC4", "#F8766D"), breaks = c("Adjusted p <= 0.05", "Adjusted p > 0.05")) +
      xlab("Spearman Correlation Estimates") + 
      ggtitle(paste0(j, " AUC correlations with DIA\n",i," Protein Expression\n(", perc.ties, "% tied)")) +
      theme(axis.text = element_text(size=12), axis.title = element_text(size=16), legend.title = element_text(size=16),
            legend.text=element_text(size=12), title = element_text(size = 24, hjust = 0.5),
            legend.position = "bottom")
    ggsave(paste0(j,"AUC_correlationsWithSortedProteomicsDIA_",i,".pdf"),
           temp.plot, width=7, height=7)
  }
}

#### run GSEA ####
setwd(file.path(base.path,"correlations"))
dir.create("GSEA")
setwd("GSEA")
for (i in names(inputs)) {
  # prep data
  setwd(file.path(base.path,"correlations"))
  ven.prot.corr <- read.csv(paste0("venAUC_correlationsWithSortedProteomicsDIA_",i,".csv"))
  av.prot.corr <- read.csv(paste0("azaVenAUC_correlationsWithSortedProteomicsDIA_",i,".csv"))
  gsea.inputs <- list("Ven" = ven.prot.corr[,c("Protein","Spearman.est")],
                      "AzaVen" = av.prot.corr[,c("Protein","Spearman.est")])
  
  # run GSEA with ties
  temp.gsea <- panSEA::mGSEA(gsea.inputs, 
                             feature.names = rep("Protein", length(gsea.inputs)),
                             rank.var = rep("Spearman.est", length(gsea.inputs)),
                             gmt=as.list(rep("msigdb_Homo sapiens_HS_H", length(gsea.inputs))),
                             ties=TRUE, types=names(gsea.inputs))
  
  # save results
  for (j in names(temp.gsea$all.results)) {
    setwd(file.path(base.path,"correlations","GSEA"))
    dir.create(j)
    setwd(j)
    write.csv(temp.gsea$all.results[[j]]$result,
              paste0(j,"AUC_gsea_WithSortedProteomicsDIA_",i,".csv"),row.names=FALSE)
    write.csv(temp.gsea$all.results[[j]]$result.w.ties,
              paste0(j,"AUC_gseaWithTies_WithSortedProteomicsDIA_",i,".csv"),row.names=FALSE)
    ggsave(paste0(j,"AUC_gseaVolcano_WithSortedProteomicsDIA_",i,".pdf"),
           temp.gsea$all.results[[j]]$volcano.plot, width=5, height=5)
    ggsave(paste0(j,"AUC_gseaBar_WithSortedProteomicsDIA_",i,".pdf"),
           temp.gsea$all.results[[j]]$bar.plot, width=5, height=5)
    ggsave(paste0(j,"AUC_gseaDot_WithSortedProteomicsDIA_",i,".pdf"),
           temp.gsea$all.results[[j]]$dot.plot, width=5, height=5)
    ggsave(paste0(j,"AUC_gseaDotSD_WithSortedProteomicsDIA_",i,".pdf"),
           temp.gsea$all.results[[j]]$dot.sd, width=5, height=5)
  }
}


#### redo plots ####
setwd(file.path(base.path,"correlations"))
dir.create("plots")
setwd("plots")
conditions <- c("Overall", "Bead", "Flow", "CD14", "CD14_Bead", "CD14_Flow", 
                "CD34", "CD34_Bead", "CD34_Flow", "MSC_Flow")
n.vals <- c(10,12,15)
for (i in conditions) {
  # get correlation results
  setwd(file.path(base.path,"correlations"))
  ven.prot.corr <- read.csv(paste0("venAUC_correlationsWithSortedProteomicsDIA_",i,".csv"))
  av.prot.corr <- read.csv(paste0("azaVenAUC_correlationsWithSortedProteomicsDIA_",i,".csv"))
  gsea.inputs <- list("Ven" = ven.prot.corr,
                      "AzaVen" = av.prot.corr)
  
  # bar plots
  for (j in names(gsea.inputs)) {
    setwd(file.path(base.path,"correlations","plots"))
    dir.create(j)
    setwd(j)
    temp.df <- gsea.inputs[[j]]
    n.total <- nrow(temp.df)
    temp.df <- temp.df[temp.df$Spearman.q<0.05,c("Protein","Spearman.est")]
    n.sig <- nrow(temp.df)
    if (nrow(temp.df)>0) {
      temp.df$Direction <- "Positive"
      if(any(temp.df$Spearman.est<0)){temp.df[temp.df$Spearman.est<0,]$Direction <- "Negative"}
      for (n in n.vals) {
        temp.df <- temp.df %>% slice_max(abs(Spearman.est), n=n)
        bar <- ggplot2::ggplot(temp.df, 
                               aes(x=Spearman.est, 
                                   y=reorder(Protein, Spearman.est), 
                                   fill = Direction)) + 
          geom_bar(stat='identity') + ggplot2::theme(
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
            axis.line = element_line(colour = "black", linewidth = 0.65),
            legend.text = element_text(size = 20),
            axis.text = element_text(size = 20),
            axis.title = element_text(size = 26, face = "bold"),
            panel.background = element_rect(
              fill = "white", colour = "white", linewidth = 0.5,
              linetype = "solid", color = "black"
            ), text = element_text(size = 20),
            legend.position = "bottom", legend.key = element_blank()
          ) + ggplot2::xlab("Spearman Correlation") + 
          theme(axis.title.y=element_blank()) + 
          ggtitle(paste(n.sig,"/", n.total, "proteins correlated")) +
          ggplot2::scale_fill_manual(
            values = c("red", "blue"), name = "Direction",
            breaks = c("Positive","Negative")
          ) 
        ggsave(paste0(j,"AUC_correlationsWithSortedProteomicsDIA_",i,"_top",n,"Sig_barPlot.pdf"),
               bar, width=7, height=7)
      }
    }
  }
  
  # get GSEA results
  setwd(file.path(base.path,"correlations","GSEA"))
  ven.gsea <- read.csv(paste0("Ven/VenAUC_gsea_WithSortedProteomicsDIA_",i,".csv"))
  av.gsea <- read.csv(paste0("AzaVen/AzaVenAUC_gsea_WithSortedProteomicsDIA_",i,".csv"))
  gsea.list <- list("Ven" = ven.gsea,
                      "AzaVen" = av.gsea)
  
  # dot plots
  for (j in names(gsea.list)) {
    setwd(file.path(base.path,"correlations","plots"))
    dir.create(j)
    setwd(j)
    bar.data <- gsea.list[[j]]
    n.total <- nrow(gsea.list[[j]])
    bar.data <- bar.data[bar.data$p_value<0.05 & bar.data$FDR_q_value<0.25,
                         c("Feature_set","NES","FDR_q_value")]
    n.sig <- nrow(bar.data)
    if (nrow(bar.data)>0) {
      for (n in n.vals) {
        bar.data <- bar.data %>% slice_max(abs(NES),n=n)
        bar.data$Gene_set <- sub("HALLMARK_","",bar.data$Feature_set)
        bar.data$minusLogFDR <- 4
        if (any(bar.data$FDR_q_value!=0)) {
          bar.data[bar.data$FDR_q_value!=0,]$minusLogFDR <- 
            -log10(bar.data[bar.data$FDR_q_value!=0,]$FDR_q_value)
        }
        dot.plt <- ggplot2::ggplot(bar.data,
                                   ggplot2::aes(x = j, y = reorder(Gene_set,NES), color = NES,
                                                size = minusLogFDR)) +
          ggplot2::geom_point() + scale_color_gradient2(low="blue",mid="grey",high="red")+
          scale_size_continuous(breaks=c(2,3,4), range=c(3,5))+
          geom_point(data = bar.data, col = "black", stroke = 1.5, shape = 21) +
          theme_classic() + ggplot2::labs(y = "Hallmark Pathways",
                                          color = "NES", size = "-log(FDR)") + 
          theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
          ggtitle(paste(n.sig,"/", n.total, "Hallmark\npathways enriched"))
        ggsave(paste0(j,"AUC_gseaDot_WithSortedProteomicsDIA_",i,"_top",n,"Sig.pdf"),
               dot.plt, width=4, height=4)  
      }
    }
  }
}


