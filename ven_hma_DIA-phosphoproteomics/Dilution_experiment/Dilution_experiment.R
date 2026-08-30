library(ggplot2)
library(dplyr)
library(tidyr)
library(openxlsx)
library(stringr)
library(ggside)

# load dilution experiment phosphosite level data generated through in silico, combined, or hybrid searches
df <- read.xlsx("PTRC_insilico_site_human.xlsx", check.names=FALSE)

## make sure columns containing intensity values are in numeric
df[,8:22] <-lapply(df[,8:22], as.numeric) 

## different uniprot protein accession IDs could be mapped to the same gene name. Create a new SITE2 column that is formatted as GeneName-Residue#, e.g.TADA2A-S6
df <- df %>% mutate(SITE = str_c(Gene.Names, "-", Residue, Site))

## subset site and sample columns
df <- df[,7:22]
## sum rows with same SITE ID
df <- df %>% group_by(SITE) %>%
  summarise(SITE=dplyr::first(SITE),
            across(everything(), sum, na.rm=TRUE))

## log2 transform
df[df == 0] <- NA
df[,2:16] <- log(df[, 2:16], 2)  
## check data distribution
boxplot(df[,2:16], cex.axis=1, las=2)

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

df[,2:16] <- Zero_Center_Norm(df[,2:16])

## check data distribution after median centering
boxplot(df[,2:16], cex.axis=1, las=2)

## subset 3 replicates from each dilution point. Focus on sites identified across all 3 replicates. Calculate the mean across replicates per site. 
df_D1 <- df[,c(1,2,7,12)]
df_D1$CountD1 <- apply(df_D1[,2:4], 1, function(x) sum(is.na(x)))
df_D1 <- filter(df_D1, CountD1 == 0)
df_D1$D1 <- rowMeans(df_D1[,2:4])


df_D2 <- df[,c(1,3,8,13)]
df_D2$CountD2 <- apply(df_D2[,2:4], 1, function(x) sum(is.na(x)))
df_D2 <- filter(df_D2, CountD2 ==0)
df_D2$D2 <- rowMeans(df_D2[,2:4])

df_D3 <- df[,c(1,4,9,14)]
df_D3$CountD3 <- apply(df_D3[,2:4], 1, function(x) sum(is.na(x)))
df_D3 <- filter(df_D3, CountD3 == 0)
df_D3$D3 <- rowMeans(df_D3[,2:4])

df_D4 <- df[,c(1,5,10,15)]
df_D4$CountD4 <- apply(df_D4[,2:4], 1, function(x) sum(is.na(x)))
df_D4 <- filter(df_D4, CountD4 == 0)
df_D4$D4 <- rowMeans(df_D4[,2:4])

df_D5 <- df[,c(1,6,11,16)]
df_D5$CountD5 <- apply(df_D5[,2:4], 1, function(x) sum(is.na(x)))
df_D5 <- filter(df_D5, CountD5 == 0)
df_D5$D5 <- rowMeans(df_D5[,2:4])

## merge sites identified across all dilution points.
M2 <- merge(df_D1[,c(1,6)], df_D2[,c(1,6)], by="SITE")
M3 <- merge(M2, df_D3[,c(1,6)], by ="SITE")
M4 <- merge(M3, df_D4[,c(1,6)], by ="SITE")
M5 <- merge(M4, df_D5[,c(1,6)], by ="SITE")

## use middle dilution point 50:50 as reference. Calculate log2FC between each dilution and the 50:50 reference.
M5$logFC_D1 <- M5$D1-M5$D3
M5$logFC_D2 <- M5$D2-M5$D3
M5$logFC_D3 <- M5$D3-M5$D3
M5$logFC_D4 <- M5$D4-M5$D3
M5$logFC_D5 <- M5$D5-M5$D3

# quick overview of log2FC
boxplot(M5[,7:11], cex.axis=1, las=2)

# pivot data longer for ggplot 
M5 <- M5[,c(1,7:11)]
dfp <- M5 %>% pivot_longer(cols= 2:6,
                           names_to="Dilution",
                           values_to="intensity")
## plot log2FC 
ggplot(dfp, aes(x = Dilution, y = intensity)) +
  geom_hline(yintercept=c(1,0.585,0,-1,-3.322), color=c("red","khaki4", "springgreen3", "deepskyblue","magenta"), linetype="dashed", linewidth=1)+
  geom_boxplot(
    aes(color=Dilution),
    fill=NA,
    position=position_dodge(width=0.6),
    outlier.shape=NA,
    width=0.5,
    lwd=1)+
  geom_point(
    aes(color=Dilution),
    position=position_jitterdodge(
      jitter.width = 0.3, 
      dodge.width = 0.6
    ), 
    alpha = 0.03, size = 0.3) +
  scale_x_discrete(
    name = "Dilution",
    labels= c("logFC_D1"= "100:0",
              "logFC_D2"= "75:25",
              "logFC_D3"="50:50",
              "logFC_D4"="25:75",
              "logFC_D5"="5:95"))+
  labs(y = "logFC",
       color = "Dilution"
  )+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border     = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color="black", size=1),
        axis.title.x=element_text(size=30, face="bold"),
        axis.title.y=element_text(size=30, face="bold"),
        axis.text.y=element_text(size=27, face="bold", color="black"),
        axis.text.x=element_text(size=30, face="bold", color="black"),
        legend.title=element_text(size=27, face="bold"),
        legend.text=element_text(size=20))

