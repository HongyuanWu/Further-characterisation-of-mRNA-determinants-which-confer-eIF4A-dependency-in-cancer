#VIOLIN PLOTS_WaALDRON VERSUS MODELSKA
#This script has been created by Rocco Roberto Penna with the support of Dr. Joseph A. Waldron

#AIM: to compare mRNA features generally linked to eIF4A-dependency with eIF4A-antidependent data
#OBJECTIVE: to compare 5'UTR length, GC content, G content and C content between Iwasaki et al. 2016 data sets (hippuristal and RocA) 


#R packages---
library(dplyr)
library(tidyverse)
library(ggplot2)


#Data import---
##Transcript features
fpUTR <- read.csv("Files/gencode_v27_pc_transcripts_fpUTRs_composition.csv", header = T)
fpUTR %>%
  select(transcript, GC_content, length, G_content, C_content) -> fpUTR #Isolating 5'UTR length, GC, G & C contents of all transcripts

CDS <- read.csv("Files/gencode_v27_pc_transcripts_CDSs_composition.csv", header = T)
CDS %>%
  select(transcript, GC_content, length, G_content, C_content) -> CDS #Isolating CDS length, GC, G & C contents of all transcripts

tpUTR <- read.csv("Files/gencode_v27_pc_transcripts_tpUTRs_composition.csv", header = T)
tpUTR %>%
  select(transcript, GC_content, length, G_content, C_content) -> tpUTR #Isolating 3'UTR length, GC, G & C contents of all transcripts

##Iwasaki eIF4A dependent data (Hippuristanol and RocA)
source("III_Box and Violin plots/Box Plots_Sources_Dependent.R") #Importing eIF4A-dep data (focus on Iwasaki data sets)


#Creating overlapping group---
overlap <- intersect(iwa_hipp_final_d, iwa_roca_final_d) #Creating overlap
iwaH_only <- setdiff(iwa_hipp_final_d, overlap) #Isolating Iwasaki-only transcripts (hippuristanol)
iwaR_only <- setdiff(iwa_roca_final_d, overlap) #Isolating Iwasaki-only transcripts (RocA)


#5'UTR---
##Data manipulation
overlap %>%
  inner_join(fpUTR, by = c("transcript_v" = "transcript")) %>% #Coupling 5'UTR features
  mutate(group = rep("Both")) %>% #Adding unique label (overlap)
  select(group, GC_content, length, G_content, C_content) -> overlap_fp #Isolating transcripts

iwaH_only %>%
  inner_join(fpUTR, by = c("transcript_v" = "transcript")) %>%
  mutate(group = rep("Iwasaki (Hipp)\net al.")) %>% #Adding unique label (Iwasaki hippuristanol)
  select(group, GC_content, length, G_content, C_content) -> iwaH_fp

iwaR_only %>%
  inner_join(fpUTR, by = c("transcript_v" = "transcript")) %>%
  mutate(group = rep("Iwasaki (RocA)\net al.")) %>% #Adding unique label (Iwasaki RocA)
  select(group, GC_content, length, G_content, C_content) -> iwaR_fp

HR_fp <- rbind(iwaH_fp, overlap_fp, iwaR_fp) #Merging data
HR_fp %>%
  mutate(group = fct_relevel(group, "Iwasaki (Hipp)\net al.", "Both", "Iwasaki (RocA)\net al.")) -> HR_fp #Ordering

##Statistical test function (ANOVA)
#Funtion of label printing
myP <- function(p) {
  if (p < 2.2e-16) {
    p_label <- "P < 2.2e-16"
  } else {
    if (p < 0.001) {
      rounded_p <- formatC(p, format = "e", digits = 2)
    } else {
      rounded_p <- round(p, digits = 3)
    }
    p_label <- paste("P =", rounded_p)
  }
  return(paste0(p_label))
}

#One-way Anova and Tukey multiple pairwise-comparisons 
res.aov_fp_GC_hr <- aov(GC_content ~ group, data = HR_fp) #GC content
summary(res.aov_fp_GC_hr)
Tuk_fp_GC_HR <- TukeyHSD(res.aov_fp_GC_hr)
Tuk_fp_GC_HR
p_Hb_fp_GC <- myP(Tuk_fp_GC_HR$group[10])
p_HR_fp_GC <- myP(Tuk_fp_GC_HR$group[11])
p_Rb_fp_GC <- myP(Tuk_fp_GC_HR$group[12])

res.aov_fp_l_hr <- aov(log10(length) ~ group, data = HR_fp) #length
summary(res.aov_fp_l_hr)
Tuk_fp_l_HR <- TukeyHSD(res.aov_fp_l_hr)
Tuk_fp_l_HR
p_Hb_fp_l <- myP(Tuk_fp_l_HR$group[10])
p_HR_fp_l <- myP(Tuk_fp_l_HR$group[11])
p_Rb_fp_l <- myP(Tuk_fp_l_HR$group[12])

res.aov_fp_G_hr <- aov(G_content ~ group, data = HR_fp) #G content
summary(res.aov_fp_G_hr)
Tuk_fp_G_HR <- TukeyHSD(res.aov_fp_G_hr)
Tuk_fp_G_HR
p_Hb_fp_G <- myP(Tuk_fp_G_HR$group[10])
p_HR_fp_G <- myP(Tuk_fp_G_HR$group[11])
p_Rb_fp_G <- myP(Tuk_fp_G_HR$group[12])

res.aov_fp_C_hr <- aov(C_content ~ group, data = HR_fp) #C content
summary(res.aov_fp_C_hr)
Tuk_fp_C_HR <- TukeyHSD(res.aov_fp_C_hr)
Tuk_fp_C_HR
p_Hb_fp_C <- myP(Tuk_fp_C_HR$group[10])
p_HR_fp_C <- myP(Tuk_fp_C_HR$group[11])
p_Rb_fp_C <- myP(Tuk_fp_C_HR$group[12])


#CDS---
##Data manipulation
overlap %>%
  inner_join(CDS, by = c("transcript_v" = "transcript")) %>% #Coupling CDS features
  mutate(group = rep("Both")) %>%
  select(group, GC_content, length, G_content, C_content) -> overlap_cds

iwaH_only %>%
  inner_join(CDS, by = c("transcript_v" = "transcript")) %>%
  mutate(group = rep("Iwasaki (Hipp)\net al.")) %>%
  select(group, GC_content, length, G_content, C_content) -> iwaH_cds

iwaR_only %>%
  inner_join(CDS, by = c("transcript_v" = "transcript")) %>%
  mutate(group = rep("Iwasaki (RocA)\net al.")) %>%
  select(group, GC_content, length, G_content, C_content) -> iwaR_cds

HR_cds <- rbind(iwaH_cds, overlap_cds, iwaR_cds)
HR_cds %>%
  mutate(group = fct_relevel(group, "Iwasaki (Hipp)\net al.", "Both", "Iwasaki (RocA)\net al.")) -> HR_cds

##Statistical test function (ANOVA)
myP <- function(p) {
  if (p < 2.2e-16) {
    p_label <- "P < 2.2e-16"
  } else {
    if (p < 0.001) {
      rounded_p <- formatC(p, format = "e", digits = 2)
    } else {
      rounded_p <- round(p, digits = 3)
    }
    p_label <- paste("P =", rounded_p)
  }
  return(paste0(p_label))
}

#One-way Anova and Tukey multiple pairwise-comparisons 
res.aov_cds_GC_hr <- aov(GC_content ~ group, data = HR_cds) #GC content
summary(res.aov_cds_GC_hr)
Tuk_cds_GC_HR <- TukeyHSD(res.aov_cds_GC_hr)
Tuk_cds_GC_HR
p_Hb_cds_GC <- myP(Tuk_cds_GC_HR$group[10])
p_HR_cds_GC <- myP(Tuk_cds_GC_HR$group[11])
p_Rb_cds_GC <- myP(Tuk_cds_GC_HR$group[12])

res.aov_cds_l_hr <- aov(log10(length) ~ group, data = HR_cds) #length
summary(res.aov_cds_l_hr)
Tuk_cds_l_HR <- TukeyHSD(res.aov_cds_l_hr)
Tuk_cds_l_HR
p_Hb_cds_l <- myP(Tuk_cds_l_HR$group[10])
p_HR_cds_l <- myP(Tuk_cds_l_HR$group[11])
p_Rb_cds_l <- myP(Tuk_cds_l_HR$group[12])

res.aov_cds_G_hr <- aov(G_content ~ group, data = HR_cds) #G content
summary(res.aov_cds_G_hr)
Tuk_cds_G_HR <- TukeyHSD(res.aov_cds_G_hr)
Tuk_cds_G_HR
p_Hb_cds_G <- myP(Tuk_cds_G_HR$group[10])
p_HR_cds_G <- myP(Tuk_cds_G_HR$group[11])
p_Rb_cds_G <- myP(Tuk_cds_G_HR$group[12])

res.aov_cds_C_hr <- aov(C_content ~ group, data = HR_cds) #C content
summary(res.aov_cds_C_hr)
Tuk_cds_C_HR <- TukeyHSD(res.aov_cds_C_hr)
Tuk_cds_C_HR
p_Hb_cds_C <- myP(Tuk_cds_C_HR$group[10])
p_HR_cds_C <- myP(Tuk_cds_C_HR$group[11])
p_Rb_cds_C <- myP(Tuk_cds_C_HR$group[12])


#3' UTR
##Data manipulation
overlap %>%
  inner_join(tpUTR, by = c("transcript_v" = "transcript")) %>% #Coupling 3'UTR features
  mutate(group = rep("Both")) %>%
  select(group, GC_content, length, G_content, C_content) -> overlap_tp

iwaH_only %>%
  inner_join(tpUTR, by = c("transcript_v" = "transcript")) %>%
  mutate(group = rep("Iwasaki (Hipp)\net al.")) %>%
  select(group, GC_content, length, G_content, C_content) -> iwaH_tp

iwaR_only %>%
  inner_join(tpUTR, by = c("transcript_v" = "transcript")) %>%
  mutate(group = rep("Iwasaki (RocA)\net al.")) %>%
  select(group, GC_content, length, G_content, C_content) -> iwaR_tp

HR_tp <- rbind(iwaH_tp, overlap_tp, iwaR_tp)
HR_tp %>%
  mutate(group = fct_relevel(group, "Iwasaki (Hipp)\net al.", "Both", "Iwasaki (RocA)\net al.")) -> HR_tp

##Statistical test function (ANOVA)
myP <- function(p) {
  if (p < 2.2e-16) {
    p_label <- "P < 2.2e-16"
  } else {
    if (p < 0.001) {
      rounded_p <- formatC(p, format = "e", digits = 2)
    } else {
      rounded_p <- round(p, digits = 3)
    }
    p_label <- paste("P =", rounded_p)
  }
  return(paste0(p_label))
}

#One-way Anova and Tukey multiple pairwise-comparisons 
res.aov_tp_GC_hr <- aov(GC_content ~ group, data = HR_tp) #GC content
summary(res.aov_tp_GC_hr)
Tuk_tp_GC_HR <- TukeyHSD(res.aov_tp_GC_hr)
Tuk_tp_GC_HR
p_Hb_tp_GC <- myP(Tuk_tp_GC_HR$group[10])
p_HR_tp_GC <- myP(Tuk_tp_GC_HR$group[11])
p_Rb_tp_GC <- myP(Tuk_tp_GC_HR$group[12])

res.aov_tp_l_hr <- aov(log10(length) ~ group, data = HR_tp) #length
summary(res.aov_tp_l_hr)
Tuk_tp_l_HR <- TukeyHSD(res.aov_tp_l_hr)
Tuk_tp_l_HR
p_Hb_tp_l <- myP(Tuk_tp_l_HR$group[10])
p_HR_tp_l <- myP(Tuk_tp_l_HR$group[11])
p_Rb_tp_l <- myP(Tuk_tp_l_HR$group[12])

res.aov_tp_G_hr <- aov(G_content ~ group, data = HR_tp) #G content
summary(res.aov_tp_G_hr)
Tuk_tp_G_HR <- TukeyHSD(res.aov_tp_G_hr)
Tuk_tp_G_HR
p_Hb_tp_G <- myP(Tuk_tp_G_HR$group[10])
p_HR_tp_G <- myP(Tuk_tp_G_HR$group[11])
p_Rb_tp_G <- myP(Tuk_tp_G_HR$group[12])

res.aov_tp_C_hr <- aov(C_content ~ group, data = HR_tp) #C content
summary(res.aov_tp_C_hr)
Tuk_tp_C_HR <- TukeyHSD(res.aov_tp_C_hr)
Tuk_tp_C_HR
p_Hb_tp_C <- myP(Tuk_tp_C_HR$group[10])
p_HR_tp_C <- myP(Tuk_tp_C_HR$group[11])
p_Rb_tp_C <- myP(Tuk_tp_C_HR$group[12])


#Plot---
##5'UTR
##GC content
GC_HR_fp <- ggplot(HR_fp, aes(x = group, y = GC_content, fill = group)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(alpha = 1, linetype = 2) +
  ylim(c(0.00, 1.00)) +
  ylab("5'UTR GC content (%)") +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        panel.grid.minor = element_blank())+
  annotate("text", x = 1.5, y = 0.95, size=4, label = p_Hb_fp_GC)+
  annotate("text", x = 2.5, y = 0.95, size=4, label = p_Rb_fp_GC)+
  annotate("text", x = 2, y = 0.95, size=4, label = p_HR_fp_GC)

##length
length_HR_fp <- ggplot(HR_fp, aes(x = group, y = length, fill = group)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(alpha = 1, linetype = 2) +
  scale_y_log10(breaks=c(10,100,1000),limits=c(10, 1500))+
  ylab("5'UTR length (nt, log10)") +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        panel.grid.minor = element_blank())+
  annotate("text", x = 1.5, y = 1500, size=4, label = p_Hb_fp_l)+
  annotate("text", x = 2.5, y = 1500, size=4, label = p_Rb_fp_l)+
  annotate("text", x = 2, y = 1500, size=4, label = p_HR_fp_l)

##G content
G_HR_fp <- ggplot(HR_fp, aes(x = group, y = G_content, fill = group)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(alpha = 1, linetype = 2) +
  ylim(c(0.00, 1.00)) +
  ylab("5'UTR G content (%)") +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        panel.grid.minor = element_blank())+
  annotate("text", x = 1.5, y = 0.95, size=4, label = p_Hb_fp_G)+
  annotate("text", x = 2.5, y = 0.95, size=4, label = p_Rb_fp_G)+
  annotate("text", x = 2, y = 0.95, size=4, label = p_HR_fp_G)

##C content
C_HR_fp <- ggplot(HR_fp, aes(x = group, y = C_content, fill = group)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(alpha = 1, linetype = 2) +
  ylim(c(0.00, 1.00)) +
  ylab("5'UTR C content (%)") +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        panel.grid.minor = element_blank())+
  annotate("text", x = 1.5, y = 0.95, size=4, label = p_Hb_fp_C)+
  annotate("text", x = 2.5, y = 0.95, size=4, label = p_Rb_fp_C)+
  annotate("text", x = 2, y = 0.95, size=4, label = p_HR_fp_C)

##CDS
##GC content
GC_HR_cds <- ggplot(HR_cds, aes(x = group, y = GC_content, fill = group)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(alpha = 1, linetype = 2) +
  ylim(c(0.00, 1.00)) +
  ylab("CDS GC content (%)") +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        panel.grid.minor = element_blank())+
  annotate("text", x = 1.5, y = 0.95, size=4, label = p_Hb_cds_GC)+
  annotate("text", x = 2.5, y = 0.95, size=4, label = p_Rb_cds_GC)+
  annotate("text", x = 2, y = 0.95, size=4, label = p_HR_cds_GC)

##length
length_HR_cds <- ggplot(HR_cds, aes(x = group, y = length, fill = group)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(alpha = 1, linetype = 2) +
  scale_y_log10(breaks=c(10,100,1000),limits=c(10, 1500))+
  ylab("CDS length (nt, log10)") +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        panel.grid.minor = element_blank())+
  annotate("text", x = 1.5, y = 1500, size=4, label = p_Hb_cds_l)+
  annotate("text", x = 2.5, y = 1500, size=4, label = p_Rb_cds_l)+
  annotate("text", x = 2, y = 40, size=4, label = p_HR_cds_l)

##G content
G_HR_cds <- ggplot(HR_cds, aes(x = group, y = G_content, fill = group)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(alpha = 1, linetype = 2) +
  ylim(c(0.00, 1.00)) +
  ylab("CDS G content (%)") +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        panel.grid.minor = element_blank())+
  annotate("text", x = 1.5, y = 0.95, size=4, label = p_Hb_cds_G)+
  annotate("text", x = 2.5, y = 0.95, size=4, label = p_Rb_cds_G)+
  annotate("text", x = 2, y = 0.95, size=4, label = p_HR_cds_G)

##C content
C_HR_cds <- ggplot(HR_cds, aes(x = group, y = C_content, fill = group)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(alpha = 1, linetype = 2) +
  ylim(c(0.00, 1.00)) +
  ylab("CDS C content (%)") +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        panel.grid.minor = element_blank())+
  annotate("text", x = 1.5, y = 0.95, size=4, label = p_Hb_cds_C)+
  annotate("text", x = 2.5, y = 0.95, size=4, label = p_Rb_cds_C)+
  annotate("text", x = 2, y = 0.95, size=4, label = p_HR_cds_C)

##3'UTR
##GC content
GC_HR_tp <- ggplot(HR_tp, aes(x = group, y = GC_content, fill = group)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(alpha = 1, linetype = 2) +
  ylim(c(0.00, 1.00)) +
  ylab("3'UTR GC content (%)") +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        panel.grid.minor = element_blank())+
  annotate("text", x = 1.5, y = 0.95, size=4, label = p_Hb_tp_GC)+
  annotate("text", x = 2.5, y = 0.95, size=4, label = p_Rb_tp_GC)+
  annotate("text", x = 2, y = 0.95, size=4, label = p_HR_tp_GC)

##length
length_HR_tp <- ggplot(HR_tp, aes(x = group, y = length, fill = group)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(alpha = 1, linetype = 2) +
  scale_y_log10(breaks=c(10,100,1000),limits=c(10, 1500))+
  ylab("3'UTR length (nt, log10)") +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        panel.grid.minor = element_blank())+
  annotate("text", x = 1.5, y = 1500, size=4, label = p_Hb_tp_l)+
  annotate("text", x = 2.5, y = 1500, size=4, label = p_Rb_tp_l)+
  annotate("text", x = 2, y = 1500, size=4, label = p_HR_tp_l)

##G content
G_HR_tp <- ggplot(HR_tp, aes(x = group, y = G_content, fill = group)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(alpha = 1, linetype = 2) +
  ylim(c(0.00, 1.00)) +
  ylab("3'UTR G content (%)") +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        panel.grid.minor = element_blank())+
  annotate("text", x = 1.5, y = 0.95, size=4, label = p_Hb_tp_G)+
  annotate("text", x = 2.5, y = 0.95, size=4, label = p_Rb_tp_G)+
  annotate("text", x = 2, y = 0.95, size=4, label = p_HR_tp_G)

##C content
C_HR_tp <- ggplot(HR_tp, aes(x = group, y = C_content, fill = group)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(alpha = 1, linetype = 2) +
  ylim(c(0.00, 1.00)) +
  ylab("3'UTR C content (%)") +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        panel.grid.minor = element_blank())+
  annotate("text", x = 1.5, y = 0.95, size=4, label = p_Hb_tp_C)+
  annotate("text", x = 2.5, y = 0.95, size=4, label = p_Rb_tp_C)+
  annotate("text", x = 2, y = 0.95, size=4, label = p_HR_tp_C)

