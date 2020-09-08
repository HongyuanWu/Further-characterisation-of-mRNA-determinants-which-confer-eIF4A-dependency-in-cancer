#VIOLIN PLOTS_WaALDRON VERSUS MODELSKA
#This script has been created by Rocco Roberto Penna with the support of Dr. Joseph A. Waldron

#AIM: to compare mRNA features generally linked to eIF4A-dependency with eIF4A-antidependent data
#OBJECTIVE: to compare 5'UTR length, GC content, G content and C content between Waldron et al. 2019 and Modelska et al. 2015 eIF4A-dependent data


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

##Rubio and Wolfe eIF4A dependent data
source("III_Violin plots/Box Plots_Sources_Dependent.R") #Importing eIF4A-dep data (focus on Waldron and Modelska data)


#Creating overlapping group---
overlap <- intersect(wal_final_d, mod_final_d) #Creating overlap
wal_only <- setdiff(wal_final_d, overlap) #Isolating Waldron-only transcripts
mod_only <- setdiff(mod_final_d, overlap) #Isolating Modelska-only transcripts


#5'UTR---
##Data manipulation
overlap %>%
  inner_join(fpUTR, by = c("transcript_v" = "transcript")) %>% #Coupling 5'UTR features
  mutate(group = rep("Both")) %>% #Adding unique label (overlap)
  select(group, GC_content, length, G_content, C_content) -> overlap_fp #Isolating transcripts

wal_only %>%
  inner_join(fpUTR, by = c("transcript_v" = "transcript")) %>%
  mutate(group = rep("Waldron\net al.")) %>% #Adding unique label (Waldron)
  select(group, GC_content, length, G_content, C_content) -> wal_fp

mod_only %>%
  inner_join(fpUTR, by = c("transcript_v" = "transcript")) %>%
  mutate(group = rep("Modelska\net al.")) %>% #Adding unique label (Modelska)
  select(group, GC_content, length, G_content, C_content) -> mod_fp

WM_fp <- rbind(wal_fp, overlap_fp, mod_fp) #Mergin data
WM_fp %>%
  mutate(group = fct_relevel(group, "Waldron\net al.", "Both", "Modelska\net al.")) -> WM_fp #Ordering

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
res.aov_fp_GC_wm <- aov(GC_content ~ group, data = WM_fp) #GC content
summary(res.aov_fp_GC_wm)
Tuk_fp_GC_WM <- TukeyHSD(res.aov_fp_GC_wm)
Tuk_fp_GC_WM
p_Wb_fp_GC <- myP(Tuk_fp_GC_WM$group[10])
p_WM_fp_GC <- myP(Tuk_fp_GC_WM$group[11])
p_Mb_fp_GC <- myP(Tuk_fp_GC_WM$group[12])

res.aov_fp_l_wm <- aov(log10(length) ~ group, data = WM_fp) #length
summary(res.aov_fp_l_wm)
Tuk_fp_l_WM <- TukeyHSD(res.aov_fp_l_wm)
Tuk_fp_l_WM
p_Wb_fp_l <- myP(Tuk_fp_l_WM$group[10])
p_WM_fp_l <- myP(Tuk_fp_l_WM$group[11])
p_Mb_fp_l <- myP(Tuk_fp_l_WM$group[12])

res.aov_fp_G_wm <- aov(G_content ~ group, data = WM_fp) #G content
summary(res.aov_fp_G_wm)
Tuk_fp_G_WM <- TukeyHSD(res.aov_fp_G_wm)
Tuk_fp_G_WM
p_Wb_fp_G <- myP(Tuk_fp_G_WM$group[10])
p_WM_fp_G <- myP(Tuk_fp_G_WM$group[11])
p_Mb_fp_G <- myP(Tuk_fp_G_WM$group[12])

res.aov_fp_C_wm <- aov(C_content ~ group, data = WM_fp) #C content
summary(res.aov_fp_C_wm)
Tuk_fp_C_WM <- TukeyHSD(res.aov_fp_C_wm)
Tuk_fp_C_WM
p_Wb_fp_C <- myP(Tuk_fp_C_WM$group[10])
p_WM_fp_C <- myP(Tuk_fp_C_WM$group[11])
p_Mb_fp_C <- myP(Tuk_fp_C_WM$group[12])


#CDS---
##Data manipulation
overlap %>%
  inner_join(CDS, by = c("transcript_v" = "transcript")) %>% #Coupling CDS features
  mutate(group = rep("Both")) %>% 
  select(group, GC_content, length, G_content, C_content) -> overlap_cds

wal_only %>%
  inner_join(CDS, by = c("transcript_v" = "transcript")) %>%
  mutate(group = rep("Waldron\net al.")) %>%
  select(group, GC_content, length, G_content, C_content) -> wal_cds

mod_only %>%
  inner_join(CDS, by = c("transcript_v" = "transcript")) %>%
  mutate(group = rep("Modelska\net al.")) %>%
  select(group, GC_content, length, G_content, C_content) -> mod_cds

WM_cds <- rbind(wal_cds, overlap_cds, mod_cds)
WM_cds %>%
  mutate(group = fct_relevel(group, "Waldron\net al.", "Both", "Modelska\net al.")) -> WM_cds

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
res.aov_cds_GC_wm <- aov(GC_content ~ group, data = WM_cds) #GC content
summary(res.aov_cds_GC_wm)
Tuk_cds_GC_WM <- TukeyHSD(res.aov_cds_GC_wm)
Tuk_cds_GC_WM
p_Wb_cds_GC <- myP(Tuk_cds_GC_WM$group[10])
p_WM_cds_GC <- myP(Tuk_cds_GC_WM$group[11])
p_Mb_cds_GC <- myP(Tuk_cds_GC_WM$group[12])

res.aov_cds_l_wm <- aov(log10(length) ~ group, data = WM_cds) #length
summary(res.aov_cds_l_wm)
Tuk_cds_l_WM <- TukeyHSD(res.aov_cds_l_wm)
Tuk_cds_l_WM
p_Wb_cds_l <- myP(Tuk_cds_l_WM$group[10])
p_WM_cds_l <- myP(Tuk_cds_l_WM$group[11])
p_Mb_cds_l <- myP(Tuk_cds_l_WM$group[12])

res.aov_cds_G_wm <- aov(G_content ~ group, data = WM_cds) #G content
summary(res.aov_cds_G_wm)
Tuk_cds_G_WM <- TukeyHSD(res.aov_cds_G_wm)
Tuk_cds_G_WM
p_Wb_cds_G <- myP(Tuk_cds_G_WM$group[10])
p_WM_cds_G <- myP(Tuk_cds_G_WM$group[11])
p_Mb_cds_G <- myP(Tuk_cds_G_WM$group[12])

res.aov_cds_C_wm <- aov(C_content ~ group, data = WM_cds) #C content
summary(res.aov_cds_C_wm)
Tuk_cds_C_WM <- TukeyHSD(res.aov_cds_C_wm)
Tuk_cds_C_WM
p_Wb_cds_C <- myP(Tuk_cds_C_WM$group[10])
p_WM_cds_C <- myP(Tuk_cds_C_WM$group[11])
p_Mb_cds_C <- myP(Tuk_cds_C_WM$group[12])


#3' UTR
##Data manipulation
overlap %>%
  inner_join(tpUTR, by = c("transcript_v" = "transcript")) %>% #Coupling 3'UTR features
  mutate(group = rep("Both")) %>%
  select(group, GC_content, length, G_content, C_content) -> overlap_tp

wal_only %>%
  inner_join(tpUTR, by = c("transcript_v" = "transcript")) %>%
  mutate(group = rep("Waldron\net al.")) %>%
  select(group, GC_content, length, G_content, C_content) -> wal_tp

mod_only %>%
  inner_join(tpUTR, by = c("transcript_v" = "transcript")) %>%
  mutate(group = rep("Modelska\net al.")) %>%
  select(group, GC_content, length, G_content, C_content) -> mod_tp

WM_tp <- rbind(wal_tp, overlap_tp, mod_tp)
WM_tp %>%
  mutate(group = fct_relevel(group, "Waldron\net al.", "Both", "Modelska\net al.")) -> WM_tp

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
res.aov_tp_GC_wm <- aov(GC_content ~ group, data = WM_tp) #GC content
summary(res.aov_tp_GC_wm)
Tuk_tp_GC_WM <- TukeyHSD(res.aov_tp_GC_wm)
Tuk_tp_GC_WM
p_Wb_tp_GC <- myP(Tuk_tp_GC_WM$group[10])
p_WM_tp_GC <- myP(Tuk_tp_GC_WM$group[11])
p_Mb_tp_GC <- myP(Tuk_tp_GC_WM$group[12])

res.aov_tp_l_wm <- aov(log10(length) ~ group, data = WM_tp) #length
summary(res.aov_tp_l_wm)
Tuk_tp_l_WM <- TukeyHSD(res.aov_tp_l_wm)
Tuk_tp_l_WM
p_Wb_tp_l <- myP(Tuk_tp_l_WM$group[10])
p_WM_tp_l <- myP(Tuk_tp_l_WM$group[11])
p_Mb_tp_l <- myP(Tuk_tp_l_WM$group[12])

res.aov_tp_G_wm <- aov(G_content ~ group, data = WM_tp) #G content
summary(res.aov_tp_G_wm)
Tuk_tp_G_WM <- TukeyHSD(res.aov_tp_G_wm)
Tuk_tp_G_WM
p_Wb_tp_G <- myP(Tuk_tp_G_WM$group[10])
p_WM_tp_G <- myP(Tuk_tp_G_WM$group[11])
p_Mb_tp_G <- myP(Tuk_tp_G_WM$group[12])

res.aov_tp_C_wm <- aov(C_content ~ group, data = WM_tp) #C content
summary(res.aov_tp_C_wm)
Tuk_tp_C_WM <- TukeyHSD(res.aov_tp_C_wm)
Tuk_tp_C_WM
p_Wb_tp_C <- myP(Tuk_tp_C_WM$group[10])
p_WM_tp_C <- myP(Tuk_tp_C_WM$group[11])
p_Mb_tp_C <- myP(Tuk_tp_C_WM$group[12])


#Plot---
##5'UTR
##GC content
GC_WM_fp <- ggplot(WM_fp, aes(x = group, y = GC_content, fill = group)) +
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
  annotate("text", x = 1.5, y = 0.95, size=4, label = p_Wb_fp_GC)+
  annotate("text", x = 2.5, y = 0.95, size=4, label = p_Mb_fp_GC)+
  annotate("text", x = 2, y = 0.95, size=4, label = p_WM_fp_GC)

#length
length_WM_fp <- ggplot(WM_fp, aes(x = group, y = length, fill = group)) +
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
  annotate("text", x = 1.5, y = 1500, size=4, label = p_Wb_fp_l)+
  annotate("text", x = 2.5, y = 1500, size=4, label = p_Mb_fp_l)+
  annotate("text", x = 2, y = 1500, size=4, label = p_WM_fp_l)

#G content
G_WM_fp <- ggplot(WM_fp, aes(x = group, y = G_content, fill = group)) +
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
  annotate("text", x = 1.5, y = 0.95, size=4, label = p_Wb_fp_G)+
  annotate("text", x = 2.5, y = 0.95, size=4, label = p_Mb_fp_G)+
  annotate("text", x = 2, y = 0.95, size=4, label = p_WM_fp_G)

#C content
C_WM_fp <- ggplot(WM_fp, aes(x = group, y = C_content, fill = group)) +
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
  annotate("text", x = 1.5, y = 0.95, size=4, label = p_Wb_fp_C)+
  annotate("text", x = 2.5, y = 0.95, size=4, label = p_Mb_fp_C)+
  annotate("text", x = 2, y = 0.95, size=4, label = p_WM_fp_C)

##CDS
#GC content
GC_WM_cds <- ggplot(WM_cds, aes(x = group, y = GC_content, fill = group)) +
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
  annotate("text", x = 1.5, y = 0.95, size=4, label = p_Wb_cds_GC)+
  annotate("text", x = 2.5, y = 0.95, size=4, label = p_Mb_cds_GC)+
  annotate("text", x = 2, y = 0.95, size=4, label = p_WM_cds_GC)

##length
length_WM_cds <- ggplot(WM_cds, aes(x = group, y = length, fill = group)) +
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
  annotate("text", x = 1.5, y = 1500, size=4, label = p_Wb_cds_l)+
  annotate("text", x = 2.5, y = 1500, size=4, label = p_Mb_cds_l)+
  annotate("text", x = 2, y = 100, size=4, label = p_WM_cds_l)

##G content
G_WM_cds <- ggplot(WM_cds, aes(x = group, y = G_content, fill = group)) +
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
  annotate("text", x = 1.5, y = 0.95, size=4, label = p_Wb_cds_G)+
  annotate("text", x = 2.5, y = 0.95, size=4, label = p_Mb_cds_G)+
  annotate("text", x = 2, y = 0.95, size=4, label = p_WM_cds_G)

##C content
C_WM_cds <- ggplot(WM_cds, aes(x = group, y = C_content, fill = group)) +
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
  annotate("text", x = 1.5, y = 0.95, size=4, label = p_Wb_cds_C)+
  annotate("text", x = 2.5, y = 0.95, size=4, label = p_Mb_cds_C)+
  annotate("text", x = 2, y = 0.95, size=4, label = p_WM_cds_C)

##3'UTR
##GC content
GC_WM_tp <- ggplot(WM_tp, aes(x = group, y = GC_content, fill = group)) +
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
  annotate("text", x = 1.5, y = 0.95, size=4, label = p_Wb_tp_GC)+
  annotate("text", x = 2.5, y = 0.95, size=4, label = p_Mb_tp_GC)+
  annotate("text", x = 2, y = 0.95, size=4, label = p_WM_tp_GC)

##length
length_WM_tp <- ggplot(WM_tp, aes(x = group, y = length, fill = group)) +
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
  annotate("text", x = 1.5, y = 1500, size=4, label = p_Wb_tp_l)+
  annotate("text", x = 2.5, y = 1500, size=4, label = p_Mb_tp_l)+
  annotate("text", x = 2, y = 1500, size=4, label = p_WM_tp_l)

##G content
G_WM_tp <- ggplot(WM_tp, aes(x = group, y = G_content, fill = group)) +
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
  annotate("text", x = 1.5, y = 0.95, size=4, label = p_Wb_tp_G)+
  annotate("text", x = 2.5, y = 0.95, size=4, label = p_Mb_tp_G)+
  annotate("text", x = 2, y = 0.95, size=4, label = p_WM_tp_G)

##C content
C_WM_tp <- ggplot(WM_tp, aes(x = group, y = C_content, fill = group)) +
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
  annotate("text", x = 1.5, y = 0.95, size=4, label = p_Wb_tp_C)+
  annotate("text", x = 2.5, y = 0.95, size=4, label = p_Mb_tp_C)+
  annotate("text", x = 2, y = 0.95, size=4, label = p_WM_tp_C)

