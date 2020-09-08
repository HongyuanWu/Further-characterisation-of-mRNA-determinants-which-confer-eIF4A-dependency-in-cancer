#VIOLIN PLOTS_WaALDRON VERSUS MODELSKA
#This script has been created by Rocco Roberto Penna with the support of Dr. Joseph A. Waldron

#AIM: to compare mRNA features generally linked to eIF4A-dependency with eIF4A-antidependent data
#OBJECTIVE: to compare 5'UTR length, GC content, G content and C content between Rubio et al. 2014 and Wolfe et al. 2014 eIF4A-dependent data


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
source("III_Box and Violin plots/Box Plots_Sources_Dependent.R") #Importing eIF4A-dep data (focus on Rubio and Wolfe data)


#Creating overlapping group---
overlap <- intersect(rub_final_d, wolf_final_d) #Creating overlap
rub_only <- setdiff(rub_final_d, overlap) #Isolating Rubio-only transcripts
wolf_only <- setdiff(wolf_final_d, overlap) #Isolating Wolfe-only transcripts


#5'UTR---
##Data manipulation
overlap %>%
  inner_join(fpUTR, by = c("transcript_v" = "transcript")) %>% #Coupling 5'UTR features
  mutate(group = rep("Both")) %>% #Adding unique label (overlap)
  select(group, GC_content, length, G_content, C_content) -> overlap_fp #Isolating transcripts

rub_only %>%
  inner_join(fpUTR, by = c("transcript_v" = "transcript")) %>%
  mutate(group = rep("Rubio\net al.")) %>% #Adding unique label (Rubio)
  select(group, GC_content, length, G_content, C_content) -> rub_fp

wolf_only %>%
  inner_join(fpUTR, by = c("transcript_v" = "transcript")) %>%
  mutate(group = rep("Wolfe\net al.")) %>% #Adding unique label (Wolfe)
  select(group, GC_content, length, G_content, C_content) -> wolf_fp

RW_fp <- rbind(rub_fp, overlap_fp, wolf_fp) #Mergin data
RW_fp %>%
mutate(group = fct_relevel(group, "Rubio\net al.", "Both", "Wolfe\net al.")) -> RW_fp #Ordering

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
res.aov_fp_GC_rw <- aov(GC_content ~ group, data = RW_fp) #GC content
summary(res.aov_fp_GC_rw)
Tuk_fp_GC_RW <- TukeyHSD(res.aov_fp_GC_rw)
Tuk_fp_GC_RW
p_Rb_fp_GC <- myP(Tuk_fp_GC_RW$group[10])
p_WR_fp_GC <- myP(Tuk_fp_GC_RW$group[11])
p_Wb_fp_GC <- myP(Tuk_fp_GC_RW$group[12])

res.aov_fp_l_rw <- aov(log10(length) ~ group, data = RW_fp) #length
summary(res.aov_fp_l_rw)
Tuk_fp_l_RW <- TukeyHSD(res.aov_fp_l_rw)
Tuk_fp_l_RW
p_Rb_fp_l <- myP(Tuk_fp_l_RW$group[10])
p_WR_fp_l <- myP(Tuk_fp_l_RW$group[11])
p_Wb_fp_l <- myP(Tuk_fp_l_RW$group[12])

res.aov_fp_G_rw <- aov(G_content ~ group, data = RW_fp) #G content
summary(res.aov_fp_G_rw)
Tuk_fp_G_RW <- TukeyHSD(res.aov_fp_G_rw)
Tuk_fp_G_RW
p_Rb_fp_G <- myP(Tuk_fp_G_RW$group[10])
p_WR_fp_G <- myP(Tuk_fp_G_RW$group[11])
p_Wb_fp_G <- myP(Tuk_fp_G_RW$group[12])

res.aov_fp_C_rw <- aov(C_content ~ group, data = RW_fp) #C content
summary(res.aov_fp_C_rw)
Tuk_fp_C_RW <- TukeyHSD(res.aov_fp_C_rw)
Tuk_fp_C_RW
p_Rb_fp_C <- myP(Tuk_fp_C_RW$group[10])
p_WR_fp_C <- myP(Tuk_fp_C_RW$group[11])
p_Wb_fp_C <- myP(Tuk_fp_C_RW$group[12])


#CDS---
##Data manipulation
overlap %>%
  inner_join(CDS, by = c("transcript_v" = "transcript")) %>% #Coupling CDS features
  mutate(group = rep("Both")) %>%
  select(group, GC_content, length, G_content, C_content) -> overlap_cds

rub_only %>%
  inner_join(CDS, by = c("transcript_v" = "transcript")) %>%
  mutate(group = rep("Rubio\net al.")) %>%
  select(group, GC_content, length, G_content, C_content) -> rub_cds

wolf_only %>%
  inner_join(CDS, by = c("transcript_v" = "transcript")) %>%
  mutate(group = rep("Wolfe\net al.")) %>%
  select(group, GC_content, length, G_content, C_content) -> wolf_cds

RW_cds <- rbind(rub_cds, overlap_cds, wolf_cds)
RW_cds %>%
  mutate(group = fct_relevel(group, "Rubio\net al.", "Both", "Wolfe\net al.")) -> RW_cds

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
res.aov_cds_GC_rw <- aov(GC_content ~ group, data = RW_cds) #GC content
summary(res.aov_cds_GC_rw)
Tuk_cds_GC_RW <- TukeyHSD(res.aov_cds_GC_rw)
Tuk_cds_GC_RW
p_Rb_cds_GC <- myP(Tuk_cds_GC_RW$group[10])
p_WR_cds_GC <- myP(Tuk_cds_GC_RW$group[11])
p_Wb_cds_GC <- myP(Tuk_cds_GC_RW$group[12])

res.aov_cds_l_rw <- aov(log10(length) ~ group, data = RW_cds) #length
summary(res.aov_cds_l_rw)
Tuk_cds_l_RW <- TukeyHSD(res.aov_cds_l_rw)
Tuk_cds_l_RW
p_Rb_cds_l <- myP(Tuk_cds_l_RW$group[10])
p_WR_cds_l <- myP(Tuk_cds_l_RW$group[11])
p_Wb_cds_l <- myP(Tuk_cds_l_RW$group[12])

res.aov_cds_G_rw <- aov(G_content ~ group, data = RW_cds) #G content
summary(res.aov_cds_G_rw)
Tuk_cds_G_RW <- TukeyHSD(res.aov_cds_G_rw)
Tuk_cds_G_RW
p_Rb_cds_G <- myP(Tuk_cds_G_RW$group[10])
p_WR_cds_G <- myP(Tuk_cds_G_RW$group[11])
p_Wb_cds_G <- myP(Tuk_cds_G_RW$group[12])

res.aov_cds_C_rw <- aov(C_content ~ group, data = RW_cds) #C content
summary(res.aov_cds_C_rw)
Tuk_cds_C_RW <- TukeyHSD(res.aov_cds_C_rw)
Tuk_cds_C_RW
p_Rb_cds_C <- myP(Tuk_cds_C_RW$group[10])
p_WR_cds_C <- myP(Tuk_cds_C_RW$group[11])
p_Wb_cds_C <- myP(Tuk_cds_C_RW$group[12])


#3' UTR
##Data manipulation
overlap %>%
  inner_join(tpUTR, by = c("transcript_v" = "transcript")) %>% #Coupling 3'UTR features
  mutate(group = rep("Both")) %>%
  select(group, GC_content, length, G_content, C_content) -> overlap_tp

rub_only %>%
  inner_join(tpUTR, by = c("transcript_v" = "transcript")) %>%
  mutate(group = rep("Rubio\net al.")) %>%
  select(group, GC_content, length, G_content, C_content) -> rub_tp

wolf_only %>%
  inner_join(tpUTR, by = c("transcript_v" = "transcript")) %>%
  mutate(group = rep("Wolfe\net al.")) %>%
  select(group, GC_content, length, G_content, C_content) -> wolf_tp

RW_tp <- rbind(rub_tp, overlap_tp, wolf_tp)
RW_tp %>%
  mutate(group = fct_relevel(group, "Rubio\net al.", "Both", "Wolfe\net al.")) -> RW_tp

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
res.aov_tp_GC_rw <- aov(GC_content ~ group, data = RW_tp) #GC content
summary(res.aov_tp_GC_rw)
Tuk_tp_GC_RW <- TukeyHSD(res.aov_tp_GC_rw)
Tuk_tp_GC_RW
p_Rb_tp_GC <- myP(Tuk_tp_GC_RW$group[10])
p_WR_tp_GC <- myP(Tuk_tp_GC_RW$group[11])
p_Wb_tp_GC <- myP(Tuk_tp_GC_RW$group[12])

res.aov_tp_l_rw <- aov(log10(length) ~ group, data = RW_tp) #length
summary(res.aov_tp_l_rw)
Tuk_tp_l_RW <- TukeyHSD(res.aov_tp_l_rw)
Tuk_tp_l_RW
p_Rb_tp_l <- myP(Tuk_tp_l_RW$group[10])
p_WR_tp_l <- myP(Tuk_tp_l_RW$group[11])
p_Wb_tp_l <- myP(Tuk_tp_l_RW$group[12])

res.aov_tp_G_rw <- aov(G_content ~ group, data = RW_tp) #G content
summary(res.aov_tp_G_rw)
Tuk_tp_G_RW <- TukeyHSD(res.aov_tp_G_rw)
Tuk_tp_G_RW
p_Rb_tp_G <- myP(Tuk_tp_G_RW$group[10])
p_WR_tp_G <- myP(Tuk_tp_G_RW$group[11])
p_Wb_tp_G <- myP(Tuk_tp_G_RW$group[12])

res.aov_tp_C_rw <- aov(C_content ~ group, data = RW_tp) #C content
summary(res.aov_tp_C_rw)
Tuk_tp_C_RW <- TukeyHSD(res.aov_tp_C_rw)
Tuk_tp_C_RW
p_Rb_tp_C <- myP(Tuk_tp_C_RW$group[10])
p_WR_tp_C <- myP(Tuk_tp_C_RW$group[11])
p_Wb_tp_C <- myP(Tuk_tp_C_RW$group[12])


#Plot---
##5'UTR
##GC content
GC_RW_fp <- ggplot(RW_fp, aes(x = group, y = GC_content, fill = group)) +
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
  annotate("text", x = 1.5, y = 0.95, size=4, label = p_Rb_fp_GC)+
  annotate("text", x = 2.5, y = 0.95, size=4, label = p_Wb_fp_GC)+
  annotate("text", x = 2, y = 0.95, size=4, label = p_WR_fp_GC)

##length
length_RW_fp <- ggplot(RW_fp, aes(x = group, y = length, fill = group)) +
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
  annotate("text", x = 1.5, y = 1500, size=4, label = p_Rb_fp_l)+
  annotate("text", x = 2.5, y = 1500, size=4, label = p_Wb_fp_l)+
  annotate("text", x = 2, y = 100, size=4, label = p_WR_fp_l)

#G content
G_RW_fp <- ggplot(RW_fp, aes(x = group, y = G_content, fill = group)) +
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
  annotate("text", x = 1.5, y = 0.95, size=4, label = p_Rb_fp_G)+
  annotate("text", x = 2.5, y = 0.95, size=4, label = p_Wb_fp_G)+
  annotate("text", x = 2, y = 0.95, size=4, label = p_WR_fp_G)

#C content
C_RW_fp <- ggplot(RW_fp, aes(x = group, y = C_content, fill = group)) +
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
  annotate("text", x = 1.5, y = 0.95, size=4, label = p_Rb_fp_C)+
  annotate("text", x = 2.5, y = 0.95, size=4, label = p_Wb_fp_C)+
  annotate("text", x = 2, y = 0.95, size=4, label = p_WR_fp_C)

##CDS
##GC content
GC_RW_cds <- ggplot(RW_cds, aes(x = group, y = GC_content, fill = group)) +
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
  annotate("text", x = 1.5, y = 0.95, size=4, label = p_Rb_cds_GC)+
  annotate("text", x = 2.5, y = 0.95, size=4, label = p_Wb_cds_GC)+
  annotate("text", x = 2, y = 0.95, size=4, label = p_WR_cds_GC)

##length
length_RW_cds <- ggplot(RW_cds, aes(x = group, y = length, fill = group)) +
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
  annotate("text", x = 1.5, y = 1500, size=4, label = p_Rb_cds_l)+
  annotate("text", x = 2.5, y = 1500, size=4, label = p_Wb_cds_l)+
  annotate("text", x = 2, y = 100, size=4, label = p_WR_cds_l)

##G content
G_RW_cds <- ggplot(RW_cds, aes(x = group, y = G_content, fill = group)) +
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
  annotate("text", x = 1.5, y = 0.95, size=4, label = p_Rb_cds_G)+
  annotate("text", x = 2.5, y = 0.95, size=4, label = p_Wb_cds_G)+
  annotate("text", x = 2, y = 0.95, size=4, label = p_WR_cds_G)

##C content
C_RW_cds <- ggplot(RW_cds, aes(x = group, y = C_content, fill = group)) +
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
  annotate("text", x = 1.5, y = 0.95, size=4, label = p_Rb_cds_C)+
  annotate("text", x = 2.5, y = 0.95, size=4, label = p_Wb_cds_C)+
  annotate("text", x = 2, y = 0.95, size=4, label = p_WR_cds_C)

##3'UTR
##GC content
GC_RW_tp <- ggplot(RW_tp, aes(x = group, y = GC_content, fill = group)) +
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
  annotate("text", x = 1.5, y = 0.95, size=4, label = p_Rb_tp_GC)+
  annotate("text", x = 2.5, y = 0.95, size=4, label = p_Wb_tp_GC)+
  annotate("text", x = 2, y = 0.95, size=4, label = p_WR_tp_GC)

##length
length_RW_tp <- ggplot(RW_tp, aes(x = group, y = length, fill = group)) +
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
  annotate("text", x = 1.5, y = 1500, size=4, label = p_Rb_tp_l)+
  annotate("text", x = 2.5, y = 1500, size=4, label = p_Wb_tp_l)+
  annotate("text", x = 2, y = 1500, size=4, label = p_WR_tp_l)

##G content
G_RW_tp <- ggplot(RW_tp, aes(x = group, y = G_content, fill = group)) +
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
  annotate("text", x = 1.5, y = 0.95, size=4, label = p_Rb_tp_G)+
  annotate("text", x = 2.5, y = 0.95, size=4, label = p_Wb_tp_G)+
  annotate("text", x = 2, y = 0.95, size=4, label = p_WR_tp_G)

##C content
C_RW_tp <- ggplot(RW_tp, aes(x = group, y = C_content, fill = group)) +
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
  annotate("text", x = 1.5, y = 0.95, size=4, label = p_Rb_tp_C)+
  annotate("text", x = 2.5, y = 0.95, size=4, label = p_Wb_tp_C)+
  annotate("text", x = 2, y = 0.95, size=4, label = p_WR_tp_C)

