#BOX PLOTS_5'UTR FEATURES
#This script has been created by Rocco Roberto Penna with the support of Dr. Joseph A. Waldron

#AIM: to compare mRNA features generally linked to eIF4A-dependency with eIF4A-antidependent data
#OBJECTIVE: to compare 5'UTR length, GC content, G content and C content between eIF4A-dep and antidep data


#R packages---
library(dplyr)
library(ggplot2)


#Data import---
##5'UTR features
fpUTR <- read.csv(path.expand("Files/gencode_v27_pc_transcripts_fpUTRs_composition.csv"), header = T)
fpUTR %>%
  select(transcript, GC_content, length, G_content, C_content) -> fpUTR #Isolating 5'UTR length, GC, G & C contents of all transcripts

##eIF4A dependent data
source("III_Box and Violin plots/Box Plots_Sources_Dependent.R") #Importing eIF4A-dep data

chan_final_d %>%
  inner_join(fpUTR, by = c("transcript_v" = "transcript")) %>% #Coupling 5'UTR features according to study
  mutate(group = rep("Chan et al."),
         translation = rep("dependent")) %>% #Adding unique label
  select(group, translation, GC_content, length, G_content, C_content) -> chan_fp_dep #Isolating transcripts. 662 transcripts

wal_final_d %>%
  inner_join(fpUTR, by = c("transcript_v" = "transcript")) %>%
  mutate(group = rep("Waldron et al."),
         translation = rep("dependent"))%>%
  select(group,translation, GC_content, length, G_content, C_content) -> wal_fp_dep #420 transcripts

mod_final_d %>%
  inner_join(fpUTR, by = c("transcript_v" = "transcript")) %>%
  mutate(group = rep("Modelska et al."),
         translation = rep("dependent")) %>%
  select(group, translation, GC_content, length, G_content, C_content) -> mod_fp_dep #84 transcripts

wolf_final_d %>%
  inner_join(fpUTR, by = c("transcript_v" = "transcript")) %>%
  mutate(group = rep("Wolfe et al."),
         translation = rep("dependent")) %>%
  select(group, translation, GC_content, length, G_content, C_content) -> wolf_fp_dep #168 transcripts

rub_final_d %>%
  inner_join(fpUTR, by = c("transcript_v" = "transcript")) %>%
  mutate(group = rep("Rubio et al."),
         translation = rep("dependent")) %>%
  select(group, translation, GC_content, length, G_content, C_content) -> rub_fp_dep #172 transcripts

iwa_hipp_final_d %>%
  inner_join(fpUTR, by = c("transcript_v" = "transcript")) %>%
  mutate(group = rep("Iwasaki et al.\n(Hippuristanol)"),
         translation = rep("dependent")) %>%
  select(group, translation, GC_content, length, G_content, C_content) -> iwa_hipp_fp_dep #258 transcripts

iwa_roca_final_d %>%
  inner_join(fpUTR, by = c("transcript_v" = "transcript")) %>%
  mutate(group = rep("Iwasaki et al.\n(RocA)"),
         translation = rep("dependent")) %>%
  select(group, translation, GC_content, length, G_content, C_content) -> iwa_roca_fp_dep #183 transcripts

##eIF4A antidependent data
source("III_Box and Violin plots/Box Plots_Sources_Antidependent.R") #Importing eIF4A-antidep data

chan_final_a %>%
  inner_join(fpUTR, by = c("transcript_v" = "transcript")) %>%
  mutate(group = rep("Chan et al."),
         translation = rep("antidependent")) %>%
  select(group, translation, GC_content, length, G_content, C_content) -> chan_fp_anti #617 transcripts

wal_final_a %>%
  inner_join(fpUTR, by = c("transcript_v" = "transcript")) %>%
  mutate(group = rep("Waldron et al."),
         translation = rep("antidependent")) %>%
  select(group, translation, GC_content, length, G_content, C_content) -> wal_fp_anti #31 transcripts

mod_final_a %>%
  inner_join(fpUTR, by = c("transcript_v" = "transcript")) %>%
  mutate(group = rep("Modelska et al."),
         translation = rep("antidependent")) %>%
  select(group, translation, GC_content, length, G_content, C_content) -> mod_fp_anti #31 transcripts

wolf_final_a %>%
  inner_join(fpUTR, by = c("transcript_v" = "transcript")) %>%
  mutate(group = rep("Wolfe et al."),
         translation = rep("antidependent")) %>%
  select(group, translation, GC_content, length, G_content, C_content) -> wolf_fp_anti #102 transcripts

rub_final_a %>%
  inner_join(fpUTR, by = c("transcript_v" = "transcript")) %>%
  mutate(group = rep("Rubio et al."),
         translation = rep("antidependent")) %>%
  select(group, translation, GC_content, length, G_content, C_content) -> rub_fp_anti #86 transcripts

iwa_hipp_final_a %>%
  inner_join(fpUTR, by = c("transcript_v" = "transcript")) %>%
  mutate(group = rep("Iwasaki et al.\n(Hippuristanol)"),
         translation = rep("antidependent")) %>%
  select(group, translation, GC_content, length, G_content, C_content) -> iwa_hipp_fp_anti #249 transcripts

iwa_roca_final_a %>%
  inner_join(fpUTR, by = c("transcript_v" = "transcript")) %>%
  mutate(group = rep("Iwasaki et al.\n(RocA)"),
         translation = rep("antidependent")) %>%
  select(group, translation, GC_content, length, G_content, C_content) -> iwa_roca_fp_anti #176 transcripts

##Background
most_abundant %>%
  inner_join(fpUTR, by = c("transcript_v" = "transcript")) %>%
  mutate(group = rep("Backgroud"),
         translation = rep("both")) %>%
  select(group, translation, GC_content, length, G_content, C_content) -> back #12406 transcripts


#Merging data---
total_data_fp <- rbind(chan_fp_dep, chan_fp_anti,
                       wal_fp_dep, wal_fp_anti,
                       mod_fp_dep, mod_fp_anti,
                       wolf_fp_dep, wolf_fp_anti,
                       rub_fp_dep, rub_fp_anti,
                       iwa_hipp_fp_dep, iwa_hipp_fp_anti,
                       iwa_roca_fp_dep, iwa_roca_fp_anti,
                       back) #Merging all data in a single data frame

chan <- rbind(chan_fp_dep, chan_fp_anti) #Merging for each study
wal <- rbind(wal_fp_dep, wal_fp_anti)
mod <- rbind(mod_fp_dep, mod_fp_anti)
wolf <- rbind(wolf_fp_dep, wolf_fp_anti)
rub <- rbind(rub_fp_dep, rub_fp_anti)
iwa_hipp <- rbind(iwa_hipp_fp_dep, iwa_hipp_fp_anti)
iwa_roc <- rbind(iwa_roca_fp_dep, iwa_roca_fp_anti)


#Removing unneccessary data---
rm(chan_final_d, chan_final_a, chan_fp_dep, chan_fp_anti,
   wal_final_d, wal_final_a, wal_fp_dep, wal_fp_anti,
   mod_final_d, mod_final_a, mod_fp_dep, mod_fp_anti,
   wolf_final_d, wolf_final_a, wolf_fp_dep, wolf_fp_anti,
   rub_final_d, rub_final_a, rub_fp_dep, rub_fp_anti,
   iwa_hipp_final_d, iwa_hipp_final_a, iwa_hipp_fp_dep, iwa_hipp_fp_anti,
   iwa_roca_final_d, iwa_roca_final_a, iwa_roca_fp_dep, iwa_roca_fp_anti,
   back, most_abundant, fpUTR)


#Statistical analysis (t test)---
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

#Two-sided Welch two sample t-test
##Chan data
t_chan_GC <- t.test(GC_content ~ translation, data = chan) #t-test. GC content
chan_label_GC <- myP(t_chan_GC$p.value) #label

t_chan_l <- t.test(log10(length) ~ translation, data = chan) #length
chan_label_l <- myP(t_chan_l$p.value)

t_chan_G <- t.test(G_content ~ translation, data = chan) #G content
chan_label_G <- myP(t_chan_G$p.value)

t_chan_C <- t.test(C_content ~ translation, data = chan) #C content
chan_label_C <- myP(t_chan_C$p.value)

##Waldron data
t_wal_GC <- t.test(GC_content ~ translation, data = wal) #GC content
wal_label_GC <- myP(t_wal_GC$p.value)

t_wal_l <- t.test(log10(length) ~ translation, data = wal) #length
wal_label_l <- myP(t_wal_l$p.value)

t_wal_G <- t.test(G_content ~ translation, data = wal) #G content
wal_label_G <- myP(t_wal_G$p.value)

t_wal_C <- t.test(C_content ~ translation, data = wal) #C content
wal_label_C <- myP(t_wal_C$p.value)

##Modelska data
t_mod_GC <- t.test(GC_content ~ translation, data = mod) #GC content
mod_label_GC <- myP(t_mod_GC$p.value)

t_mod_l <- t.test(log10(length) ~ translation, data = mod) #length
mod_label_l <- myP(t_mod_l$p.value)

t_mod_G <- t.test(G_content ~ translation, data = mod) #G content
mod_label_G <- myP(t_mod_G$p.value)

t_mod_C <- t.test(C_content ~ translation, data = mod) #C content
mod_label_C <- myP(t_mod_C$p.value)

##Rubio data
t_rub_GC <- t.test(GC_content ~ translation, data = rub) #GC content
rub_label_GC <- myP(t_rub_GC$p.value)

t_rub_l <- t.test(log10(length) ~ translation, data = rub) #length
rub_label_l <- myP(t_rub_l$p.value)

t_rub_G <- t.test(G_content ~ translation, data = rub) #G content
rub_label_G <- myP(t_rub_G$p.value)

t_rub_C <- t.test(C_content ~ translation, data = rub) #C content
rub_label_C <- myP(t_rub_C$p.value)

##Wolfe data
t_wolf_GC <- t.test(GC_content ~ translation, data = wolf) #GC content
wolf_label_GC <- myP(t_wolf_GC$p.value)

t_wolf_l <- t.test(log10(length) ~ translation, data = wolf) #length
wolf_label_l <- myP(t_wolf_l$p.value)

t_wolf_G <- t.test(G_content ~ translation, data = wolf) #G content
wolf_label_G <- myP(t_wolf_G$p.value)

t_wolf_C <- t.test(C_content ~ translation, data = wolf) #C content
wolf_label_C <- myP(t_wolf_C$p.value)

##Iwasaki data (Hippuristanol)
t_iwaH_GC <- t.test(GC_content ~ translation, data = iwa_hipp) #GC content
iwaH_label_GC <- myP(t_iwaH_GC$p.value)

t_iwaH_l <- t.test(log10(length) ~ translation, data = iwa_hipp) #length
iwaH_label_l <- myP(t_iwaH_l$p.value)

t_iwaH_G <- t.test(G_content ~ translation, data = iwa_hipp) #G content
iwaH_label_G <- myP(t_iwaH_G$p.value)

t_iwaH_C <- t.test(C_content ~ translation, data = iwa_hipp) #C content
iwaH_label_C <- myP(t_iwaH_C$p.value)

##Iwasaki data (RocA)
t_iwaR_GC <- t.test(GC_content ~ translation, data = iwa_roc) #GC content
iwaR_label_GC <- myP(t_iwaR_GC$p.value)

t_iwaR_l <- t.test(log10(length) ~ translation, data = iwa_roc) #length
iwaR_label_l <- myP(t_iwaR_l$p.value)

t_iwaR_G <- t.test(G_content ~ translation, data = iwa_roc) #G content
iwaR_label_G <- myP(t_iwaR_G$p.value)

t_iwaR_C <- t.test(C_content ~ translation, data = iwa_roc) #C content
iwaR_label_C <- myP(t_iwaR_C$p.value)


#Plot---
##GC content
Split_GC <- ggplot(total_data_fp, aes(x = group, y = GC_content, fill = translation))+
  geom_boxplot()+
  ylim(c(0.0, 1.1)) +
  ylab("5'UTR GC content (%)") +
  scale_fill_discrete(name = "Group", labels = c("eIF4A-antidependent", "Background", "eIF4A-dependent"))+
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14, face = "italic", hjust = 0.5),
        panel.grid.minor = element_blank())+
  annotate("text", x = 2, y = 1.07, size=4, label = chan_label_GC)+
  annotate("text", x = 3, y = 1.07, size=4, label = iwaH_label_GC)+
  annotate("text", x = 4, y = 1.07, size=4, label = iwaR_label_GC)+
  annotate("text", x = 5, y = 1.07, size=4, label = mod_label_GC)+
  annotate("text", x = 6, y = 1.07, size=4, label = rub_label_GC)+
  annotate("text", x = 7, y = 1.07, size=4, label = wal_label_GC)+
  annotate("text", x = 8, y = 1.07, size=4, label = wolf_label_GC)
  
##Length
Split_length <- ggplot(total_data_fp, aes(x = group, y = length, fill = translation))+
  geom_boxplot()+
  ylab("5'UTR length (nt, log10)") +
  scale_y_log10(breaks=c(10,100,1000),limits=c(10, 1500))+
  scale_fill_discrete(name = "Group", labels = c("eIF4A-antidependent", "Background", "eIF4A-dependent"))+
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14, face = "italic", hjust = 0.5),
        panel.grid.minor = element_blank())+
  annotate("text", x = 2, y = 1500, size=4, label = chan_label_l)+
  annotate("text", x = 3, y = 1500, size=4, label = iwaH_label_l)+
  annotate("text", x = 4, y = 1500, size=4, label = iwaR_label_l)+
  annotate("text", x = 5, y = 1500, size=4, label = mod_label_l)+
  annotate("text", x = 6, y = 1500, size=4, label = rub_label_l)+
  annotate("text", x = 7, y = 1500, size=4, label = wal_label_l)+
  annotate("text", x = 8, y = 1500, size=4, label = wolf_label_l)

##G content
Split_G <- ggplot(total_data_fp, aes(x = group, y = G_content, fill = translation)) +
  geom_boxplot() +
  ylim(c(0.00, 1.1)) +
  ylab("5'UTR G content (%)") +
  scale_fill_discrete(name = "Group", labels = c("eIF4A-antidependent", "Background", "eIF4A-dependent"))+
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        panel.grid.minor = element_blank())+
  annotate("text", x = 2, y = 1.07, size=4, label = chan_label_G)+
  annotate("text", x = 3, y = 1.07, size=4, label = iwaH_label_G)+
  annotate("text", x = 4, y = 1.07, size=4, label = iwaR_label_G)+
  annotate("text", x = 5, y = 1.07, size=4, label = mod_label_G)+
  annotate("text", x = 6, y = 1.07, size=4, label = rub_label_G)+
  annotate("text", x = 7, y = 1.07, size=4, label = wal_label_G)+
  annotate("text", x = 8, y = 1.07, size=4, label = wolf_label_G)

##C content
Split_C <- ggplot(total_data_fp, aes(x = group, y = C_content, fill = translation)) +
  geom_boxplot()+
  ylim(c(0.00, 1.1)) +
  ylab("5'UTR C content (%)") +
  scale_fill_discrete(name = "Group", labels = c("eIF4A-antidependent", "Background", "eIF4A-dependent"))+
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        panel.grid.minor = element_blank())+
  annotate("text", x = 2, y = 1.07, size=4, label = chan_label_C)+
  annotate("text", x = 3, y = 1.07, size=4, label = iwaH_label_C)+
  annotate("text", x = 4, y = 1.07, size=4, label = iwaR_label_C)+
  annotate("text", x = 5, y = 1.07, size=4, label = mod_label_C)+
  annotate("text", x = 6, y = 1.07, size=4, label = rub_label_C)+
  annotate("text", x = 7, y = 1.07, size=4, label = wal_label_C)+
  annotate("text", x = 8, y = 1.07, size=4, label = wolf_label_C)
