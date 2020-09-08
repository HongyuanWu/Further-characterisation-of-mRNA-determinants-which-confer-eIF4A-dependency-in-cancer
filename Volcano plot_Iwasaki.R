#VOLCANO PLOT_IWASAKI
#This script has been created by Rocco Roberto Penna with the support of Dr. Joseph A. Waldron

#AIM: to show data distribution of data sets for each study
#OBJECTIVE: To create volcano plot depicting Iwasaki et al. 2016 data sets


#R packages---
library(dplyr)
library(readxl)
library(ggplot2)


#Data import---
##RocA
iwa_dep_roca <- read_xlsx("Files/Iwasaki_et al._2016_RocA.xlsx", sheet = 1, col_names = T) #eIF4A-dep data
iwa_anti_roca <- read_xlsx("Files/Iwasaki_et al._2016_RocA.xlsx", sheet = 2, col_names = T) #eIF4A-antidep data

iwa_data_roca <- rbind(iwa_dep_roca, iwa_anti_roca) #merging data
rm(iwa_dep_roca, iwa_anti_roca) #removing unneccessary data

##Hippuristanol
iwa_dep_hipp <- read_xlsx("Files/Iwasaki_et al._2016_Hipp.xlsx", sheet = 1, col_names = T)
iwa_anti_hipp <- read_xlsx("Files/Iwasaki_et al._2016_Hipp.xlsx", sheet = 2, col_names = T)

iwa_data_hipp <- rbind(iwa_dep_hipp, iwa_anti_hipp) 
rm(iwa_dep_hipp, iwa_anti_hipp)


#Data manipulation---
##RocA
iwa_data_roca %>%
  rename(logFC = `Translation fold change to mean [log2]`, adjp = `q value`, gene_name = Gene) %>% #renaming column names
  mutate(translation = factor(case_when(adjp < 0.01 & logFC < 0 ~ "4AI-dep",
                                        adjp < 0.01 & logFC > 0 ~ "4AI-antidep"),
                              levels = c("4AI-dep", "4AI-antidep"),
                              labels = c("4AI-dependent", "4AI-antidependent"),
                              ordered = T),
         alpha_score = rep(1),
         translation_drug = factor(case_when(adjp < 0.01 & logFC < 0 ~ "4AI-depR",
                                             adjp < 0.01 & logFC > 0 ~ "4AI-antidepR"),
                                   levels = c("4AI-depR", "4AI-antidepR"),
                                   labels = c("4AI-dependent \n(RocA)", "4AI-antidependent \n(RocA)"),
                                   ordered = T) ) %>% #creating new columns
  select(logFC, adjp, translation, alpha_score, gene_name, translation_drug) -> iwaR_plot_data #isolating necessary data

##Hippuristanol
iwa_data_hipp %>%
  rename(logFC = `Translation fold change to mean [log2]`, adjp = `q value`, gene_name = Gene) %>%
  mutate(translation = factor(case_when(adjp < 0.01 & logFC < 0 ~ "4AI-dep",
                                        adjp < 0.01 & logFC > 0 ~ "4AI-antidep"),
                              levels = c("4AI-dep", "4AI-antidep"),
                              labels = c("4AI-dependent", "4AI-antidependent"),
                              ordered = T),
         alpha_score = rep(1),
         translation_drug = factor(case_when(adjp < 0.01 & logFC < 0 ~ "4AI-depH",
                                             adjp < 0.01 & logFC > 0 ~ "4AI-antidepH"),
                                   levels = c("4AI-depH", "4AI-antidepH"),
                                   labels = c("4AI-dependent \n(Hippuristanol)", "4AI-antidependent \n(Hippuristanol)"),
                                   ordered = T)) %>%
  select(logFC, adjp, translation, alpha_score, gene_name, translation_drug) -> iwaH_plot_data


#Plot---
##Theme
myTheme <- theme_bw()+
  theme(axis.text = element_text(size=18), 
        axis.title = element_text(size=18), 
        legend.title = element_blank(),
        legend.text = element_text(size=18))

scatter_iwa_R <- ggplot(iwaR_plot_data, aes(x = logFC, y =-log10(adjp), colour = translation, alpha = alpha_score))+
  geom_point(size = 1)+
  scale_colour_manual(values=c("#74add1", "#fdae61"))+
  scale_alpha(guide = F)+
  xlab("Translation Efficiency \n(log2 fold-change)")+
  ylab("Adjusted p value \n(-log10)")+
  xlim(c(-3,3))+
  ylim(c(0, 30))+
  geom_vline(xintercept = 0, lty = 2)+
  geom_hline(yintercept = -log10(0.01), lty = 2)+
  myTheme

scatter_iwa_H <- ggplot(iwaH_plot_data, aes(x = logFC, y = -log10(adjp), colour = translation, alpha = alpha_score))+
  geom_point(size = 1)+
  scale_colour_manual(values=c("#74add1", "#fdae61"))+
  scale_alpha(guide = F)+
  xlab("Translation Efficiency \n(log2 fold-change)")+
  ylab("Adjusted p value \n(-log10)")+
  xlim(c(-3,3))+
  ylim(c(0, 30))+
  geom_vline(xintercept = 0, lty = 2)+
  geom_hline(yintercept = -log10(0.01), lty = 2)+
  myTheme

#Plot tot data (merged)---
iwa_tot <- rbind(iwaH_plot_data, iwaR_plot_data)
iw_part<- inner_join(iwaH_plot_data, iwaR_plot_data, by = "gene_name")
  
scatter_iwa_Tot <- ggplot(iwa_tot, aes(x = logFC, y = -log10(adjp), colour = translation_drug, alpha = alpha_score))+
  geom_point(size = 1)+
  scale_alpha(guide = F)+
  scale_colour_manual(values = c("red", "yellow", "olivedrab3", "blue"))+
  xlab("Translation Efficiency \n(log2 fold-change)")+
  ylab("Adjusted p value \n(-log10)")+
  xlim(c(-3,3))+
  ylim(c(0, 30))+
  geom_vline(xintercept = 0, lty = 2)+
  geom_hline(yintercept = -log10(0.01), lty = 2)+
  theme_bw()+
  myTheme

