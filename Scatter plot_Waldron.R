#SCATTER PLOT_WALDRON
#This script has been created by Rocco Roberto Penna with the support of Dr. Joseph A. Waldron

#AIM: to show data distribution of data sets for each study
#OBJECTIVE: To create SCATTER plot depicting Chan et al. 2019 data sets


#R packages---
library(dplyr)
library(ggplot2)
library(tidyverse)


#Data import---
wal_data <- read.csv("Files/Wal_penn-DOD-gene.mmdiff90.csv", header = T) #all transcripts
hom <- read_tsv("Files/homo_sapiens_core_37_75_genes_canonical_transcript.tsv", col_names =T) #homo sapiens transcripts


#Data manipulation---
hom %>%
  filter(biotype == "protein_coding") %>% #filtering protein coding transcripts
  mutate(external_name = factor(toupper(external_name)))%>% #converting to upper case
  select(external_name, gene_stable_id, canonical_transcript_stable_id) -> homo_sapiens #isolating data

#Set thresholds
positive_change <- 0.25
no_change <- 0.02

wal_data %>%
  mutate(DOD = eta1_1 - eta1_2,
         translation = factor(case_when(posterior_probability > positive_change & DOD < 0 ~ "4A-dep",
                                        posterior_probability > positive_change & DOD > 0 ~ "4A-antidep",
                                        posterior_probability < no_change ~ "4A-indep",
                                        posterior_probability <= positive_change & posterior_probability >= no_change ~ "NA"),
                              levels = c("4A-dep", "4A-indep", "4A-antidep", "NA"),
                              labels = c("4A-dependent", "4A-independent", "4A-antidependent", "not assigned"), order = T),
         alpha_score = case_when(posterior_probability > positive_change & DOD < 0 ~ 1,
                                 posterior_probability > positive_change & DOD > 0 ~ 1,
                                 posterior_probability < no_change ~ 1,
                                 posterior_probability <= positive_change & posterior_probability >= no_change ~ 0.5)) %>% #creating new columns
  rename(polys_logFC = eta1_1,
         subs_logFC = eta1_2) %>% #renaming column names
  inner_join(homo_sapiens, by = c("external_gene_name"="external_name")) %>% #filtering for protein coding transcripts
  select(polys_logFC, subs_logFC, translation, DOD, feature_id, alpha_score)-> wal_plot_data #selecting necessary data


#Scatter plot---
##Theme
myTheme <- theme_bw()+
  theme(axis.text = element_text(size=18), 
        axis.title = element_text(size=20), 
        legend.title = element_blank(),
        legend.text = element_text(size=18))

scatter_wal <- ggplot(data = wal_plot_data, aes(x = subs_logFC, y = polys_logFC, colour = translation, alpha = alpha_score))+
  geom_point(size = 1)+
  scale_colour_manual(values=c("#74add1", "#fdae61", "#A468D5", "grey"))+
  scale_alpha(guide = F)+
  geom_abline(lty=2)+
  geom_hline(yintercept = 0, lty=2)+
  geom_vline(xintercept = 0, lty=2)+
  xlim(c(-1.5, 1.5))+
  ylim(c(-1.5, 1.5))+
  xlab("Sub-polysomes log-fold change")+
  ylab("Polysomes log-fold change")+
  myTheme

