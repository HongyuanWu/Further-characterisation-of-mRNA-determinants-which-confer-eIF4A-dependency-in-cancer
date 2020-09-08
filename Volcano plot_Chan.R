#VOLCANO PLOT_CHAN
#This script has been created by Rocco Roberto Penna with the support of Dr. Joseph A. Waldron

#AIM: to show data distribution of data sets for each study
#OBJECTIVE: To create volcano plot depicting Chan et al. 2019 data sets


#R packages---
library(dplyr)
library(readxl)
library(ggplot2)


#Data import---
chan_data_dep <- read_xlsx("Files/Chan_et al._2019.xlsx", sheet = 2, col_names = T) #eIF4A-dep data
chan_data_antidep <- read_xlsx("Files/Chan_et al._2019.xlsx", sheet = 1, col_names = T) #eIF4A-antidep data

chan_data <- rbind(chan_data_dep, chan_data_antidep) #merging data


#Data manipulation---
chan_data %>%
  rename(gene_name = "Gene ID", log2FC = apvEff, adjpv = apvRvmPAdj) %>% #renaming column names
  mutate(translation = factor(case_when(log2FC < 0 ~ "4AI-dep",
                                        log2FC > 0 ~"4AI-antidep"),
                              levels = c("4AI-dep", "4AI-antidep"),
                              labels = c("4AI-dependent", "4AI-antidependent"),
                              ordered = T),
         alpha_score = rep(1)) %>% #creating new columns
  select(gene_name, log2FC, adjpv, translation, alpha_score) -> chan_plot_data #isolating necessary data


#Plot---
##Theme
myTheme <- theme_bw()+
  theme(axis.text = element_text(size=18), 
        axis.title = element_text(size=20), 
        legend.title = element_blank(),
        legend.text = element_text(size=18))

scatter_chan <- ggplot(chan_plot_data, aes(x = log2FC, y= -log10(adjpv), colour = translation, alpha = alpha_score))+
  geom_point(size = 1)+
  scale_colour_manual(values=c("#fdae61", "#74add1"))+
  scale_alpha(guide = F)+
  xlab("Translation Efficiency \n (log2 fold-change)")+
  ylab("Adjusted p-value \n (-log10)")+
  xlim(c(-3,3))+
  geom_hline(yintercept = -log10(0.25),  lty = 2)+
  myTheme

