#SCATTER PLOT_MODELSKA
#This script has been created by Rocco Roberto Penna with the support of Dr. Joseph A. Waldron

#AIM: to show data distribution of data sets for each study
#OBJECTIVE: To create scatter plot depicting Modelska et al. 2015 data sets


#R packages---
library(dplyr)
library(ggplot2)


#Data import---
mod_data_dep <- read.table("Files/Mod_4A1-dep.txt", header = T) #eIF4A-dep data
mod_data_antidep <- read.table("Files/Mod_4A1-indep.txt", header = T) #eIF4A-antidep data

mod_data <- rbind(mod_data_dep, mod_data_antidep) #merging data
rm(mod_data_antidep, mod_data_dep) #removing unnecessary data


#Data manipulation---
##Set threshold
pos_change_mod <- 0.2

mod_data %>%
  filter(Posterior_probability >= pos_change_mod) %>% #filtering for significance
  mutate(alpha_score = rep(1),
         DOD = eta1_2 - eta1_1,
         translation = factor(case_when(DOD < 0 ~ "4AI-dep",
                                        DOD > 0 ~ "4AI-antidep"),
                              levels = c("4AI-dep", "4AI-antidep"),
                              labels = c("4AI-dependent", "4AI-antidependent"), order = T)) %>% #creating new columns
  rename(transcript_id = Transcript_ID, polys_logFC = eta1_2, subs_logFC = eta1_1) %>% #renaming columns 
  select(Gene_name, transcript_id, polys_logFC, subs_logFC, translation, alpha_score, DOD, Posterior_probability) -> mod_plot_data #selecting necessary data
  

#Plot---
##Theme
myTheme <- theme_bw()+
  theme(axis.text = element_text(size=18), 
        axis.title = element_text(size=20), 
        legend.title = element_blank(),
        legend.text = element_text(size=18),
        text = element_text(color = "black"))


scatter_mod <- ggplot(mod_plot_data, aes(x = subs_logFC, y = polys_logFC, colour = translation, alpha = alpha_score))+
  geom_point(size = 1)+
  scale_colour_manual(values= c("#74add1", "#fdae61"))+
  scale_alpha(guide = F)+
  geom_abline(lty=2)+
  geom_hline(yintercept = 0, lty=2)+
  geom_vline(xintercept = 0, lty=2)+
  xlim(c(-1.5, 1.5))+
  ylim(c(-1.5, 1.5))+
  xlab("Sub-polysomes log-fold change")+
  ylab("Polysomes log-fold change")+
  myTheme

