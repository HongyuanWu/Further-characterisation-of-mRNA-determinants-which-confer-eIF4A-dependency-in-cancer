#DENSITY PLOT_RUBIO
#This script has been created by Rocco Roberto Penna with the support of Dr. Joseph A. Waldron

#AIM: to show data distribution of data sets for each study
#OBJECTIVE: To create density plot depicting Rubio et al. 2014 data sets


#R packages---
library(dplyr)
library(ggplot2)


#Data import---
rub_dep <- read.table("Files/Rubio_dep.txt", header = T) #eIF4A-dep
rub_antidep <- read.table("Files/Rubio_antidep.txt", header = T) #eIF4A-antidep

rub_data <- rbind(rub_antidep, rub_dep) #merging data
rm(rub_dep, rub_antidep) #remove unnecessary data


#Data manipulation---
rub_data %>%
  mutate(alpha_score = rep(1),
         translation = factor(case_when(D.TE < 0 ~ "4AI-dep",
                                        D.TE > 0 ~ "4AI-antidep"),
                              levels = c("4AI-dep", "4AI-antidep"),
                              labels = c("4AI-dependent", "4AI-antidependent"),
                              ordered = T)) %>% #creating new columns
  rename(gene_name = Refseq_name) %>%
  select(gene_name, D.TE, z.score, translation, alpha_score)-> rub_plot_data #isolated necessary info for the plot

dep_mean <- mean(rub_plot_data$D.TE[rub_plot_data$translation == "4AI-dependent"]) #mean translational efficiency eIF4A-dep
antidep_mean <- mean(rub_plot_data$D.TE[rub_plot_data$translation == "4AI-antidependent"]) #mean translational efficiency eIF4A-antidep


#Plot---
##Theme
myTheme <- theme_bw()+
  theme(axis.text = element_text(size=18), 
        axis.title = element_text(size=20), 
        legend.title = element_blank(),
        legend.text = element_text(size=18))

density_rub <- ggplot(rub_plot_data, aes(x = D.TE, colour = translation, fill = translation))+
  geom_density(alpha = 0.5)+
  scale_fill_manual(values=c("#74add1", "#fdae61"))+
  scale_colour_manual(values=c("#74add1", "#fdae61"))+
  xlab("Translation Efficiency \n(log2 fold-change)") +
  ylab("Density")+
  xlim(-10, 10)+
  geom_vline(xintercept = 0, lty = 2)+
  geom_vline(xintercept = dep_mean, colour = "#74add1", linetype="dashed", size = 2)+
  geom_vline(xintercept = antidep_mean, colour = "#fdae61", linetype="dashed", size = 2)+
  geom_hline(yintercept = 1)+
  #geom_text(aes(label=ifelse(D.TE > 3,yes = gene_name,'')),hjust=0,vjust=0)+
  myTheme

