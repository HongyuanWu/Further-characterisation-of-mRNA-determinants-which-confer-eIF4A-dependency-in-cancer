#VOLCANO PLOT_WOLFE
#This script has been created by Rocco Roberto Penna with the support of Dr. Joseph A. Waldron

#AIM: to show data distribution of data sets for each study
#OBJECTIVE: To create volcano plot depicting Wolfe et al. 2014 data sets


#R packages---
library(dplyr)
library(ggplot2)
library(readxl)


#Data import---
wolf_TE_down <- read_excel("Files/Wolfe et al._2014.xlsx", range = "A5:D286") #eIF4A-dep data
wolf_TE_up <- read_excel("Files/Wolfe et al._2014.xlsx", range = "A290:D480") #eIF4A-antidep data
wolf_TE_ind <- read_excel("Files/Wolfe et al._2014.xlsx", range = "A484:D5069") #eIF4A-indep data

wolf_TE_ind %>%
  rename(`Log2(Translational Efficiency)` = `log2(Translational Efficiency)`) -> wolf_TE_ind #changing name column

wolf_data <- rbind(wolf_TE_down, wolf_TE_up, wolf_TE_ind) #merging data
rm(wolf_TE_down, wolf_TE_ind, wolf_TE_up) #removing unneccesary data


#Data manipulation---
wolf_data %>%
  rename(gene_name = `Gene name`, logFC = `Log2(Translational Efficiency)`, p_value = `Translational Efficiency (p-value)`) %>%
  mutate(alpha_score = rep(1),
         translation = factor(case_when(p_value < 0.03 & logFC < 0 ~ "4AI-dep",
                                        p_value < 0.03 & logFC > 0 ~ "4AI-antidep",
                                        p_value >= 0.03 ~ "4AI-indep"),
                              levels = c("4AI-dep", "4AI-indep", "4AI-antidep"),
                              labels = c("4AI-dependent", "4AI-independent", "4AI-antidependent"),
                              ordered = T)) %>% #creating new columns
  na.omit()%>% #omission NA data
  select(gene_name, logFC, p_value, alpha_score, translation) -> wolfe_plot_data #selecting necessary data


#Plot---
myTheme <- theme_bw()+
  theme(axis.text = element_text(size=18), 
        axis.title = element_text(size=20), 
        legend.title = element_blank(),
        legend.text = element_text(size=18))

scatter_wolfe <- ggplot(wolfe_plot_data, aes(x = logFC, y = -log10(p_value), colour = translation, alpha = alpha_score))+
  geom_point(size = 1)+
  scale_colour_manual(values=c("#74add1", "#fdae61", "#A468D5"))+
  scale_alpha(guide = F)+
  geom_hline(yintercept = -log10(0.03), lty=2)+
  geom_vline(xintercept = 0, lty=2)+
  xlim(c(-2, 2))+
  ylim(0, 6)+
  xlab("Translation Efficiency \n(log2 fold-change)")+
  ylab("P-value \n(-log10)")+
  myTheme
  
