#VENN DIAGRAMS_POLYSOME PROFILING
#This script has been created by Rocco Roberto Penna with the support of Dr. Joseph A. Waldron

#AIM: to detect overlaps among studies
#OBJECTIVE: to plot in Venn diagrams data obtained using polysome profiling method

#Studies that used Polysome Profiling (PP) method: Waldron et al. (2019), Modelska et al. (2015) and Chan et al. (2019)

#R packages---
library(tidyverse)
library(VennDiagram)

#Data---
source("I_Venn Diagrams/Venn Diagrams_Sources.R")

#eIF4A-dependent data plot---
overlap_123_dep_pp <- factor(intersect(intersect(wal_dep_genes, chan_dep_genes), mod_dep_genes)) #Overlap among Waldron, Chan and Modelska data
overlap_12_dep_pp <- factor(intersect(wal_dep_genes, chan_dep_genes)) #Overlap between Waldron and Chan data
overlap_13_dep_pp <- factor(intersect(wal_dep_genes, mod_dep_genes)) #Overlap between Waldron and Modelska data
overlap_23_dep_pp <- factor(intersect(chan_dep_genes, mod_dep_genes)) #Overlap between Chan and Modelska data

grid.newpage() #reset "Plots" screen
venn_plot_dep_pp <- draw.triple.venn(area1 = length(wal_dep_genes),
                              area2 = length(chan_dep_genes),
                              area3 = length(mod_dep_genes),
                              n123 = length(overlap_123_dep_pp),
                              n12 = length(overlap_12_dep_pp),
                              n13 = length(overlap_13_dep_pp),
                              n23 = length(overlap_23_dep_pp),
                              category = c("Waldron et al.", "Chan et al.", "Modelska et al."),
                              lty = rep("blank", 3),
                              fill = c("light blue", "red", "yellow"),
                              alpha = rep(0.3, 3), 
                              euler.d = T,
                              scaled = T,
                              cex = 2,
                              cat.cex = 2,
                              overrideTriple = T)


#eIF4A-antidependent data plot---
overlap_123_antidep_pp <- intersect(intersect(wal_antidep_genes, chan_antidep_genes), mod_antidep_genes) #Overlap among Waldron, Chan and Modelska data
overlap_12_antidep_pp <- factor(intersect(wal_antidep_genes, chan_antidep_genes)) #Overlap between Waldron and Chan data
overlap_13_antidep_pp <- factor(intersect(wal_antidep_genes, mod_antidep_genes)) #Overlap between Waldron and Modelska data
overlap_23_antidep_pp <- factor(intersect(chan_antidep_genes, mod_antidep_genes)) #Overlap between Chan and Modelska data

grid.newpage()
venn_plot_antidep_pp <- draw.triple.venn(area1 = length(wal_antidep_genes),
                                     area2 = length(chan_antidep_genes),
                                     area3 = length(mod_antidep_genes),
                                     n123 = length(overlap_123_antidep_pp),
                                     n12 = length(overlap_12_antidep_pp),
                                     n13 = length(overlap_13_antidep_pp),
                                     n23 = length(overlap_23_antidep_pp),
                                     category = c("Waldron et al.", "Chan et al.", "Modelska et al."),
                                     lty = rep("blank",3),
                                     fill = c("green", "orange3", "purple"),
                                     alpha = rep(0.3, 3), 
                                     euler.d = T,
                                     scaled = T,
                                     cex = 2,
                                     cat.cex = 2,
                                     overrideTriple = T,
                                     cat.pos = c(0, 2, 0))
