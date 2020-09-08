#VENN DIAGRAMS_RIBOSOME FOOTPRINTING
#This script has been created by Rocco Roberto Penna with the support of Dr. Joseph A. Waldron

#AIM: to detect overlaps among studies
#OBJECTIVE: to plot in Venn diagrams data obtained using ribosome footprinting method

#Studies that used Ribosome Footprinting (RF) method: Rubio et al. (2014), Wolfe et al. (2014) and Iwasaki et al. (2016)
#Iwasaki study obtained two data sets from Hippuristanol and RocA. 3 type of plots have been created: Hipp only, RocA only and all datasets

#R packages---
library(tidyverse)
library(VennDiagram)

#Data---
source("I_Venn Diagrams/Venn Diagrams_Sources.R")

#Hippuristanol only plot---
#eIF4A-dependent data
overlap_123_dep_rfH <- factor(intersect(intersect(iwa_hipp_dep_genes, wolfe_dep_genes), rub_dep_genes)) #Overlap among Iwasaki (Hipp), Wolfe and Rubio data
overlap_12_dep_rfH <- factor(intersect(iwa_hipp_dep_genes, wolfe_dep_genes)) #Overlap between Iwasaki (Hipp) and Wolfe data
overlap_13_dep_rfH <- factor(intersect(iwa_hipp_dep_genes, rub_dep_genes)) #Overlap between Iwasaki (Hipp) and Rubio data
overlap_23_dep_rfH <- factor(intersect(wolfe_dep_genes, rub_dep_genes)) #Overlap between Wolfe and Rubio data


grid.newpage() #reset "Plots" screen
venn_plot_dep_rfH <- draw.triple.venn(area1 = length(iwa_hipp_dep_genes),
                                         area2 = length(wolfe_dep_genes),
                                         area3 = length(rub_dep_genes),
                                         n123 = length(overlap_123_dep_rfH),
                                         n12 = length(overlap_12_dep_rfH),
                                         n13 = length(overlap_13_dep_rfH),
                                         n23 = length(overlap_23_dep_rfH),
                                         category = c("Iwasaki et al.\n(Hipp)", "Wolfe et al.", "Rubio et al."),
                                         lty = rep("blank",3),
                                         fill = c("light blue", "red", "yellow"),
                                         alpha = rep(0.3, 3), 
                                         euler.d = T,
                                         scaled = T,
                                         cex = 2,
                                         cat.cex = 2,
                                         overrideTriple = T)

#eIF4A-antidependent data
overlap_123_antidep_rfH <- factor(intersect(intersect(iwa_hipp_antidep_genes, wolfe_antidep_genes), rub_antidep_genes)) #Overlap among Iwasaki (Hipp), Wolfe and Rubio data
overlap_12_antidep_rfH <- factor(intersect(iwa_hipp_antidep_genes, wolfe_antidep_genes)) #Overlap between Iwasaki (Hipp) and Wolfe data
overlap_13_antidep_rfH <- factor(intersect(iwa_hipp_antidep_genes, rub_antidep_genes)) #Overlap between Iwasaki (Hipp) and Rubio data
overlap_23_antidep_rfH <- factor(intersect(wolfe_antidep_genes, rub_antidep_genes)) #Overlap between Wolfe and Rubio data

grid.newpage()
venn_plot_antidep_rfH <- draw.triple.venn(area1 = length(iwa_hipp_antidep_genes),
                                     area2 = length(wolfe_antidep_genes),
                                     area3 = length(rub_antidep_genes),
                                     n123 = length(overlap_123_antidep_rfH),
                                     n12 = length(overlap_12_antidep_rfH),
                                     n13 = length(overlap_13_antidep_rfH),
                                     n23 = length(overlap_23_antidep_rfH),
                                     category = c("Iwasaki et al.\n(Hipp)", "Wolfe et al.", "Rubio et al."),
                                     lty = rep("blank",3),
                                     fill = c("green", "orange3", "purple"),
                                     alpha = rep(0.3, 3), 
                                     euler.d = T,
                                     scaled = T,
                                     cex = 2,
                                     cat.cex = 2,
                                     overrideTriple = T)

#####################################################################################################################################

#RocA only plot---
#eIF4A-dependent data
overlap_123_dep_rfR <- factor(intersect(intersect(iwa_roca_dep_genes, wolfe_dep_genes), rub_dep_genes)) #Overlap among Iwasaki (RocA), Wolfe and Rubio data
overlap_12_dep_rfR <- factor(intersect(iwa_roca_dep_genes, wolfe_dep_genes)) #Overlap between Iwasaki (RocA) and Wolfe data
overlap_13_dep_rfR <- factor(intersect(iwa_roca_dep_genes, rub_dep_genes)) #Overlap between Iwasaki (RocA) and Rubio data
overlap_23_dep_rfR <- factor(intersect(wolfe_dep_genes, rub_dep_genes)) #Overlap between Wolfe and Rubio data


grid.newpage() #reset "Plots" screen
venn_plot_dep_rfR <- draw.triple.venn(area1 = length(iwa_roca_dep_genes),
                                      area2 = length(wolfe_dep_genes),
                                      area3 = length(rub_dep_genes),
                                      n123 = length(overlap_123_dep_rfR),
                                      n12 = length(overlap_12_dep_rfR),
                                      n13 = length(overlap_13_dep_rfR),
                                      n23 = length(overlap_23_dep_rfR),
                                      category = c("Iwasaki et al.\n(RocA)", "Wolfe et al.", "Rubio et al."),
                                      lty = rep("blank",3),
                                      fill = c("light blue", "red", "yellow"),
                                      alpha = rep(0.3, 3), 
                                      euler.d = T,
                                      scaled = T,
                                      cex = 2,
                                      cat.cex = 2,
                                      overrideTriple = T)

#eIF4A-antidependent data
overlap_123_antidep_rfR <- factor(intersect(intersect(iwa_roca_antidep_genes, wolfe_antidep_genes), rub_antidep_genes)) #Overlap among Iwasaki (RocA), Wolfe and Rubio data
overlap_12_antidep_rfR <- factor(intersect(iwa_roca_antidep_genes, wolfe_antidep_genes)) #Overlap between Iwasaki (RocA) and Wolfe data
overlap_13_antidep_rfR <- factor(intersect(iwa_roca_antidep_genes, rub_antidep_genes)) #Overlap between Iwasaki (RocA) and Rubio data
overlap_23_antidep_rfR <- factor(intersect(wolfe_antidep_genes, rub_antidep_genes)) #Overlap between Wolfe and Rubio data

grid.newpage()
venn_plot_antidep_rfR <- draw.triple.venn(area1 = length(iwa_roca_antidep_genes),
                                          area2 = length(wolfe_antidep_genes),
                                          area3 = length(rub_antidep_genes),
                                          n123 = length(overlap_123_antidep_rfR),
                                          n12 = length(overlap_12_antidep_rfR),
                                          n13 = length(overlap_13_antidep_rfR),
                                          n23 = length(overlap_23_antidep_rfR),
                                          category = c("Iwasaki et al.\n(RocA)", "Wolfe et al.", "Rubio et al."),
                                          lty = rep("blank",3),
                                          fill = c("green", "orange3", "purple"),
                                          alpha = rep(0.3, 3), 
                                          euler.d = T,
                                          scaled = T,
                                          cex = 2,
                                          cat.cex = 2,
                                          overrideTriple = T)

#####################################################################################################################################

#All data sets plot---
#eIF4A-dependent data
over_1234_d <- intersect(intersect(intersect(iwa_hipp_dep_genes, iwa_roca_dep_genes), rub_dep_genes), wolfe_dep_genes)
over_12_d <- intersect(iwa_hipp_dep_genes, iwa_roca_dep_genes)
over_13_d <- intersect(iwa_hipp_dep_genes, rub_dep_genes)
over_14_d <- intersect(iwa_hipp_dep_genes, wolfe_dep_genes)
over_23_d <- intersect(iwa_roca_dep_genes, rub_dep_genes)
over_24_d <- intersect(iwa_roca_dep_genes, wolfe_dep_genes)
over_34_d <- intersect(rub_dep_genes, wolfe_dep_genes)
over_123_d <- intersect(intersect(iwa_hipp_dep_genes, iwa_roca_dep_genes), rub_dep_genes)
over_124_d <- intersect(intersect(iwa_hipp_dep_genes, iwa_roca_dep_genes), wolfe_dep_genes)
over_134_d <- intersect(intersect(iwa_hipp_dep_genes, rub_dep_genes), wolfe_dep_genes)
over_234_d <- intersect(intersect(iwa_roca_dep_genes, rub_dep_genes), wolfe_dep_genes)

grid.newpage()
RF_dep <- draw.quad.venn(area1 = length(iwa_hipp_dep_genes),
                         area2 = length(iwa_roca_dep_genes),
                         area3 = length(rub_dep_genes),
                         area4 = length(wolfe_dep_genes),
                         n12 = length(over_12_d),
                         n13 = length(over_13_d),
                         n14 = length(over_14_d),
                         n23 = length(over_23_d),
                         n24 = length(over_24_d),
                         n34 = length(over_34_d),
                         n123 = length(over_123_d),
                         n124 = length(over_124_d),
                         n134 = length(over_134_d),
                         n234 = length(over_234_d),
                         n1234 = length(over_1234_d),
                         category = c("Iwasaki et al.\n(Hipp)", "Iwasaki et al.\n(RocA)", "Rubio et al.", "Wolfe et al."),
                         lty = rep("blank",4),
                         fill = c("light blue", "red", "yellow", "green"),
                         alpha = rep(0.3, 0.6, 1, 2), 
                         euler.d = T,
                         scaled = T,
                         cex = 2,
                         cat.cex = 2,
                         overrideTriple = T)


#eIF4A-antidependent data
over_1234_a <- intersect(intersect(intersect(iwa_hipp_antidep_genes, iwa_roca_antidep_genes), rub_antidep_genes), wolfe_antidep_genes)
over_12_a <- intersect(iwa_hipp_antidep_genes, iwa_roca_antidep_genes)
over_13_a <- intersect(iwa_hipp_antidep_genes, rub_antidep_genes)
over_14_a <- intersect(iwa_hipp_antidep_genes, wolfe_antidep_genes)
over_23_a <- intersect(iwa_roca_antidep_genes, rub_antidep_genes)
over_24_a <- intersect(iwa_roca_antidep_genes, wolfe_antidep_genes)
over_34_a <- intersect(rub_antidep_genes, wolfe_antidep_genes)
over_123_a <- intersect(intersect(iwa_hipp_antidep_genes, iwa_roca_antidep_genes), rub_antidep_genes)
over_124_a <- intersect(intersect(iwa_hipp_antidep_genes, iwa_roca_antidep_genes), wolfe_antidep_genes)
over_134_a <- intersect(intersect(iwa_hipp_antidep_genes, rub_antidep_genes), wolfe_antidep_genes)
over_234_a <- intersect(intersect(iwa_roca_antidep_genes, rub_antidep_genes), wolfe_antidep_genes)


grid.newpage()
RF_antidep <- draw.quad.venn(area1 = length(iwa_hipp_antidep_genes),
                         area2 = length(iwa_roca_antidep_genes),
                         area3 = length(rub_antidep_genes),
                         area4 = length(wolfe_antidep_genes),
                         n12 = length(over_12_a),
                         n13 = length(over_13_a),
                         n14 = length(over_14_a),
                         n23 = length(over_23_a),
                         n24 = length(over_24_a),
                         n34 = length(over_34_a),
                         n123 = length(over_123_a),
                         n124 = length(over_124_a),
                         n134 = length(over_134_a),
                         n234 = length(over_234_a),
                         n1234 = length(over_1234_a),
                         category = c("Iwasaki et al.\n(Hipp)", "Iwasaki et al.\n(RocA)", "Rubio et al.", "Wolfe et al."),
                         lty = rep("blank",4),
                         fill = c("green", "orange3", "purple", "blue"),
                         alpha = rep(0.3, 0.6, 1, 2), 
                         euler.d = T,
                         scaled = T,
                         cex = 2,
                         cat.cex = 2,
                         overrideTriple = T)

