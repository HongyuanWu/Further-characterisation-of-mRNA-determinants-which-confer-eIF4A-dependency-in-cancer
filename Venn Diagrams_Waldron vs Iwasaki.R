#VENN DIAGRAMS_SOURCES
#This script has been created by Rocco Roberto Penna with the support of Dr. Joseph A. Waldron

#AIM: to confirm that the use of same cell line and technique lead to better overlaps
#OBJECTIVE: to plot in Venn diagrams data obtained by Waldron et al. (2019) and Iwasaki et al. (2016)

#R packages---
library(tidyverse)
library(VennDiagram)

#Data---
source("I_Venn Diagrams/Venn Diagrams_Sources.R")

#eIF4A-dependent data plot---
overlap_123_d <- factor(intersect(intersect(iwa_hipp_dep_genes, iwa_roca_dep_genes), wal_dep_genes))
overlap_12_d <- factor(intersect(iwa_hipp_dep_genes, iwa_roca_dep_genes))
overlap_13_d <- factor(intersect(iwa_hipp_dep_genes, wal_dep_genes))
overlap_23_d <- factor(intersect(iwa_roca_dep_genes, wal_dep_genes))

grid.newpage()
venn_plot_dep <- draw.triple.venn(area1 = length(iwa_hipp_dep_genes),
                                      area2 = length(iwa_roca_dep_genes),
                                      area3 = length(wal_dep_genes),
                                      n123 = length(overlap_123_d),
                                      n12 = length(overlap_12_d),
                                      n13 = length(overlap_13_d),
                                      n23 = length(overlap_23_d),
                                      category = c("Iwasaki et al. (Hipp)", "Iwasaki et al. (RocA)", "Waldron et al."),
                                      lty = rep("blank",3),
                                      fill = c("light blue", "red", "yellow"),
                                      alpha = rep(0.3, 3), 
                                      euler.d = T,
                                      scaled = T,
                                      cex = 2,
                                      cat.cex = 2,
                                      overrideTriple = T,
                                      cat.dist = c(0.035, 0.035, 0.04),
                                      cat.pos = c(350, 10, 180))


#eIF4A-antidependent data plot---
overlap_123_a <- factor(intersect(intersect(iwa_hipp_antidep_genes, iwa_roca_antidep_genes), wal_antidep_genes))
overlap_12_a <- factor(intersect(iwa_hipp_antidep_genes, iwa_roca_antidep_genes))
overlap_13_a <- factor(intersect(iwa_hipp_antidep_genes, wal_antidep_genes))
overlap_23_a <- factor(intersect(iwa_roca_antidep_genes, wal_antidep_genes))

grid.newpage()
venn_plot_antidep <- draw.triple.venn(area1 = length(iwa_hipp_antidep_genes),
                                          area2 = length(iwa_roca_antidep_genes),
                                          area3 = length(wal_antidep_genes),
                                          n123 = length(overlap_123_a),
                                          n12 = length(overlap_12_a),
                                          n13 = length(overlap_13_a),
                                          n23 = length(overlap_23_a),
                                          category = c("Iwasaki et al. (Hipp)", "Iwasaki et al. (RocA)", "Waldron et al."),
                                          lty = rep("blank",3),
                                          fill = c("green", "orange3", "purple"),
                                          alpha = rep(0.3, 3),
                                          euler.d = T,
                                          scaled = T,
                                          cex = 2,
                                          cat.cex = 2,
                                          overrideTriple = T,
                                          cat.dist = c(0.035, 0.035, 0.03),
                                          cat.pos = c(0, 0, 180))

