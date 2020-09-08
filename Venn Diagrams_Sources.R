#VENN DIAGRAMS_SOURCES
#This script has been created by Rocco Roberto Penna with the support of Dr. Joseph A. Waldron

#AIM: to detect overlaps among studies
#OBJECTIVE: to obtain the list of transcripts for each study


#R packages---
library(dplyr)
library(readxl)

#Waldron Data---
wal_data <- read.csv("Files/Wal_penn-DOD-gene.mmdiff90.csv", header = T) #importing data

wal_data %>%
  mutate(external_gene_name = factor(toupper(external_gene_name)))%>% #converting transcript names to uppercase
  filter(posterior_probability > 0.25 & eta1_1 - eta1_2 < 0) %>% #filtering for eIF4AI-dependence conditions (post prob > 0.25 & DOD < 0)
  rename(gene_name = external_gene_name)%>% #changing column name
  mutate(gene_name = factor(gene_name))%>% #mutate column in factor
  pull(gene_name) %>% #isolating transcript names
  unique()-> wal_dep_genes #removing copies. 713 transcript names

wal_data %>%
  mutate(external_gene_name = factor(toupper(external_gene_name)))%>% #converting transcript names to uppercase
  filter(posterior_probability > 0.25 & eta1_1 - eta1_2 > 0) %>% #filtering for eIF4AI-antidependency conditions (post prob > 0.25 & DOD > 0)
  rename(gene_name = external_gene_name) %>% #changing column name
  mutate(gene_name = factor(gene_name))%>% #mutate column in factor
  pull(gene_name) %>% #isolating transcrip names
  unique()-> wal_antidep_genes #removing copies. 88 transcript names


#Modelska Data---
mod_dep <- read.table("Files/Mod_4A1-dep.txt", header = T) #eIF4A-dependent data import
mod_antidep <- read.table("Files/Mod_4A1-indep.txt", header = T) #eIF4A-antidependent data import

mod_dep %>%
  rename(gene_name = Gene_name) %>%
  mutate(gene_name = factor(toupper(gene_name)))%>%
  pull(gene_name) %>%
  unique()-> mod_dep_genes #161 transcript names

mod_antidep %>%
  rename(gene_name = Gene_name) %>%
  mutate(gene_name = factor(toupper(gene_name)))%>%
  pull(gene_name) %>%
  unique()-> mod_antidep_genes #49 transcript names


#Chan Data---
chan_dep <- read_xlsx("Files/Chan_et al._2019.xlsx", sheet = 2, col_names = T) #eIF4A-dependent data import
chan_antidep <- read_xlsx("Files/Chan_et al._2019.xlsx", sheet = 1, col_names = T) #eIF4A-antidependent data import

chan_dep %>%
  rename(gene_name = "Gene ID") %>%
  mutate(gene_name = factor(toupper(gene_name)))%>%
  pull(gene_name) %>%
  unique()-> chan_dep_genes #1168 transcript names


chan_antidep %>%
  rename(gene_name = "Gene ID") %>%
  mutate(gene_name = factor(toupper(gene_name)))%>%
  pull(gene_name) %>%
  unique()-> chan_antidep_genes #1060 transcript names


#Iwasaki Data---
#RocA Treatment
iwa_roca_dep <- read_xlsx("Files/Iwasaki_et al._2016_RocA.xlsx", sheet = 1, col_names = T) #eIF4A-dependent data import
iwa_roca_antidep <- read_xlsx("Files/Iwasaki_et al._2016_RocA.xlsx", sheet = 2, col_names = T) #eIF4A-antidependent data import

iwa_roca_dep %>%
  rename(gene_name = Gene) %>%
  mutate(gene_name = factor(toupper(gene_name)))%>%
  pull(gene_name) %>%
  unique()-> iwa_roca_dep_genes #290 transcript names

iwa_roca_antidep %>%
  rename(gene_name = Gene) %>%
  mutate(gene_name = factor(toupper(gene_name)))%>%
  pull(gene_name) %>%
  unique()-> iwa_roca_antidep_genes #307 transcript names

#Hippuristanol Treatment
iwa_hipp_dep <- read_xlsx("Files/Iwasaki_et al._2016_Hipp.xlsx", sheet = 1, col_names = T) #eIF4A-dependent data import
iwa_hipp_antidep <- read_xlsx("Files/Iwasaki_et al._2016_Hipp.xlsx", sheet = 2, col_names = T) #eIF4A-antidependent data import

iwa_hipp_dep %>%
  rename(gene_name = Gene) %>%
  mutate(gene_name = factor(toupper(gene_name)))%>%
  pull(gene_name) %>%
  unique()-> iwa_hipp_dep_genes #415 transcript names

iwa_hipp_antidep %>%
  rename(gene_name = Gene) %>%
  mutate(gene_name = factor(toupper(gene_name)))%>%
  pull(gene_name) %>%
  unique()-> iwa_hipp_antidep_genes #383 transcript names


#Rubio Data---
rub_dep <- read.table("Files/Rubio_dep.txt", header = T) #eIF4A-dependent data import
rub_antidep <- read.table("Files/Rubio_antidep.txt", header = T) #eIF4A-antidependent data import

rub_dep %>%
  rename(gene_name = Refseq_name) %>%
  mutate(gene_name = factor(toupper(gene_name)))%>%
  pull(gene_name) %>%
  unique()-> rub_dep_genes #284 transcript names

rub_antidep %>%
  rename(gene_name = Refseq_name) %>%
  mutate(gene_name = factor(toupper(gene_name)))%>%
  pull(gene_name) %>%
  unique()-> rub_antidep_genes #146 transcript names


#Wolfe Data---
wolfe_dep <- read_excel("Files/Wolfe et al._2014.xlsx", range = "A5:D286", col_names = T) #eIF4A-dependent data import
wolfe_antidep <- read_excel("Files/Wolfe et al._2014.xlsx", range = "A290:D480", col_names = T) #eIF4A-antidependent data import

wolfe_dep %>%
  rename(gene_name = "Gene name") %>%
  mutate(gene_name = factor(toupper(gene_name)))%>%
  pull(gene_name) %>%
  unique()-> wolfe_dep_genes #281 transcript names

wolfe_antidep %>%
  rename(gene_name = "Gene name") %>%
  mutate(gene_name = factor(toupper(gene_name)))%>%
  pull(gene_name) %>%
  unique()-> wolfe_antidep_genes #190 transcript names

#Removing unneccessary data (optional)
rm(chan_antidep, chan_dep, wal_data, mod_antidep, mod_dep, iwa_hipp_antidep, iwa_hipp_dep, iwa_roca_antidep, iwa_roca_dep, rub_antidep, rub_dep, wolfe_dep, wolfe_antidep)

