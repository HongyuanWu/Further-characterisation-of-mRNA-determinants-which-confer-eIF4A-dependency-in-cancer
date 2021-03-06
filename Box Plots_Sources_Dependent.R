#BOX PLOTS_EIF4A-DEPENDENT DATA
#This script has been created by Rocco Roberto Penna with the support of Dr. Joseph A. Waldron

#AIM: to compare mRNA features generally linked to eIF4A-dependency with eIF4A-antidependent data
#OBJECTIVE: to import eIF4A-dependent data from all studies


#R packages---
library(dplyr)
library(readxl)


#Data import and manipulation---
#Obtained list of transcript names for each eIF4A dependent data set

##Chan data
chan_dep <- read_xlsx("Files/Chan_et al._2019.xlsx", sheet = 2, col_names = T) #eIF4A-dep data import. 1168 transcripts

chan_dep %>%
  rename(gene_name = `Gene ID`) %>% #renaming gene name column
  mutate_all(funs(toupper)) %>% #converting gene names to upper case
  select(gene_name) -> chan_filtered #isolating list of gene names. 1168 transcripts

##Waldron data
wal_data <- read.csv("Files/Wal_penn-DOD-gene.mmdiff90.csv", header = T) #57177 transcripts

wal_data %>%
  filter(posterior_probability > 0.25 & eta1_1 - eta1_2 <0) %>%  #filtering for 4AI-dependent values
  rename(gene_id_v = feature_id, gene_name = external_gene_name) %>%
  select(gene_name) -> wal_filtered #713 transcripts

##Modelska data
mod_dep <- read.table("Files/Mod_4A1-dep.txt", header = T) #175 transcripts

mod_dep %>%
  group_by(Gene_name) %>% #for the repeated genes, one transcript is chosen
  sample_n(size = 1) %>%
  ungroup() %>%
  rename(transcript = Transcript_ID, gene_name = Gene_name) %>%  
  select(gene_name) -> mod_filtered #161 transcripts

##Wolfe data
wolf_dep <- read_xlsx("Files/Wolfe et al._2014.xlsx", range = "A5:D286", col_names = T) #281 transcripts

wolf_dep %>%
  rename(gene_name = `Gene name`) %>%
  select(gene_name) -> wolf_filtered #281 transcripts

##Rubio data
rub_dep <- read.table("Files/Rubio_dep.txt", header = T) #284 transcripts

rub_dep %>%
  rename(gene_name = Refseq_name) %>%
  select(gene_name) -> rub_filtered #284 transcripts

##Iwasaki data
iwa_hipp_dep <- read_xlsx("Files/Iwasaki_et al._2016_Hipp.xlsx", sheet = 1, col_names = T) #1039 transcripts

iwa_hipp_dep %>%
  rename(gene_name = Gene) %>%
  group_by(gene_name) %>%
  sample_n(size = 1) %>%
  ungroup() %>%
  select(gene_name) -> iwa_hipp_filtered #415 transcripts

iwa_roca_dep <- read_xlsx("Files/Iwasaki_et al._2016_RocA.xlsx", sheet = 1, col_names = T) #697 transcripts

iwa_roca_dep %>%
  rename(gene_name = Gene) %>%
  group_by(gene_name) %>%
  sample_n(size = 1) %>%
  ungroup() %>%
  select(gene_name) -> iwa_roca_filtered #290 transcripts


#Most abundant---
#Data from MCF7 cells
#Differential expression (DE) analysis and gene list files are loaded on R

##Data import
totals_data <- read.csv("Files/penn-DE.mmdiff90.csv", header = T) #DE file (by Waldron et al. 2019)
genes <- read.table("Files/martquery_0610114411_335.txt", header = T, fill = T)  #gene to transcript file

##Most abundant data isolation
totals_data %>%
  mutate(abundance = case_when(posterior_probability > 0.25 ~ alpha1,
                               posterior_probability < 0.25 ~ alpha0)) %>%
  rename(transcript_v = feature_id) %>%
  select(transcript_v, abundance, external_gene_name) -> abundance_data
rm(totals_data)

##Renaming variables in "genes" data frame
genes %>%
  rename(gene_id = Gene_stable_ID, gene_id_v = Gene_stable_ID_version, transcript = Transcript_stable_ID, transcript_v = Transcript_stable_ID_version, gene_name = Gene_name) %>%
  select(gene_id, gene_id_v, transcript, transcript_v, gene_name) -> genes

##Isolating most abundant transcript in MCF7 cells
abundance_data %>%
  group_by(external_gene_name) %>% #transcripts are grouped according to the common gene
  top_n(n = 1, wt = abundance) %>% # most abundant transcript is selected
  ungroup() -> most_abundant #most abundant transcript isolated
rm(abundance_data, genes)

most_abundant <- most_abundant[order(most_abundant$external_gene_name),] #the most_abundant data is alphabetically ordered


#ENSEMBL transcripts---
chan_filtered %>%
  left_join(most_abundant, by = c("gene_name" = "external_gene_name")) %>% #filtering for most abundant data
  na.omit() %>% #omitting NA values. 117 NA values
  select(transcript_v) -> chan_final_d #isolating ENSEMBL transcripts. 1051 transcripts

wal_filtered %>%
  left_join(most_abundant, by = c("gene_name" = "external_gene_name")) %>%
  na.omit() %>% #0 NA values
  select(transcript_v) -> wal_final_d #713 transcripts

mod_filtered %>%
  left_join(most_abundant, by = c("gene_name" = "external_gene_name")) %>%
  na.omit() %>% #8 NA values
  select(transcript_v) -> mod_final_d #153 transcripts

wolf_filtered %>%
  left_join(most_abundant, by = c("gene_name" = "external_gene_name")) %>%
  na.omit() %>% #15 NA values
  select(transcript_v) -> wolf_final_d #266 transcripts

rub_filtered %>%
  left_join(most_abundant, by = c("gene_name" = "external_gene_name")) %>%
  na.omit() %>% #27 NA values
  select(transcript_v) -> rub_final_d #257 transcripts

iwa_hipp_filtered %>%
  left_join(most_abundant, by = c("gene_name" = "external_gene_name")) %>%
  na.omit() %>% #10 NA values
  select(transcript_v) -> iwa_hipp_final_d #405 transcripts

iwa_roca_filtered %>%
  left_join(most_abundant, by = c("gene_name" = "external_gene_name")) %>%
  na.omit() %>% #6 NA values
  select(transcript_v) -> iwa_roca_final_d #284 transcripts


#Unneccessary files---
rm(chan_dep, chan_filtered)
rm(wal_data, wal_filtered)
rm(mod_dep, mod_filtered)
rm(wolf_dep, wolf_filtered)
rm(rub_dep, rub_filtered)
rm(iwa_hipp_dep, iwa_hipp_filtered, iwa_roca_dep, iwa_roca_filtered)

