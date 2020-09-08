#CORRELATION_WOLFE
#This script has been created by Rocco Roberto Penna with the support of Dr. Joseph A. Waldron

#AIM: to correlate 5'UTR GC content and length against changes in translational efficiency
#OBJECTIVE: to correlate 5'UTR GC content and length against changes in translational efficiency of data obtained by Wolfe et al. 2014


#R packages---
library(dplyr)
library(ggplot2)
library(readxl)


#Data import---
##Wolfe data
wolf_TE_down <- read_excel("Files/Wolfe et al._2014.xlsx", range = "A5:D286") #eIF4A-dep
wolf_TE_up <- read_excel("Files/Wolfe et al._2014.xlsx", range = "A290:D480") #eIF4A-antidep
wolf_TE_ind <- read_excel("Files/Wolfe et al._2014.xlsx", range = "A484:D5069") #eIF4A-indep

##5'UTR
fpUTR <- read.csv("Files/gencode_v27_pc_transcripts_fpUTRs_composition.csv", header = T) #5'UTR features


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


#Data manipulation---
wolf_TE_ind %>%
  rename(`Log2(Translational Efficiency)` = `log2(Translational Efficiency)`) -> wolf_TE_ind #Renaming column name
wolf_data <- rbind(wolf_TE_down, wolf_TE_up, wolf_TE_ind) #Merging data
rm(wolf_TE_down, wolf_TE_ind, wolf_TE_up) #Removing unneccessary data
wolf_data$`Gene name`<- factor(toupper(wolf_data$`Gene name`)) #Converting names to uppercase

wolf_data %>%
  left_join(most_abundant, by = c("Gene name" = "external_gene_name")) %>% #filtering for most abundant transcripts
  na.omit() %>%
  inner_join(fpUTR, by = c("transcript_v" = "transcript"))%>% #coupling 5'UTR features
  rename(logFC = `Log2(Translational Efficiency)`)%>%
  select(transcript_v, logFC, GC_content, length) -> wolfe_plot_data

#Correlation test---
#Funtion of label printing
myR <- function(x) {
  r <- round(as.numeric(x$estimate), digits = 2)
  p <- as.numeric(x$p.value)
  if (p > 0.001) {
    rounded_p <- round(p, digits = 3)
    p_label <- paste('P =', rounded_p)
  }else{
    if (p<2.2E-16) {
      p_label <- 'P < 2.2e-16'
    } else {
      rounded_p <- formatC(p, format = "e", digits = 2)
      p_label <- paste('P =', rounded_p)
    }
  }
  return(paste0('r = ', r, '\n', p_label))
}

#Two-sided Pearson's product-moment correlation
r_GC <- cor.test(x = wolfe_plot_data$logFC, y = wolfe_plot_data$GC_content) #GC content
r_label_GC <- myR(x = r_GC)

r_len <- cor.test(x = wolfe_plot_data$logFC, y = log10(wolfe_plot_data$length)) #length
r_label_len <- myR(x = r_len)

#Plot---
##Theme
myTheme <- theme_bw()+
  theme(axis.text = element_text(size=18), 
        axis.title = element_text(size=20), 
        legend.title = element_blank(),
        legend.text = element_text(size=18),
        plot.title = element_text(hjust = 1))

##GC content
cor_GC_wolf <- ggplot(wolfe_plot_data, aes(x = logFC, y = GC_content))+
  geom_point(size = 1)+
  ylim(0,1)+
  xlim(-1, 1)+
  xlab("Translation Efficiency\n(log-fold change)")+
  ylab("5'UTR GC content (%)")+
  myTheme+
  ggtitle(r_label_GC)+
  geom_smooth(method = "lm", se = F)

##length
cor_length_wolf <- ggplot(wolfe_plot_data, aes(x = logFC, y = length))+
  geom_point(size = 1)+
  scale_y_log10(breaks=c(10,100,1000),limits=c(10, 1500))+
  ylab("5'UTR length (nt, log10)") +
  xlab("Translation Efficiency\n(log-fold change)")+
  xlim(-1, 1)+
  myTheme+
  ggtitle(r_label_len)+
  geom_smooth(method = "lm", se = F)
  
