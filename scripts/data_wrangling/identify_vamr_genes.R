"""
This script detects viral antimicrobial resistance genes (vAMRs) by linking virus data with AMR gene data.
It processes both relaxed and strict quality viral sequences and generates various summary tables for further analysis.

Usage: source('identify_vamr_genes.R')
"""

library(tidyr)
library(ggplot2)
library(dplyr)
library(networkD3)
library(reshape2)
library(stringr)

# Set working directory
setwd("data_directory_path")

# Load virus data
virus <- read.csv("virus/overview_header_clean.tsv", sep="\t")
virus$contig_id <- as.character(virus$contig_id)
virus$sample <- as.character(virus$sample)
virus <- virus %>% filter(quality %in% c('Medium-quality', 'High-quality', 'Complete'))
virus_strict <- virus %>% filter(quality %in% c('High-quality', 'Complete'))

# Load AMRPP mapping of samples, genes, and contigs
amrpp_contigs_mapping <- read.csv("detecting_vAMRs/contigs_amrpp-complete-hit.csv", header=FALSE)
colnames(amrpp_contigs_mapping) <- c("sample", "amrpp_gene", "contig_id")
amrpp_contigs_mapping$contig_id <- gsub(" ", "", as.character(amrpp_contigs_mapping$contig_id))
amrpp_contigs_mapping$sample <- as.character(amrpp_contigs_mapping$sample)
amrpp_contigs_mapping <- unique(amrpp_contigs_mapping)
amrpp <- amrpp_contigs_mapping

# Load actual AMRPP hits
amrpp_actual_hits <- read.csv("amr/tools/amrpp/new/AMR_analytic_matrix.csv")
amrpp_actual_hits_melt <- melt(amrpp_actual_hits)
amrpp_actual_hits_melt <- filter(amrpp_actual_hits_melt, value > 0)

# Filter AMRPP data
filt_amrpp <- data.frame(sample = character(), amr = character(), contig = character())
Sys.time()
for (i in 1:nrow(amrpp_actual_hits_melt)) {
  sample_aux <- amrpp_actual_hits_melt$variable[i]
  amr <- amrpp_actual_hits_melt$gene_accession[i]
  chunk <- amrpp %>% filter(sample == sample_aux & amrpp_gene == amr)
  filt_amrpp <- rbind(filt_amrpp, chunk)
}
Sys.time()
amrpp <- filt_amrpp
rm(amrpp_actual_hits, amrpp_actual_hits_melt, chunk)

# Link AMRPP with virus data (relaxed)
amrpp_link_amrpp_virus <- merge(amrpp, virus, by=c("sample", "contig_id"))
amrpp_sample_amr_city <- amrpp_link_amrpp_virus[c('sample', 'amrpp_gene')]
amrpp_sample_amr_city$city <- sapply(strsplit(amrpp_sample_amr_city$sample, "_"), function(x) x[2])

# Generate summary tables
amrpp_all_amrs <- data.frame(as.character(unique(amrpp_sample_amr_city$amrpp_gene)))
amrpp_amrs_city_wise <- amrpp_sample_amr_city %>%
  group_by(city) %>%
  summarise(amrpp_genes = toString(unique(amrpp_gene)))
amrpp_amrs_sample_wise <- amrpp_sample_amr_city %>%
  group_by(sample) %>%
  summarise(amrpp_genes = toString(unique(amrpp_gene)))

dir.create('detecting_vAMRs/amrpp')
write.table(amrpp_all_amrs, "detecting_vAMRs/amrpp/vAMRS_all.csv", sep=',', row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(amrpp_amrs_city_wise, "detecting_vAMRs/amrpp/vAMRs_city_wise.csv", sep=',', row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(amrpp_amrs_sample_wise, "detecting_vAMRs/amrpp/vAMRs_sample_wise.csv", sep=',', row.names=FALSE, col.names=FALSE, quote=FALSE)

# Link AMRPP with virus data (strict)
virus <- virus_strict
amrpp_link_amrpp_virus <- merge(amrpp, virus, by=c("sample", "contig_id"))
amrpp_sample_amr_city <- amrpp_link_amrpp_virus[c('sample', 'amrpp_gene')]
amrpp_sample_amr_city$city <- sapply(strsplit(amrpp_sample_amr_city$sample, "_"), function(x) x[2])

# Generate strict summary tables
amrpp_all_amrs <- data.frame(as.character(unique(amrpp_sample_amr_city$amrpp_gene)))
colnames(amrpp_all_amrs) <- 'vamrs'
vamrs <- amrpp_all_amrs %>%
  separate(vamrs, into = c("meg", "class", "subclass", "descr", "gene", "conf"), sep = "\\|")
vamrs_list <- data.frame(unique(vamrs$gene))

write.table(vamrs_list, "detecting_vAMRs/vAMRS_gene_list_flexible.txt", row.names = FALSE, col.names = FALSE)

amrpp_amrs_city_wise <- amrpp_sample_amr_city %>%
  group_by(city) %>%
  summarise(amrpp_genes = toString(unique(amrpp_gene)))
amrpp_amrs_sample_wise <- amrpp_sample_amr_city %>%
  group_by(sample) %>%
  summarise(amrpp_genes = toString(unique(amrpp_gene)))

write.table(amrpp_all_amrs, "detecting_vAMRs/amrpp/vAMRS_all_strict.csv", sep=',', row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(amrpp_amrs_city_wise, "detecting_vAMRs/amrpp/vAMRs_city_wise_strict.csv", sep=',', row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(amrpp_amrs_sample_wise, "detecting_vAMRs/amrpp/vAMRs_sample_wise_strict.csv", sep=',', row.names=FALSE, col.names=FALSE, quote=FALSE)