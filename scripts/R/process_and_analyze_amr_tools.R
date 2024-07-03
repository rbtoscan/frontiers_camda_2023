"""

This script processes and analyzes AMR data from various tools, calculates exploratory statistics, 
generates figures, performs AUC calculations, and normalizes the data. It also predicts isolates 
using the Boruta feature selection method and saves the results.

Usage: source('process_and_analyze_amr_tools.R')
"""

library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(knitr)
library(kableExtra)

# Set working directory
setwd('data_directory_path/amr')

########################################### Read AMR tools and parse them ########################################### 

# Function to read and process RGI data
rgi_read <- function() {
  rgi <- read.csv('tools/rgi/rgi_all.txt', sep='\t')
  colnames(rgi)[colnames(rgi) == 'gCSD16_DEN_10'] <- 'sample'
  rgi <- subset(rgi, Best_Hit_Bitscore != "Best_Hit_Bitscore")
  rgi <- rgi[, c("sample", "Best_Hit_ARO")]
  colnames(rgi)[colnames(rgi) == 'Best_Hit_ARO'] <- 'gene'
  rgi <- rgi %>%
    group_by(sample, gene) %>%
    summarise(count = n()) %>%
    ungroup()
  df <- rgi %>%
    pivot_wider(names_from = gene, values_from = count, values_fill = list(count = 0))
  matrix_result <- as.matrix(df[,-1])
  rownames(matrix_result) <- df$sample
  matrix_rgi <- data.frame(matrix_result)
  return (matrix_rgi)
}

# Function to read and process AMRFinder data
amrfinder_read <- function() {
  amrfinder <- read.csv('tools/amrfinder/amrfinder_all_sort_uniq_header.txt', sep='\t')
  colnames(amrfinder)[colnames(amrfinder) == 'gCSD16_DEN_10'] <- 'sample'
  amrfinder <- amrfinder[, c("sample", "Gene.symbol")]
  colnames(amrfinder)[colnames(amrfinder) == 'Gene.symbol'] <- 'gene'
  amrfinder <- amrfinder %>%
    group_by(sample, gene) %>%
    summarise(count = n()) %>%
    ungroup()
  df <- amrfinder %>%
    pivot_wider(names_from = gene, values_from = count, values_fill = list(count = 0))
  matrix_result <- as.matrix(df[,-1])
  rownames(matrix_result) <- df$sample
  matrix_amrfinder <- data.frame(matrix_result)
  return (matrix_amrfinder)
}

# Function to read and process AMR++ data
amrpp_read <- function() {
  amrpp <- read.csv('tools/amrpp/new/AMR_analytic_matrix.csv')
  amrpp <- amrpp %>% separate(gene_accession, into = c("meg", "class", "subclass", "descr", "gene", "conf"), sep = "\\|")
  amrpp <- subset(amrpp, select = -c(meg, class, subclass, descr, conf))
  amrpp <- t(amrpp)
  colnames(amrpp) <- as.character(amrpp[1, ])
  amrpp <- data.frame(amrpp[-1, ])
  return (amrpp)
}

# Function to read and process Bowtie data
bowtie_read <- function() {
  bowtie <- read.csv('tools/bowtie/amrs_samples_bowtie-hits.tsv', sep='\t')
  rownames(bowtie) <- bowtie$X
  bowtie$X <- NULL
  bowtie <- data.frame(t(bowtie))
  return (bowtie)
}

# Digest RGI, AMRFinder, AMR++, and Bowtie data
rgi <- rgi_read()
amrfinder <- amrfinder_read()
amrpp <- amrpp_read()
bowtie <- bowtie_read()

################################ Calculating Exploratory Statistics ################################################## 

# Function to calculate gene statistics
calculate_gene_stats <- function(df) {
  df <- df %>% mutate(across(everything(), as.numeric))
  binary_df <- df %>% mutate(across(everything(), ~ ifelse(. != 0, 1, 0)))
  total_unique_genes <- sum(colSums(binary_df) > 0)
  unique_genes_per_sample <- rowSums(binary_df)
  avg_unique_genes <- mean(unique_genes_per_sample)
  sd_unique_genes <- sd(unique_genes_per_sample)
  return (data.frame(list(
    total_unique_genes = total_unique_genes,
    avg_unique_genes = avg_unique_genes,
    sd_unique_genes = sd_unique_genes
  )))
}

amrpp_stats <- calculate_gene_stats(amrpp)
bowtie_stats <- calculate_gene_stats(bowtie)
rgi_stats <- calculate_gene_stats(rgi)
amrfinder_stats <- calculate_gene_stats(amrfinder)

# Function to create a summary dataframe with tool names
create_summary_df <- function(stats, tool_name) {
  data.frame(
    tool = tool_name,
    total_unique_genes = stats$total_unique_genes,
    avg_unique_genes = stats$avg_unique_genes,
    sd_unique_genes = stats$sd_unique_genes
  )
}

amrpp_summary <- create_summary_df(amrpp_stats, "AMR++")
bowtie_summary <- create_summary_df(bowtie_stats, "Bowtie")
rgi_summary <- create_summary_df(rgi_stats, "RGI")
amrfinder_summary <- create_summary_df(amrfinder_stats, "AMRFinder")

# Combine all summary dataframes into one
amr_tools_statistics <- bind_rows(amrpp_summary, bowtie_summary, rgi_summary, amrfinder_summary)
rm(amrpp_summary, bowtie_summary, rgi_summary, amrfinder_summary, bowtie_stats, amrfinder_stats, rgi_stats, amrpp_stats)

########################################### Figures using AMR++ ###################################################### 

# Function to generate figure showing AMR++ class distribution city-wise normalized by sample reads
make_fig_amrpp_class_city_norm_reads <- function() {
  source('normalization/normalize_metrics.R')
  setwd('data_directory_path/amr')
  amrpp <- read.csv('tools/amrpp/new/AMR_analytic_matrix.csv')
  amrpp <- amrpp %>% separate(gene_accession, into = c("meg", "class", "subclass", "descr", "gene", "conf"), sep = "\\|")
  amrpp <- subset(amrpp, select = -c(meg, conf, descr))
  long_df <- amrpp %>%
    pivot_longer(
      cols = starts_with("gCSD"),
      names_to = "sample",
      values_to = "value"
    )
  long_df_enriched <- merge(long_df, norm_metrics, by='sample')
  long_df_enriched$gene_count_norm <- long_df_enriched$value / long_df_enriched$reads
  long_df <- long_df_enriched
  df <- long_df %>%
    separate(sample, into = c("year", "city", "number"), sep = "_")
  gene_counts <- df %>%
    group_by(city, class) %>%
    summarise(total_genes = sum(gene_count_norm, na.rm = TRUE)) %>%
    ungroup()
  plot <- ggplot(gene_counts, aes(x = city, y = total_genes, fill = class)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Total Number of Genes per Class in Each City",
         x = "City",
         y = "Total Genes",
         fill = "Class") +
    theme_minimal()
  return (plot)
}

fig_amrpp_class_city_norm_reads <- make_fig_amrpp_class_city_norm_reads()
print(fig_amrpp_class_city_norm_reads)
rm(fig_amrpp_class_city_norm_reads)

# Function to generate figure showing AMR++ subclass distribution city-wise normalized by sample reads
make_fig_amrpp_subclass_city_norm_reads <- function() {
  source('normalization/normalize_metrics.R')
  setwd('data_directory_path/amr')
  amrpp <- read.csv('tools/amrpp/new/AMR_analytic_matrix.csv')
  amrpp <- amrpp %>% separate(gene_accession, into = c("meg", "class", "subclass", "descr", "gene", "conf"), sep = "\\|")
  amrpp <- subset(amrpp, select = -c(meg, conf, descr))
  long_df <- amrpp %>%
    pivot_longer(
      cols = starts_with("gCSD"),
      names_to = "sample",
      values_to = "value"
    )
  long_df_enriched <- merge(long_df, norm_metrics, by='sample')
  long_df_enriched$gene_count_norm <- long_df_enriched$value / long_df_enriched$reads
  long_df <- long_df_enriched
  df <- long_df %>%
    separate(sample, into = c("year", "city", "number"), sep = "_")
  gene_counts <- df %>%
    group_by(city, subclass) %>%
    summarise(total_genes = sum(gene_count_norm, na.rm = TRUE)) %>%
    ungroup()
  plot <- gg