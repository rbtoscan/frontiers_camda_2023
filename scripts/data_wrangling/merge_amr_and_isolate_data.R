"""
This script joins AMR data with isolates data to create a combined dataset for further analysis. 
It processes both datasets, transposes and collapses them, and writes the resulting datasets to CSV files.

Usage: source('merge_amr_and_isolate_data.R')
"""

library(dplyr)
library(stringr)
library(tidyr)

setwd('data_directory_path/join_amr_isolates')

# Load data
isolates_matrix <- read.csv('data/isolates/isolates_matrix.csv')
amrpp_isolates_mapping_file <- read.csv('data/amr/amr_mapping_file.csv', sep=';')

# Transpose isolates matrix
isolates_matrix_t <- t(isolates_matrix)
colnames(isolates_matrix_t) <- as.character(isolates_matrix_t[1, ])  # Convert the first row to character and assign as column names
isolates_matrix_t <- data.frame(isolates_matrix_t[-1, ])
isolates_matrix_t$isolates <- rownames(isolates_matrix_t)

# Merge with AMRPP isolates mapping
isolates_matrix_amrpp <- merge(isolates_matrix_t, amrpp_isolates_mapping_file, by='isolates')
isolates_matrix_amrpp <- subset(isolates_matrix_amrpp, select = -isolates)
colnames(isolates_matrix_amrpp)[colnames(isolates_matrix_amrpp) == 'amrpp'] <- 'gene_accession'
isolates_matrix_amrpp <- isolates_matrix_amrpp %>% select(gene_accession, everything())
isolates_matrix_amrpp_genes <- isolates_matrix_amrpp

# Clean up
rm(isolates_matrix_amrpp, isolates_matrix_t, isolates_matrix, amrpp_isolates_mapping_file)

# Load AMR data
amr_raw <- read.csv('data/amr/AMR_analytic_matrix.csv')

# Function to transpose and collapse data
transpose_collapse <- function(df) {
  df <- df %>% separate(gene_accession, into = c("meg", "class", "subclass", "descr", "gene", "conf"), sep = '\\|')
  df <- df[, !colnames(df) %in% c("meg", "class", "subclass", "descr", "conf")]
  
  # Convert all columns to numeric except the 'gene' column
  df[, -which(names(df) == "gene")] <- lapply(df[, -which(names(df) == "gene")], function(x) as.numeric(as.character(x)))
  
  # Group by 'gene' and sum all other columns
  df_summed <- df %>%
    group_by(gene) %>%
    summarise(across(where(is.numeric), sum, na.rm = TRUE))
  
  # Binarize the summed data
  df_binarized <- df_summed %>%
    mutate(across(where(is.numeric), ~ ifelse(. > 0, 1, 0)))
  
  df_sum_bin_t <- t(df_binarized)
  colnames(df_sum_bin_t) <- as.character(df_sum_bin_t[1, ])
  df_sum_bin_t <- data.frame(df_sum_bin_t[-1, ])
  return (df_sum_bin_t)
}

# Function to add city information
add_city <- function(df) {
  df$info <- rownames(df)
  df <- df %>% separate(info, into = c("year", "city", "num"), sep = '_')
  df <- df[, !colnames(df) %in% c("year", "num")]
  df <- df %>%
    select(city, everything())
}

# Collapse isolate and AMR data
collap_isolate <- transpose_collapse(isolates_matrix_amrpp_genes)
collap_amr_city <- transpose_collapse(amr_raw)

# Write results to CSV
write.csv(collap_isolate, 'isolates.csv')
write.csv(collap_amr_city, 'training.csv')