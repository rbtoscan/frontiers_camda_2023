"""
This script parses AMR analytic matrices, transforms them by separating and reorganizing columns, 
and writes the processed data to CSV files for further analysis.

Usage: source('parse_amr_analytic_matrices.R')
"""

library(dplyr)
library(tidyr)
library(stringr)
library(readr)

# Set working directory
setwd('data_directory_path/comparing_amrpp_with_and_without_viral_reads')

# File paths
AMR_analytic_matrix_camda_all_reads <- "AMR_analytic_matrix_camda_all_reads.csv"
AMR_analytic_matrix_camda_nonviral_reads <- "AMR_analytic_matrix_camda_non-viral_reads.csv"

# Function to parse and transform the table
parse <- function(table) {
  df <- t(read_csv(table))
  colnames(df) <- as.character(df[1, ])  # Convert the first row to character and assign as column names
  df <- df[-1, ]  # Remove the first row from the data frame
  df <- data.frame(df)
  df$city <- row.names(df)
  df <- df %>% separate(city, into = c("year", "city", "index"), sep = "_")
  df <- subset(df, select = -c(year, index))
  df <- df[c(ncol(df), 1:(ncol(df) - 1))]
  return(df)
}

# Parse tables
table_full <- parse(AMR_analytic_matrix_camda_all_reads)
table_non_virus <- parse(AMR_analytic_matrix_camda_nonviral_reads)

# Write parsed tables to CSV
write.csv(table_full, "amrpp_full.csv", row.names = FALSE)
write.csv(table_non_virus, "amrpp_non-virus.csv", row.names = FALSE)