"""
This script processes multiple input files using Boruta for feature selection and processes results for prediction of isolate origins.
It sources necessary functions from external scripts, processes each file in the input directory, and gathers results for further analysis.

Usage: source('run_isolate_prediction.R')
"""

# Load necessary libraries
library(dplyr)
library(stringr)
library(tidyr)

# Source required scripts
source('woktek_boruta_results.r')
source('wojtek_mdfs_results.r')
source('results.R')

# Set the path to the input files
files_path <- "inputs_no_mge"

# List all files in the directory
files <- list.files(path = files_path, full.names = TRUE)

# Loop through the files and process each one
for (file_path in files) {
  boruta_wojtek(file_path)
  # Uncomment the following line if you also want to use MDFS
  # mdfs_wojtek(file_path)
}

# Set the path to the result directories
result_dirs <- "results"

# List all subdirectories in the result directory
dirs <- list.dirs(path = result_dirs, full.names = TRUE, recursive = TRUE)
file_paths <- c()

# Loop through each subdirectory and gather result files
for (subdir in dirs) {
  print(subdir)
  files <- list.files(subdir, full.names = TRUE)
  prob_files <- files[grep("^prob", basename(files))]
  file_paths <- c(file_paths, prob_files)
}

# Example of how to use the results
result_example <- "results/amrpp_norm_apr15_count_std/boruta/prob_prediction_result_boruta_amrpp_norm_apr15_count_std.csv"
# result_1(result_example)

# Assuming 'prob_numeric' is a data frame loaded from some file
onevsone <- prob_numeric[, 7:21]

# Calculate the sum of each numeric column
onevsone_column_sums <- data.frame(colSums(onevsone, na.rm = TRUE))

# Prepare data for further analysis
df <- onevsone_column_sums
df$pair <- rownames(df)
df$pair <- gsub("X\\.+", "", df$pair)
df$pair <- gsub("\\.\\.+", "", df$pair)

# Assuming 'df' is your data frame and 'pair' is the column you want to clean
# df$pair <- gsub("[^[:graph:][:al​⬤