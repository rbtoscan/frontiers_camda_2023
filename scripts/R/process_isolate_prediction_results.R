"""
This script processes the results of isolate prediction using Boruta feature selection method. 
It calculates the mean probabilities for predictions, identifies the top predictions, and prints them.

Usage: Rscript process_isolate_prediction_results.R <path_to_prediction_file>
"""

args <- commandArgs(trailingOnly = TRUE)
file_prob <- args[1]

# Function to process the first set of results
result_1 <- function(file_prob) {
  prob <- read.csv(file_prob)
  
  # Identify numeric columns and filter them
  numeric_columns <- sapply(prob, is.numeric)
  prob_numeric <- prob[, numeric_columns]
  
  vsall <- prob_numeric[, 1:6]
  vsall_column_means <- data.frame(colMeans(vsall, na.rm = TRUE))
  
  df <- vsall_column_means
  df$pair <- rownames(df)
  
  # Find the row with the highest mean probability
  highest_row <- df[which.max(df$colMeans.vsall..na.rm...TRUE.), ]
  
  # Get the top 3 rows with the highest mean probabilities
  top3_rows <- df[order(df$colMeans.vsall..na.rm...TRUE., decreasing = TRUE)[1:3], ]
  colnames(top3_rows)[colnames(top3_rows) == 'colMeans.vsall..na.rm...TRUE.'] <- 'prob'
  
  # Clean up the 'pair' column values
  top3_rows$pair <- gsub("\\.versus\\.all", "", top3_rows$pair)
  
  # Print the top 3 pairs and their probabilities
  apply(top3_rows[, 2:1], 1, function(x) cat(x[1], x[2], sep = "\t", "\n"))
}

# Function to process the second set of results
result_2 <- function(file_prob) {
  prob <- read.csv(file_prob)
  
  # Identify numeric columns and filter them
  numeric_columns <- sapply(prob, is.numeric)
  prob_numeric <- prob[, numeric_columns]
  
  onevsone <- prob_numeric[, 7:21]
  onevsone_column_means <- data.frame(colMeans(onevsone, na.rm = TRUE))
  
  df <- onevsone_column_means
  df$pair <- rownames(df)
  
  # Clean up the 'pair' column values
  df$pair <- gsub("X\\.+", "", df$pair)
  df$pair <- gsub("\\.\\.+", "", df$pair)
  
  # Separate the 'pair' column into two columns 'a' and 'b'
  df2 <- df %>% separate(pair, into = c("a", "b"), sep = "\\.")
}

# Process the first set of results
result_1(file_prob)