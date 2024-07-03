"""
This script normalizes various metrics from different input tables by merging them into a single dataframe. 
It removes unnecessary columns and calculates normalization metrics for further analysis.

Usage: source('normalize_metrics.R')
"""

library(dplyr)

# Set the working directory
setwd('data_directory_path/normalization')

# Load data
ssu_std <- read.csv("table_std.tsv", sep="\t")
ssu_cov50 <- read.csv("table_mincov_50.tsv", sep="\t")
assemblies <- read.csv("assemblies_stats_full-output.tsv", sep="\t")
reads <- distinct(read.csv("camda_reads_bases.tsv", sep="\t"))

# Merge data
df <- merge(merge(merge(ssu_std, ssu_cov50, by='sample'), assemblies, by='sample'), reads, by='sample')

# Remove unnecessary columns
df <- subset(df, select = -c(A..., T..., G..., C..., N50, Max, Min, N..., N....1, Gap..., Uncertain.bp., Ncount))

# Normalize metrics
norm_metrics <- subset(df, select = -c(Contigs, bases))

# Set row names and remove sample column
rownames(df) <- df$sample
df <- subset(df, select = -sample)

# Clean up
rm(assemblies, reads, ssu_cov50, ssu_std, df)