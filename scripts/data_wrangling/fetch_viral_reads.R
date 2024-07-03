"""
This script processes mapped read data to fetch reads that are mapped to viral contigs. 
It merges these reads with a master table of viral contigs and writes the resulting data to a CSV file.

Usage: source('fetch_viral_reads.R')
"""

# Set working directory
setwd("data_directory_path/fetching_read_mapped_to_virus")

# Load the viral contigs master table
virus_mastertable <- read.csv("data_directory_path/publication/viral_investigation/overview_header_clean.tsv", sep = "\t")

# Initialize an empty dataframe to store viral reads
viral_reads <- data.frame(sample = character(0), contig = character(0), reads = character(0))

# Get the list of files in the 'mapped_reads' directory
file_list <- list.files("mapped_reads")

# Process each file in the list
for (file in file_list) {
  sample <- sub("_aligned_reads\\.txt$", "", file)  # Extract sample name from the file name
  data <- read.csv(paste0("mapped_reads/", file), sep = '\t')  # Read the file
  colnames(data) <- c("read", "contig_id")  # Rename columns
  data$sample <- sample  # Add sample column to the data
  
  # Merge with the viral contigs master table
  chunk <- merge(virus_mastertable, data, by = c('contig_id', 'sample'))
  chunk_clean <- chunk[, c("sample", "contig_id", "read")]  # Select relevant columns
  
  # Append the chunk to the viral_reads dataframe
  viral_reads <- rbind(viral_reads, chunk_clean)
}

# Write the viral reads data to a CSV file
write.csv(viral_reads, "viral_reads.csv", row.names = FALSE)