#!/usr/bin/env Rscript
# Summary: Build or append all_samples_annotation.txt from metadata and sample lists.
# This file for generating all_samples_annotation.txt, while existed, new information will be appended
# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check the number of arguments
if (length(args) != 4) {
  stop("Usage: Rscript 1.1_prepare_sample_annotation.R <current_path> <meta> <tumor> <bgfile>")
}

# Assign arguments to variables
current_path <- args[1]  # current_path = "DeepLB path:/home/yinliang/PROJECT/DeepLB"
meta <- args[2]         # all_WGBS_sample_metadata.xlsx
tumor <- args[3]
bgfile <- args[4]

# Set directory paths
dir1 <- paste0(current_path, "/Predata/metadata/")
dir2 <- paste0(current_path, "/Result/1.1_450K_reformate_data/")

# Load necessary libraries
library("readxl")
library(dplyr)

# Read WGBS_meta and background sample list
WGBS_meta <- read_excel(paste0(dir1, meta))
bg <- read.table(paste0(dir1, bgfile))

# Read index files
tumor_index <- read.table(paste0(dir2, "index_", tolower(tumor), "_tumor"))[, -1]
normal_index <- read.table(paste0(dir2, "index_", tolower(tumor), "_normal"))[, -1]

# View data headers
head(WGBS_meta)
head(tumor_index)
head(normal_index)
head(bg)

# Prepare background sample list
bg$V2 <- "plasma"
bg$V3 <- "background"
bg$V4 <- "normal"

# Prepare test sample list
test_sample <- WGBS_meta[, c("Run", "Sample type", "Group")]
test_sample$`Sample type` <- case_when(test_sample$`Sample type`  == "cfDNA" ~ "plasma", TRUE ~ "tissue")
test_sample$Group <- gsub(" ","_",test_sample$Group)
test_sample$V4 <- case_when(test_sample$Group == "Healthy" ~ "normal", TRUE ~ "non-healthy")
test_sample <- test_sample[!(test_sample$Run %in% bg$V1), ]

# Set column names
colnames(test_sample) <- c("V1", "V2", "V3", "V4")
colnames(tumor_index) <- c("V1", "V2", "V3", "V4")
colnames(normal_index) <- c("V1", "V2", "V3", "V4")
head(test_sample)
# Print dimensions of data
print(dim(test_sample))
print(dim(tumor_index))
print(dim(normal_index))
print(dim(bg))

# Combine all data
annotation <- rbind(bg, tumor_index, normal_index, test_sample)

# Output file path
output_file <- paste0(dir1, "all_samples_annotation.txt")

if (file.exists(output_file)) {
  # Read existing file content
  existing_annotation <- read.table(output_file, header = FALSE, stringsAsFactors = FALSE)
  
  # Convert new generated content to data frame
  new_annotation <- as.data.frame(annotation)
  
  # Compare the two versions and find new rows
  new_rows <- anti_join(new_annotation, existing_annotation, by = colnames(new_annotation))
  
  # Append new rows to the existing file
  if (nrow(new_rows) > 0) {
    write.table(new_rows, file = output_file,sep = '	', quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    cat("Added", nrow(new_rows), "new rows to", output_file, "\n")
  } else {
    cat("No new rows to add to", output_file, "\n")
  }
} else {
  # If the file does not exist, create a new file with the new content
  write.table(annotation,output_file,sep = '	', quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
  cat("Created new file:", output_file, "\n")
}
print(output_file)

