#!/usr/bin/env Rscript

library(Biostrings)
library(ggplot2)
library(dplyr)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Ensure at least two arguments are provided
if (length(args) < 2) {
    stop("Not enough arguments. Please provide at least one FASTA file path and an output file name.")
}

# Last argument is the output file name
outputName <- args[length(args)]

# All other arguments are FASTA files
fasta_files <- args[-length(args)]

# Function to read a FASTA file and return a data frame with contig lengths and file name
read_contig_lengths <- function(file_path) {
  fasta <- readDNAStringSet(file_path)
  contig_lengths <- width(fasta)
  data.frame(length = contig_lengths, file = basename(file_path))
}

breaks <- c(0, 150, 250, 350, 450, 550, 650, 750, 850, 950, 1050, 1150, 1250, 1350, 1450, 1550,
1650, 1750, 1850, 1950, 2050, 2150, 2250, 2350, 2450, 2550, 2650, 2750, 2850, 2950, 3050, 3150,
 3250, 3350, 3450, 3550, 3650, 3750, 3850, 3950, 4050, 4150, 4250, 4350, 4450, 4550, 4650, 4750, 4850, 4950, 5050)
# Read all files and combine into one data frame
all_lengths <- do.call(rbind, lapply(fasta_files, read_contig_lengths))

all_lengths$first_bin <- ifelse(all_lengths$length >= breaks[1] & all_lengths$length <= breaks[2], "First Bin", "Other Bins")

# Plot the histogram with facets
ggplot(all_lengths, aes(x = length, fill = first_bin, color=first_bin)) +
    geom_histogram(color = "black", alpha = 0.6, breaks = breaks) +
    scale_fill_manual(values = c("First Bin" = "red", "Other Bins" = "blue")) +
   # scale_y_log10() +
    facet_wrap(~ file, scales = "fixed") +
    xlab("Contig Length") +
    ylab("Frequency") +
    ggtitle("Contig Length Distribution by Assembler") +
    theme(legend.position = "none")

ggsave(outputName, width = 10, height = 6)


ggplot(all_lengths, aes(x = length, fill = first_bin, color=first_bin)) +
    geom_histogram(color = "black", alpha = 0.6, breaks = breaks) +
    scale_fill_manual(values = c("First Bin" = "red", "Other Bins" = "blue")) +
    scale_y_log10() +
    facet_wrap(~ file, scales = "fixed") +
    xlab("Contig Length") +
    ylab("Frequency") +
    ggtitle("Contig Length Distribution by Assembler") +
    theme(legend.position = "none")

ggsave(paste0("log10_",outputName), width = 10, height = 6)
