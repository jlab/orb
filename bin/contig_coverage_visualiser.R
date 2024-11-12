#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

# Assuming the last argument is the output file name
outputName <- args[length(args)]

# Read each file, add a column to identify the file, and combine into one data frame
all_data <- lapply(args[-length(args)], function(file) {
  data <- read.csv(file, sep = "\t", header = TRUE)
  data$Name <- as.integer(as.character(data$Name)) # Convert factor to integer if needed
  data$file <- basename(file) # Add a column with the file name
  return(data)
}) %>% bind_rows()

p <- ggplot(all_data, aes(x = TPM)) +
    geom_histogram(fill = "blue", color = "black", alpha=0.6) +
    facet_wrap(~ file, scales = "fixed") +
    theme_minimal() +
    scale_y_log10() +
    #xlim(c(NA, 10000)) +
    #ylim(c(NA, 100)) +
    labs(title = "Histogram of TPM",
        x = "TPM",
        y = "Frequency")



# Save the plot
ggsave(paste0("y_log10_", outputName), plot = p, width = 15, height = 10, dpi = 600)


# Plot with facet wrap

p <- ggplot(all_data, aes(x = TPM)) +
    geom_histogram(fill = "blue", color = "black", alpha=0.6) +
    facet_wrap(~ file, scales = "fixed") +
    theme_minimal() +
    scale_y_log10() +
    xlim(c(NA, 10000)) +
    #ylim(c(NA, 100)) +
    labs(title = "Histogram of TPM",
        x = "TPM",
        y = "Frequency")

# Save the plot
ggsave(paste0("y_log_10_x_limit_10000_", outputName), plot = p, width = 15, height = 10, dpi = 600)

p <- ggplot(all_data, aes(x = TPM)) +
    geom_histogram(fill = "blue", color = "black", alpha=0.6) +
    facet_wrap(~ file, scales = "fixed") +
    theme_minimal() +
    scale_y_log10() +
    xlim(c(NA, 1000)) +
    #ylim(c(NA, 100)) +
    labs(title = "Histogram of TPM",
        x = "TPM",
        y = "Frequency")

# Save the plot
ggsave(paste0("y_log_10_x_limit_1000_", outputName), plot = p, width = 15, height = 10, dpi = 600)

p <- ggplot(all_data, aes(x = TPM)) +
    geom_histogram(fill = "blue", color = "black", alpha=0.6) +
    facet_wrap(~ file, scales = "fixed") +
    theme_minimal() +
    xlim(c(NA, 100)) +
    ylim(c(NA, 100)) +
    labs(title = "Histogram of TPM",
        x = "TPM",
        y = "Frequency")

# Save the plot
ggsave(paste0("y_limit_100_x_limit_100_", outputName), plot = p, width = 15, height = 10, dpi = 600)
