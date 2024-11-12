#!/usr/bin/env Rscript

library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

inputFilePath <-  args[1]
outputName <- args[2]

data <- read.csv(inputFilePath, sep = "\t", header = TRUE)

p <- ggplot(data, aes(x = TPM)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha=0.6) +
  theme_minimal() +
  scale_y_log10() +
  labs(title = "Histogram of TPM",
       x = "TPM",
       y = "Frequency")

ggsave(outputName, plot = p, width = 15, height = 10, dpi = 600)
