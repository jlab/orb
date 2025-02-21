library(DESeq2)

args <- commandArgs(trailingOnly = TRUE)

counts_file <- args[1]
prefix <- args[2]

group1_name <- args[3]
group1_replicates <- as.numeric(args[4])

group2_name <- args[5]
group2_replicates <- as.numeric(args[6])
separator <- args[7]

counts <- read.csv(counts_file, sep = separator, skip = 0, row.names = 1)
print(head(counts))

group1_cols <- names(counts)[grepl(group1_name, names(counts))][1:group1_replicates]
group2_cols <- names(counts)[grepl(group2_name, names(counts))][1:group2_replicates]

coldata <- data.frame(
  Name = c(group1_cols, group2_cols),
  Group = c(rep("group1", 10), rep("group2", 10)),
  row.names = c(group1_cols, group2_cols)
)

deseq <- DESeqDataSetFromMatrix(round(counts[c(group1_cols, group2_cols)]), coldata, ~Group)
dds <- DESeq(deseq)

res <- results(dds)

res$dispersion <- dispersions(dds)
res$deviance <- mcols(dds)$deviance
res$Intercept <- mcols(dds)$Intercept

res_df <- as.data.frame(res)
res_df <- cbind(ContigName = rownames(res_df), res_df)

write.table(res_df, file = paste0(prefix, "_DESeq2_full_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

deseq2_version <- packageVersion("DESeq2")
cat(sprintf('DESeq2: "%s"\n', deseq2_version), file = "versions.yml")
