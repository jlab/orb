library(edgeR)

args <- commandArgs(trailingOnly = TRUE)

counts_file <- args[1]
prefix <- args[2]

group1_name <- args[3]
group1_replicates <- as.numeric(args[4])

group2_name <- args[5]
group2_replicates <- as.numeric(args[6])

separator <- args[7]

counts <- read.csv(counts_file, sep = separator, skip = 0, row.names = 1)

group1_cols <- names(counts)[grepl(group1_name, names(counts))][1:group1_replicates]
group2_cols <- names(counts)[grepl(group2_name, names(counts))][1:group2_replicates]

groups <- factor(c(rep(group1_name, group1_replicates), rep(group2_name, group1_replicates)))

dge <- DGEList(counts = counts[c(group1_cols, group2_cols)], group = groups)

dge <- calcNormFactors(dge)

dge <- estimateDisp(dge)

fit <- glmFit(dge)
lrt <- glmLRT(fit)

# lrt$table["AveLogCPM"] <- lrt$AveLogCPM
# lrt$table["dispersion"] <- lrt$dispersion
# lrt$table["residual"] <- lrt$residual
# lrt$table["coefficients"] <- lrt$coefficients
# lrt$table["deviance"] <- lrt$deviance
# lrt$table <- cbind(ContigName = rownames(lrt$table), lrt$table)

res <- topTags(lrt, n = Inf)

res_table <- res$table
res_table <- cbind(ContigName = rownames(res_table), res_table)

write.table(res_table, file = paste0(prefix, "_edgeR_full_table.tsv"), sep = "\t", quote = F, row.names = F)

edgeR_version <- packageVersion("edgeR")
cat(sprintf('edgeR: "%s"\n', edgeR_version), file = "versions.yml")