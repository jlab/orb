library(edgeR)

args <- commandArgs(trailingOnly = TRUE)

counts_file <- args[1]
prefix <- args[2]

group1_name <- args[3]
group1_replicates <- as.numeric(args[4])

group2_name <- args[5]
group2_replicates <- as.numeric(args[6])

counts <- read.csv(counts_file, sep = "\t", skip = 1, row.names = 1)

group1_cols <- names(counts)[grepl(group1_name, names(counts))]
group2_cols <- names(counts)[grepl(group2_name, names(counts))]

groups <- factor(c(rep(group1_name, group1_replicates), rep(group2_name, group1_replicates)))

dge <- DGEList(counts = counts[c(group1_cols, group2_cols)], group = groups)

dge <- calcNormFactors(dge)

dge <- estimateDisp(dge)

fit <- glmFit(dge)
lrt <- glmLRT(fit)

topTags(lrt)

estimated_fold_changes <- 2 ^ lrt$table$logFC

lrt$table["AveLogCPM"] <- lrt$AveLogCPM
lrt$table["dispersion"] <- lrt$dispersion
lrt$table["residual"] <- lrt$residual
lrt$table["coefficients"] <- lrt$coefficients
lrt$table["deviance"] <- lrt$deviance
lrt$table <- cbind(ContigName = rownames(lrt$table), lrt$table)

write.table(lrt$table, file = paste0(prefix, "_edgeR_full_table.tsv"), sep = "\t", quote = F, row.names = F)

edgeR_version <- packageVersion("edgeR")
cat(sprintf('edgeR: "%s"\n', edgeR_version), file = "versions.yml")