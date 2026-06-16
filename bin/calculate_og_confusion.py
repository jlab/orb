#!/usr/bin/env python

import pandas as pd
import json
import sys
from sklearn.metrics import confusion_matrix

assembler_dge_file = sys.argv[1]
gene_summary_file = sys.argv[2]
prefix = sys.argv[3]

pval_name = sys.argv[4]
log_name = sys.argv[5]
row_prefix = sys.argv[6]

gene_summary = pd.read_csv(gene_summary_file, sep="\t")
assembler_dge = pd.read_csv(assembler_dge_file, sep="\t")

gene_summary["cds_de"] = (gene_summary[pval_name] < 0.05) & (gene_summary[log_name].abs() > 1)

assembler_dge["ContigName"] = assembler_dge["ContigName"].astype(str)

assembler_dge["de"] = (assembler_dge[pval_name] < 0.05) & (assembler_dge[log_name].abs() > 1)

cds_contig_dict = dict(zip(assembler_dge["ContigName"], assembler_dge["de"]))

gene_summary["assembler_de"] = gene_summary["ContigName"].apply(lambda x: cds_contig_dict.get(x, False))

true_labels = (gene_summary["cds_de"]).astype(int)
pred_labels = (gene_summary["assembler_de"]).astype(int)


cm = confusion_matrix(true_labels, pred_labels, labels=[0, 1])

linear_cm = {
    f"{row_prefix}DE_TN": cm[0, 0],
    f"{row_prefix}DE_FP": cm[0, 1],
    f"{row_prefix}DE_FN": cm[1, 0],
    f"{row_prefix}DE_TP": cm[1, 1],
}

pd.DataFrame.from_dict(linear_cm, orient="index", columns=[prefix]).to_csv(
    f"{row_prefix}{prefix}_linear_cm.tsv", sep="\t"
)
