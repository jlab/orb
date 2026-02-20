#!/usr/bin/env python

import polars as pl
import re
import sys

gene_summary_file = sys.argv[1]
prefix = sys.argv[2]

# gene sumamry for og mapping
mt_genes = pl.read_csv(gene_summary_file)

# get sample cols
samples = [col for col in mt_genes.columns if re.search(r"_sample_", col)]

# aggregate by orthogroups
mt_ogs_aggregated = mt_genes.group_by("orthogroup").agg(
    (pl.col(samples).sum())
)
mt_ogs_aggregated = mt_ogs_aggregated.rename({"orthogroup": "Name"})

args = "group_1 10 group_2 10 \t".split(" ")

mt_og_count_file = f"{prefix}_merged_orthogroups.tsv"
mt_ogs_aggregated.write_csv(mt_og_count_file, separator="\t")
