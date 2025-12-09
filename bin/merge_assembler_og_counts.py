#!/usr/bin/env python

import polars as pl
import json
import re
import sys


assembler_count_dir_path = sys.argv[1] # mergequantsffiles
assembler_mapping_dir = sys.argv[2] # minimap2filter #{assembler}_mapping_dedup.json"
gene_summary_file = sys.argv[3]
assembler = sys.argv[4] # "rnaspades"

assembler_counts = pl.read_csv(assembler_count_dir_path, separator="\t")

assembler_mapping_dict = json.load(open(assembler_mapping_dir))

ref_genes = pl.read_csv(gene_summary_file)

# create a dictionary for mapping gene -> og
gene_to_og = dict(zip(ref_genes["gene_name"], ref_genes["orthogroup"]))

assembler_samples = [col for col in assembler_counts.columns if re.search(r"_sample_", col)]

# now i need to map the assemblers contig -> gene than gene -> og than aggregate 
assembler_with_og = assembler_counts.with_columns(
    pl.col("Name").cast(pl.Utf8),
).with_columns(
    pl.col("Name").map_elements(lambda x: assembler_mapping_dict.get(str(x), None), return_dtype=pl.Utf8).alias("gene_block_name"),
).with_columns(
    pl.col("gene_block_name").str.replace(r"_block\d+$", "", literal=False).alias("gene_name")
).with_columns(
    pl.col("gene_name").map_elements(lambda x: gene_to_og.get(str(x), None), return_dtype=pl.Utf8).alias("orthogroups")
).group_by("orthogroups").agg((pl.col(assembler_samples).sum()))

# prepare count file for dge tools
assembler_with_og.rename({"orthogroups": "Name"})
assembler_count_file = f"{assembler}_merged_orthogroups.tsv"
assembler_with_og.write_csv(assembler_count_file, separator="\t")


