#!/usr/bin/env python

"""
This script is for categorising the contigs into multiple classes. Will produce a tsv file with contig id and the corresponding class.

The dependency files are:
all_contig_ids_fl: all contig ids of assembler
mapped_ids_fl: the contigs mapped to a block
chimeric_mapped_ids_fl: the contigs mapped to an overlap block
length_filtered_ids_fl: contigs filtered by length
assembler_mapping_fl: the reads mapped onto the contigs
gene_sumamry: the Marbel gene summary file
prefix: assembler name

Result file will have format:

contigs\tcategory
...
"""

import sys
import re

import polars as pl

all_contig_ids_fl = sys.argv[
    1
]  # "/vol/jlab/tlin/no_backup/nextflow_workdir/00/40ed8ff06c8b8b7cc7997b34a4eb99/rnaspades_contigs_ids.txt"

minimap2_categories_fl = sys.argv[
    2
]  # "/vol/jlab/tlin/no_backup/nextflow_workdir/42/62113b3dfc94ec10187ea81edf48a0/rnaspades_values.txt"

length_filtered_ids_fl = sys.argv[
    3
]  # "/vol/jlab/tlin/no_backup/nextflow_workdir/00/40ed8ff06c8b8b7cc7997b34a4eb99/rnaspades_length_filtered_ids.txt"

assembler_mapping_fl = sys.argv[
    4
]  # "/vol/jlab/tlin/no_backup/nextflow_workdir/3c/f92325783c8fceb5ac01245489cc77/rnaspades_1000000.tsv"

gene_summary_path = sys.argv[
    5
]  # "/vol/jlab/tlin/marbel_benchmarking_integration/benchmarking_sets_all_sparse_fixed_libsize/moss_microbiome/summary/gene_summary.csv"

prefix = sys.argv[6]  # "rnaspades"

all_contigs_ids = pl.read_csv(all_contig_ids_fl, has_header=False, new_columns=["contigs"])

try:
    minimap2_categories = pl.read_csv(minimap2_categories_fl, separator="\t")
    mapped_ids = minimap2_categories["contig"].to_list()
except (pl.exceptions.NoDataError, OSError):
    mapped_ids = []

try:
    length_filtered_ids = pl.read_csv(length_filtered_ids_fl, has_header=False, new_columns=["l_filtered"])[
        "l_filtered"
    ].to_list()
except (pl.exceptions.NoDataError, OSError):
    length_filtered_ids = []


assembler_mapping = pl.read_csv(
    assembler_mapping_fl, separator="\t", has_header=False, new_columns=["contigs", "start", "end", "read_id"]
)

gene_summary = pl.read_csv(gene_summary_path).select(["gene_name", "orthogroup"])

gene_og_dict = dict(zip(gene_summary["gene_name"], gene_summary["orthogroup"]))

# Get original gene name and orthogroup from read name
assembler_mapping = assembler_mapping.with_columns(
    [pl.col("read_id").str.extract(r"^([^_]+_[^_]+)", 1).alias("origin_gene")]
)

# Map gene_og_dict
assembler_mapping = assembler_mapping.with_columns(
    [
        pl.col("origin_gene")
        .map_elements(lambda x: gene_og_dict.get(x, None), return_dtype=pl.String)
        .alias("origin_orthogroup")
    ]
)

contig_aggregation = assembler_mapping.group_by("contigs").agg(
    [
        pl.col("origin_gene").n_unique().alias("origin_gene_nunique"),
        pl.col("origin_orthogroup").n_unique().alias("origin_orthogroup_nunique"),
    ]
)

single_mapped_contigs = contig_aggregation.filter(pl.col("origin_gene_nunique") == 1)["contigs"].to_list()

# Step 1: Filter for contigs with multiple origin genes
multi_mapped = contig_aggregation.filter(pl.col("origin_gene_nunique") != 1)

# multi mapped with multiple ogs
multi_mapped_contigs_multi_og = multi_mapped.filter(pl.col("origin_orthogroup_nunique") != 1)["contigs"].to_list()

# multimapped but single og
multi_mapped_contigs_single_og = multi_mapped.filter(pl.col("origin_orthogroup_nunique") == 1)["contigs"].to_list()

unassigned_contigs = all_contigs_ids.filter(~pl.col("contigs").is_in(mapped_ids))

assembler_contigs = assembler_mapping["contigs"].to_list()

unmapped_contigs = unassigned_contigs.filter(~pl.col("contigs").is_in(assembler_contigs))["contigs"].to_list()

contigs_with_cat_no_mapped = unassigned_contigs.with_columns(
    [
        pl.when(pl.col("contigs").is_in(length_filtered_ids))
        .then(pl.lit("length_filtered_contigs"))
        .when(pl.col("contigs").is_in(unmapped_contigs))
        .then(pl.lit("unmapped_contigs"))
        .when(pl.col("contigs").is_in(single_mapped_contigs))
        .then(pl.lit("single_mapped_contigs"))
        .when(pl.col("contigs").is_in(multi_mapped_contigs_multi_og))
        .then(pl.lit("multi_mapped_contigs_multi_og"))
        .when(pl.col("contigs").is_in(multi_mapped_contigs_single_og))
        .then(pl.lit("multi_mapped_contigs_single_og"))
        .otherwise(pl.lit("no_category"))
        .alias("category")
    ]
)

print("before concat")

result = pl.concat(
    [contigs_with_cat_no_mapped, minimap2_categories.select(["contig", "category"]).rename({"contig": "contigs"})],
    how="vertical",
)

result.write_csv(f"{prefix}_contigs_categorised.tsv", separator="\t")
