#!/usr/bin/env python

import pandas as pd
import sys
import json
from Bio import SeqIO


def deduplicate_by_key(df, group_key, seed):
    if df.empty:
        return df
    minimap_mapping = df.copy()
    group = [group_key]

    max_cols = ["Mapping quality", "Number of matching bases in the mapping"]
    new_max_cols = ["Max Mapping quality", "Max Number of matching bases in the mapping"]

    minimap_mapping[new_max_cols] = (minimap_mapping.groupby(group)[max_cols].transform("max"))

    min_col = ["Query sequence length"]
    new_min_col = ["Min Query sequence length"]

    minimap_mapping[new_min_col] = minimap_mapping.groupby(group)[min_col].transform("min")

    # calculate a cumulative score to include hits with less then optimal matches
    minimap_mapping["cumulative_score"] = 5 * (minimap_mapping[max_cols[0]] == minimap_mapping[new_max_cols[0]]) + 3 * (minimap_mapping[max_cols[1]] == minimap_mapping[new_max_cols[1]]) + 1* (minimap_mapping[min_col[0]] == minimap_mapping[new_min_col[0]])

    minimap_mapping["max_cumulative_score"] = minimap_mapping.groupby(group)["cumulative_score"].transform("max")

    minimap_mapping["optimal_hit"] = (minimap_mapping["cumulative_score"] == minimap_mapping["max_cumulative_score"])

    cooptimal_bools = (minimap_mapping[minimap_mapping["optimal_hit"]].groupby(group).transform("size") > 1)

    if sum(cooptimal_bools):
        minimap_mapping.loc[cooptimal_bools.index, "cooptimal_hits"] = cooptimal_bools

        minimap_mapping["cooptimal_hits"] = (
                minimap_mapping["cooptimal_hits"]
                .astype("boolean")
                .fillna(False)
                .astype(bool) 
            )
        cooptimal_hits_selection_idx = minimap_mapping[minimap_mapping["cooptimal_hits"]].groupby(group).sample(n=1, random_state=seed).index

        minimap_mapping.loc[cooptimal_hits_selection_idx, "selected_cooptimal_hit"] = True

        minimap_mapping["selected_cooptimal_hit"] = (
                minimap_mapping["selected_cooptimal_hit"]
                .astype("boolean")
                .fillna(False)
                .astype(bool)
            )
        selection_idx = minimap_mapping[minimap_mapping["optimal_hit"] & (~minimap_mapping["cooptimal_hits"] | minimap_mapping["selected_cooptimal_hit"])].index
        return df.loc[selection_idx]
    
    else:
        optimal_selection = minimap_mapping[minimap_mapping["optimal_hit"]].index 
        return df.loc[optimal_selection]



minimap_path = sys.argv[1]
seed = int(sys.argv[2])
prefix = sys.argv[3]

mapping_col_names = [
    "contig",  # Query sequence name
    "Query sequence length",
    "Query start coordinate",
    "Query end coordinate",
    "same strand",
    "block_id",  # "Target sequence name",
    "Target sequence length",
    "Target start coordinate on the original strand",
    "Target end coordinate on the original strand",
    "Number of matching bases in the mapping",
    "Number bases, including gaps, in the mapping",
    "Mapping quality", "1", "2", "3", "4", "5", "6"
]

json_cols = ["contig", "block_id"]

minimap_mapping = pd.read_csv(
        minimap_path, sep="\t", header=None, names=mapping_col_names, dtype={'contig': str}
)

minimap_mapping["correctly_mapped_bases"] = minimap_mapping["Number of matching bases in the mapping"] / minimap_mapping["Target sequence length"]
# -> "Number of matching bases in the mapping"/ "Target sequence length" (blocks are targets)
# Rationale: Number of matching bases/ length of block -> block bases should be reconstructed
minimap_mapping["covered_reference"] = minimap_mapping["Number bases, including gaps, in the mapping"] / minimap_mapping["Target sequence length"]

pre_filter_mappings = minimap_mapping.shape[0]

minimap_mapping.to_csv("ratioed_filtered_df.tsv", sep="\t")

minimap_mapping = minimap_mapping[minimap_mapping["correctly_mapped_bases"] >= 0.95]
minimap_mapping = minimap_mapping[
    (minimap_mapping["covered_reference"] <= 1.05) &
    (minimap_mapping["covered_reference"] >= 0.95)
]

minimap_mapping.to_csv("filtered_df.tsv", sep="\t")

minimap_mapping = minimap_mapping.reset_index(drop=True)

minimap_mapping["block_size"] = (
    minimap_mapping
    .groupby(["block_id"])
    .transform("size")
)

minimap_mapping["contig_size"] = (
    minimap_mapping
    .groupby(["contig"])
    .transform("size")
)

minimap_mapping = minimap_mapping.sort_values(
    by=["Mapping quality", "Number of matching bases in the mapping", "Query sequence length"], ascending=[False, False, True]
)

minimap_mapping_chosen = minimap_mapping.copy()

minimap_mapping_chosen = deduplicate_by_key(minimap_mapping_chosen, "block_id", seed)

minimap_mapping_chosen = deduplicate_by_key(minimap_mapping_chosen, "contig", seed)

minimap_mapping_chosen["category"] = "overlap_block"

minimap_mapping_chosen[["contig", "block_id", "category"]].to_csv(f"{prefix}_mapped_overlap_blocks.tsv", sep="\t")