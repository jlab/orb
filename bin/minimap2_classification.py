#!/usr/bin/env python

import pandas as pd
import sys
import json
from Bio import SeqIO


def deduplicate_by_key(df, group_key, seed):
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


def save_dict(filtered_df, output_file_name, cols):
    key_val_dict = filtered_df[json_cols].set_index(json_cols[0])[json_cols[1]].to_dict()
    with open(f"{output_file_name}.json", "w") as f:
        json.dump(key_val_dict, f)


def count_ns_in_contigs(contigs_path):
    all_entries = []

    for record in SeqIO.parse(contigs_path, "fasta"):
        total_ns = record.seq.count("N")
        max_run = 0
        current = 0
        for c in record.seq:
            if c == "N":
                current += 1
                max_run = max(max_run, current)
            else:
                current = 0
        all_entries.append([record.id, total_ns, max_run])

    return pd.DataFrame(all_entries, columns=["contig_id", "total_Ns", "max_N_run"])


minimap_path = sys.argv[1]
contig_path = sys.argv[2]
gene_summary_path = sys.argv[3]
seed = int(sys.argv[4])
prefix = sys.argv[5]

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

gene_summary = pd.read_csv(gene_summary_path)
gene_to_og = gene_summary.set_index("gene_name")["orthogroup"].to_dict()

ns_per_contig_df = count_ns_in_contigs(contig_path)

minimap_mapping["correctly_mapped_bases"] = minimap_mapping["Number of matching bases in the mapping"] / minimap_mapping["Target sequence length"]
# -> "Number of matching bases in the mapping"/ "Target sequence length" (blocks are targets)
# Rationale: Number of matching bases/ length of block -> block bases should be reconstructed
minimap_mapping["covered_reference"] = minimap_mapping["Number bases, including gaps, in the mapping"] / minimap_mapping["Target sequence length"]

pre_filter_mappings = minimap_mapping.shape[0]

minimap_mapping = minimap_mapping[minimap_mapping["correctly_mapped_bases"] >= 0.95]
minimap_mapping = minimap_mapping[
    (minimap_mapping["covered_reference"] <= 1.05) &
    (minimap_mapping["covered_reference"] >= 0.95)
]

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

# add random seeds for cooptimal results

minimap_mapping.loc[minimap_mapping_chosen.index, "assigned"] = True

minimap_mapping["assigned"] = (
    minimap_mapping["assigned"]
    .astype("boolean")
    .fillna(False)
    .astype(bool)
)

minimap_mapping["is_duplicated"] = (~(minimap_mapping["assigned"])) & (minimap_mapping["block_size"] > 1)

minimap_mapping["is_joined"] = (minimap_mapping["contig_size"] > 1)

minimap_mapping["blocks_aggregated"] = (
    minimap_mapping
    .groupby("contig")["block_id"]
    .transform(lambda x: [list(x)] * len(x))
)

minimap_mapping["gene"] = minimap_mapping["block_id"].str.replace(r"^((?:[^_]+_){1}[^_]+).*", r"\1", regex=True)

minimap_mapping["genes_aggregated"] = (
    minimap_mapping
    .groupby("contig")["gene"]
    .transform(lambda x: [list(x)] * len(x))
)

minimap_mapping["is_one_gene"] = minimap_mapping["genes_aggregated"].apply(lambda x: len(set(x)) == 1)

minimap_mapping["og"] = minimap_mapping["gene"].map(gene_to_og)

minimap_mapping["ogs_aggregated"] = (
    minimap_mapping
    .groupby("contig")["og"]
    .transform(lambda x: [list(x)] * len(x))
)

minimap_mapping["is_one_og"] = minimap_mapping["ogs_aggregated"].apply(lambda x: len(set(x)) == 1)

minimap_mapping["is_recovered"] = (minimap_mapping["assigned"] & (~minimap_mapping["is_joined"]))

minimap_mapping["is_block_recovered"] = minimap_mapping.groupby("block_id")["is_recovered"].transform("max")

block_is_super_recovered = minimap_mapping[minimap_mapping["assigned"] & (minimap_mapping["is_joined"]) & (minimap_mapping["is_one_gene"])]["blocks_aggregated"].explode()

minimap_mapping["is_block_super_recovered"] = minimap_mapping["block_id"].isin(block_is_super_recovered)

minimap_mapping["is_block_super_recovered_or_recovered"] = minimap_mapping["is_block_recovered"] | minimap_mapping["is_block_super_recovered"]

# calculate blocks but subtract the blocks already assigned

super_recovered_amount_of_blocks = minimap_mapping[minimap_mapping["assigned"] & (minimap_mapping["is_joined"]) & (minimap_mapping["is_one_gene"])]["contig_size"].sum() - (minimap_mapping["assigned"] & (minimap_mapping["is_joined"]) & ~(minimap_mapping["is_one_gene"]) & (minimap_mapping["is_one_og"]) & minimap_mapping["is_block_recovered"]).sum()

og_recovered_amount_of_blocks = minimap_mapping[minimap_mapping["assigned"] & (minimap_mapping["is_joined"]) & ~(minimap_mapping["is_one_gene"]) & (minimap_mapping["is_one_og"])]["contig_size"].sum() - (minimap_mapping["assigned"] & (minimap_mapping["is_joined"]) & ~(minimap_mapping["is_one_gene"]) & minimap_mapping["is_one_og"] & minimap_mapping["is_block_super_recovered_or_recovered"]).sum()

dedupped_rest = minimap_mapping[~minimap_mapping["contig"].isin(minimap_mapping[minimap_mapping["assigned"]]["contig"])].groupby("contig").head(1)

duplicated_non_chimeric_sum = (dedupped_rest["is_duplicated"] & ~dedupped_rest["is_joined"]).sum()
duplicated_chimeric_sum = (dedupped_rest["is_duplicated"] & dedupped_rest["is_joined"]).sum()

recovered = (
    minimap_mapping[(minimap_mapping["assigned"] & (~minimap_mapping["is_joined"])) | (minimap_mapping["assigned"] & (minimap_mapping["is_joined"]) & (minimap_mapping["is_one_gene"]))]
)

# add final contig classification

og_recovered = (
    minimap_mapping[(minimap_mapping["assigned"] & (minimap_mapping["is_joined"]) & ~(minimap_mapping["is_one_gene"]) & (minimap_mapping["is_one_og"]))]
)

chimeric = (
    minimap_mapping[minimap_mapping["assigned"] & (minimap_mapping["is_joined"]) & ~(minimap_mapping["is_one_gene"]) & ~(minimap_mapping["is_one_og"])]
)

duplicated_non_chimeric = (
    dedupped_rest[(dedupped_rest["is_duplicated"] & ~dedupped_rest["is_joined"])]
)

duplicated_chimeric = (
    dedupped_rest[(dedupped_rest["is_duplicated"] & dedupped_rest["is_joined"])]
)

minimap_mapping.loc[
    (minimap_mapping["assigned"] & (~minimap_mapping["is_joined"])), "category"
] = "minimap2_single_recovered" 

minimap_mapping.loc[
    (minimap_mapping["assigned"] & (minimap_mapping["is_joined"]) & (minimap_mapping["is_one_gene"])), "category"
] = "minimap2_merged_recovered" 

#super_recovered_amount_of_blocks = minimap_mapping[
# minimap_mapping["assigned"] & (minimap_mapping["is_joined"]) & (minimap_mapping["is_one_gene"])]["contig_size"].sum() - 
# (minimap_mapping["assigned"] & (minimap_mapping["is_joined"]) & ~(minimap_mapping["is_one_gene"]) & (minimap_mapping["is_one_og"]) & minimap_mapping["is_block_recovered"]).sum()


minimap_mapping.loc[
    (minimap_mapping["assigned"] & (minimap_mapping["is_joined"]) & ~(minimap_mapping["is_one_gene"]) & (minimap_mapping["is_one_og"])), "category"
] = "minimap2_orthologous_recovered" 

minimap_mapping.loc[
    (minimap_mapping["assigned"] & (minimap_mapping["is_joined"]) & ~(minimap_mapping["is_one_gene"]) & ~(minimap_mapping["is_one_og"])), "category"
] = "minimap2_chimeric"

minimap_mapping.loc[
    (dedupped_rest[dedupped_rest["is_duplicated"] & ~dedupped_rest["is_joined"]]).index, "category"
] = "minimap2_duplicated_non_chimeric"

minimap_mapping.loc[
    (dedupped_rest[dedupped_rest["is_duplicated"] & dedupped_rest["is_joined"]]).index, "category"
] = "minimap2_duplicated_chimeric"

minimap_mapping.to_csv("debug.tsv", sep="\t")


minimap_mapping.dropna(
    subset=["category"], inplace=True
)

# calculate how many Ns are in the assigned contigs
contigs_n_df = pd.merge(
    minimap_mapping,
    ns_per_contig_df,
    how='left',
    right_on='contig_id',
    left_on='contig'
)

recovered_with_n = (minimap_mapping["category"].isin(["minimap2_single_recovered", "minimap2_merged_recovered", "minimap2_orthologous_recovered"]) & (contigs_n_df["total_Ns"] > 0)).sum()

chimeric_with_n = ((minimap_mapping["category"] == "minimap2_chimeric") & (contigs_n_df["total_Ns"] > 0)).sum()

minimap_mapping[["contig", "block_id", "category", "blocks_aggregated"]].to_csv(f"{prefix}_minimap2_categories.tsv", sep="\t")

save_dict(recovered, f"{prefix}_recovered", json_cols)

categories = minimap_mapping["category"].value_counts()
scores = pd.Series([super_recovered_amount_of_blocks, og_recovered_amount_of_blocks, recovered_with_n, chimeric_with_n], index=["minimap2_merged_recovered_blocks", "minimap2_orthologous_recovered_blocks", "minimap2_recovered_with_n", "minimap2_chimeric_with_n"])

categories.add(scores, fill_value=0).sort_index().to_csv(f"{prefix}_minimap2_category_counts.tsv", sep="\t", header=[prefix])

ns_per_contig_df.to_csv(f"{prefix}_contig_n_stats.tsv", sep="\t", index=False)
