import edlib
from itertools import combinations
import sys
from Bio import SeqIO
import json
import os
import pandas as pd
from tqdm import tqdm


def edlib_identity(seq1, seq2):
    if seq1 == seq2:
        return 100.0
    result = edlib.align(seq1, seq2, mode="NW", task="path")
    edit_dist = result["editDistance"]
    return (1 - edit_dist / max(len(seq1), len(seq2))) * 100


def compute_edlib_pairwise_identities(genes, sequence_db):
    seqs = {g: sequence_db[g].seq for g in genes if g in sequence_db}
    if len(seqs) != len(genes):
        missing_genes = set(genes) - set(seqs.keys())
        print(f"Warning: Some sequences not found in the database for genes: {missing_genes}")
    total = 0
    count = 0
    keys = list(seqs.keys())
    for g1, g2 in combinations(keys, 2):
        identity = edlib_identity(seqs[g1], seqs[g2])
        total += identity
        count += 1
    mean_identity = total / count if count else 100.0
    return mean_identity


def calculate_og_seq_ident(summary_df, sequences_file, env_name, clear_cache=False):
    summary = pd.read_csv(summary_df)
    sequence_db = SeqIO.index_db(sequences_file)
    mean_ids = []
    og_names = []

    if not os.path.exists("data_cache"):
        os.mkdir("data_cache")

    og_identity_file = f"{env_name}_og_mean_ident_values.csv"

    if os.path.exists(f"data_cache/{og_identity_file}"):
        return pd.read_csv(f"data_cache/{og_identity_file}")

    for og, genes in tqdm(summary.groupby("orthogroup")["gene_name"]):
        mean_id = compute_edlib_pairwise_identities(genes, sequence_db)
        og_names.append(og)
        mean_ids.append(mean_id)

    df = pd.DataFrame(
        {
            "og_name": og_names,
            "mean_identity": mean_ids,
        }
    )

    df.to_csv(f"data_cache/{og_identity_file}", index=False)
    return df
