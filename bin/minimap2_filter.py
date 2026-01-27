#!/usr/bin/env python

import pandas as pd
import sys
import os
import json

# maximum matching bases
# if its the same which is unlikely randomly select one

# Minimap2 output format
# https://lh3.github.io/minimap2/minimap2.html
# Matching bases in col 10 (9 with index 0)

file_path = sys.argv[1] # "/home/tlin/Projects/bmp_analysis/data/refbased/velvet_mapping.tsv"

minimap_file_basename = os.path.basename(file_path)
minimap_file_basename = os.path.splitext(minimap_file_basename)[0]


if os.path.exists(file_path) and os.path.getsize(file_path) == 0:
    with open(f"{minimap_file_basename}_dedup.json", 'w') as json_file:
        json.dump({}, json_file)

minimap_results = pd.read_csv(file_path, sep="\t", header=None, usecols=[0, 1, 5, 6, 9, 10, 11])

minimap_results["correctly_mapped_bases"] = minimap_results[9] / minimap_results[6]
# -> "Number of matching bases in the mapping"/ "Target sequence length" (blocks are targets)
# Rationale: Number of matching bases/ length of block -> block bases should be reconstructed
minimap_results["covered_reference"] = minimap_results[10] / minimap_results[6]
# -> "Number bases, including gaps, in the mapping"/ "Target sequence length" (blocks are targets)
# Rationale: Length of alignment/ length of block -> length of block should be covered and not too many gaps

minimap_results = minimap_results[minimap_results["correctly_mapped_bases"] >= 0.95]
minimap_results = minimap_results[
    (minimap_results["covered_reference"] <= 1.05) &
    (minimap_results["covered_reference"] >= 0.95)
]

# Select longest alignment by: Number of matching bases in the mapping

minimap_filtered = minimap_results.sort_values(by=[11, 9, 1], ascending=[False, False, True]  # sort hits by mapping quality AND number of involved matching bases AND query length (descending)
            ).groupby([5], sort=False).head(1).groupby([0], sort=False).head(1) #first choose the most suitable contig for each block, then find the most suitable block for each contig

if minimap_filtered[0].duplicated().any():
    raise ValueError("Duplicate keys")

key_val_dict = minimap_filtered.set_index(0)[5].to_dict()

with open(f"{minimap_file_basename}_dedup.json", "w") as f:
    json.dump(key_val_dict, f)
