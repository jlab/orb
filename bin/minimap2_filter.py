import pandas as pd
import sys
import os
import json

# maximum matching bases
# if its the same which is unlikely randomly select one

# Minimap2 output format
# https://lh3.github.io/minimap2/minimap2.html
# Matching bases in col 10 (9 with index 0)

file_path = sys.argv[1]
minimap_file_basename = os.path.basename(file_path)
minimap_file_basename = os.path.splitext(minimap_file_basename)[0]

if os.path.exists(file_path) and os.path.getsize(file_path) == 0:
    with open(f"{minimap_file_basename}_dedup.json", 'w') as json_file:
        json.dump({}, json_file)

minimap_results = pd.read_csv(file_path, sep="\t", header=None)

indices = pd.Index([])

for cds in minimap_results[5].unique():
    condition = (minimap_results[5] == cds)
    max_value = max(minimap_results[condition][9])
    new_indices = minimap_results.loc[condition & (minimap_results[9] < max_value)].index
    if not new_indices.empty:
        indices = indices.append(new_indices)
    max_vals = (condition) & (minimap_results[9] == max_value)
    logic_vals_sum = max_vals.sum()
    if logic_vals_sum > 1:
        indices = indices.append(minimap_results[max_vals].sample(n=logic_vals_sum - 1).index)

minimap_filtered = minimap_results.drop(indices)
minimap_filtered["correctly_mapped_bases"] = minimap_filtered[9] / minimap_filtered[6]
minimap_filtered["covered_reference"] = minimap_filtered[10] / minimap_filtered[6]

minimap_filtered = minimap_filtered[minimap_filtered["correctly_mapped_bases"] >= 0.95]
minimap_filtered = minimap_filtered[
    (minimap_filtered["covered_reference"] <= 1.05) |
    (minimap_filtered["covered_reference"] >= 0.95)
]

key_val_dict = dict(zip(minimap_filtered[0], minimap_filtered[5]))

with open(f"{minimap_file_basename}_dedup.json", "w") as f:
    json.dump(key_val_dict, f)