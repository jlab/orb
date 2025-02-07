import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

prefix = sys.argv[1]


files = sys.argv[2:]

dfs = []

for file in files:
    dfs.append(pd.read_csv(file, sep="\t"))

first_columns = [df.iloc[:, 0] for df in dfs]
merged_df = pd.concat(first_columns, axis=1)

global_min = merged_df.min().min()
global_max = merged_df.max().max()

axes = merged_df.hist(bins=100, figsize=(12, 8))

for ax in axes.flatten():
    ax.set_xlim(global_min, global_max)

plt.tight_layout()
plt.savefig(f"{prefix}_histograms.png", dpi=900)
plt.close() 

global_min = merged_df.min().min()
global_max = merged_df.max().max()

bin_count = 100

bin_edges = np.linspace(global_min, global_max, bin_count + 1)

global_y_max = 0
for col in merged_df.columns:
    counts, _ = np.histogram(merged_df[col].dropna(), bins=bin_edges)
    global_y_max = max(global_y_max, counts.max())

axes = merged_df.hist(bins=bin_edges, figsize=(12, 8), log=True)

for ax in axes.flatten():
    ax.set_xlim(global_min, global_max)
    ax.set_ylim(0, global_y_max)

plt.tight_layout()
plt.savefig(f"{prefix}_log_histograms.png", dpi=900)
plt.close()

