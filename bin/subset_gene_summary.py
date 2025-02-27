import pandas as pd
import sys

summary_file = sys.argv[1]
prefix = sys.argv[2]

gene_summary = pd.read_csv(summary_file, index_col=0)

gene_counts = gene_summary[[col for col in list(gene_summary.columns) if col.startswith("sample")]]

gene_counts.to_csv(f"{prefix}_count.tsv", sep="\t", index=True)
