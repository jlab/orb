import sys
import pandas as pd

df_file = sys.argv[1]
axis = int(sys.argv[2])
prefix = sys.argv[3]

df = pd.read_csv(df_file, sep='\t', index_col=0)

df = df.sort_index(axis=axis)

df.to_csv(f"{prefix}_sorted.tsv", sep="\t", index=True)
