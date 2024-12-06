import sys
import pandas as pd


def merge_dfs(df_files):
    first = True
    for df_file in df_files:
        if first:
            dfs = pd.read_csv(df_file, sep="\t", index_col=0)
            first = False
        else:
            df = pd.read_csv(df_file, sep="\t", index_col=0)
            dfs = pd.merge(dfs, df, left_index=True, right_index=True)
    dfs.to_csv(sys.stdout, sep="\t", index=True)


len_args = len(sys.argv)
df_files = sys.argv[1:len_args]
merge_dfs(df_files)
