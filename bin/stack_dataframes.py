import sys
import pandas as pd


def merge_dfs(df_files):
    first = True
    for df_file in df_files:
        if first:
            dfs = pd.read_csv(df_file, sep="\t", index_col=False, header=True)
            first = False
        else:
            df = pd.read_csv(df_file, sep="\t", index_col=False, header=True)
            dfs = pd.concat([dfs, df])
    dfs.to_csv(sys.stdout, sep="\t", index=False, header=True)


len_args = len(sys.argv)
df_files = sys.argv[1:len_args]
merge_dfs(df_files)
