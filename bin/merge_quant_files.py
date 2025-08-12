import sys
import pandas as pd
from pathlib import Path


def merge_dfs(quant_files):
    first = True
    for quant_file in quant_files:
        quant_file = quant_file + "/quant.sf"
        if first:
            dfs = pd.read_csv(quant_file, sep="\t", index_col=0, usecols=["Name", "NumReads"])
            sample_name = Path(quant_file).parent
            dfs.rename(columns={"NumReads": sample_name}, inplace=True)
            first = False
        else:
            df = pd.read_csv(quant_file, sep="\t", index_col=0, usecols=["Name", "NumReads"])
            sample_name = Path(quant_file).parent
            df.rename(columns={"NumReads": sample_name}, inplace=True)
            dfs = pd.merge(dfs, df, left_index=True, right_index=True)
    dfs.to_csv(sys.stdout, sep="\t", index=True)


len_args = len(sys.argv)
quant_files = sys.argv[1:len_args]
merge_dfs(quant_files)
