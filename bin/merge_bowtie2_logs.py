import sys
import pandas as pd


def merge_bowtie2_logs(log_files, column_name):
    first = True
    for log_file in log_files:
        if first:
            log_dfs = pd.read_csv(log_file, sep="\t", index_col=0)
            first = False
        else:
            log_df = pd.read_csv(log_file, sep="\t", index_col=0)
            log_dfs = pd.merge(log_dfs, log_df, left_index=True, right_index=True)

    merged_df = log_dfs.mean(axis=1)
    pd.DataFrame(merged_df, columns=[f"{column_name}_mean"]).to_csv(f"{column_name}_mean_logs.tsv", sep="\t", index=True)
    log_dfs.to_csv(f"{column_name}_merged_logs.tsv", sep="\t", index=True)


len_args = len(sys.argv)
log_files = sys.argv[1:len_args-1]
col_name = sys.argv[len_args-1]
merge_bowtie2_logs(log_files, col_name)
