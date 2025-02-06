import pandas as pd
import sys

score_over_view = sys.argv[1]
blocks = sys.argv[2]
chimeric_blocks = sys.argv[3]

scores_df = pd.read_csv(score_over_view, sep='\t', index_col=0)
blocks_df = pd.read_csv(blocks, sep='\t', index_col=0)
chimeric_blocks_df = pd.read_csv(chimeric_blocks, sep='\t', index_col=0)

n_blocks = blocks_df.shape[0]
n_chimeric_blocks = chimeric_blocks_df.shape[0]

scores_df.loc["mapped_blocks_ratio"] = scores_df.loc["mapped_contigs"] / n_blocks
scores_df.loc["mapped_chim_blocks_ratio"] = scores_df.loc["chimeric_mapped_contigs"] / n_chimeric_blocks
scores_df.loc["mapped_contigs_ratio"] = (scores_df.loc["mapped_contigs"] + scores_df.loc["chimeric_mapped_contigs"])/ scores_df.loc["total_contigs"]
scores_df.loc["mapped_f1_score"] = (2 * (scores_df.loc["mapped_blocks_ratio"] + scores_df.loc["mapped_chim_blocks_ratio"]) * scores_df.loc["mapped_contigs_ratio"]) / (scores_df.loc["mapped_blocks_ratio"] + scores_df.loc["mapped_chim_blocks_ratio"] + scores_df.loc["mapped_contigs_ratio"])

scores_df.to_csv(sys.stdout, sep='\t')