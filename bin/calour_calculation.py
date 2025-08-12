import calour as ca
import numpy as np
import sys

ca.set_log_level(11)
np.random.seed(2018)

matrix = sys.argv[1]
min_abundance = int(sys.argv[2])
prefix = sys.argv[3]

exp = ca.read(data_file=matrix, data_file_type="tsv", normalize=None)
exp = exp.filter_sum_abundance(min_abundance)
exp.sample_metadata['Group'] = exp.sample_metadata['_sample_id'].apply(lambda x: [part for part in x.split('_') if 'group' in part][0])
res = exp.diff_abundance("Group", "group1", "group2")

res.feature_metadata.to_csv(f"{prefix}_calour_full_table.tsv", sep='\t', index=True)
