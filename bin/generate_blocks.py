import pandas as pd
from Bio import SeqIO
import sys

bed_file = sys.argv[1]
ref_transcriptome = sys.argv[2]
dup_ids_group = sys.argv[3]
prefix = sys.argv[4]

dict_val = {}
for line in open(dup_ids_group):
    line = line.strip()
    ids = line.split("\t")[1]
    ids_list = ids.split(", ")
    dict_val[ids_list[0]] = ids_list

bed_df = pd.read_csv(bed_file, sep="\t", header=None)

#the read name contains the read numbered, so I need to truncate to see if it is the orginal cds
bed_df['truncated_red_name'] = bed_df[3].str.split('_', expand=True)[[0,1]].agg('_'.join, axis=1)
#we deduplicated the sequences, but we need to include the reads of sequences that were deduplicated, so map to all the original reads
bed_df["mapped_group"] = bed_df['truncated_red_name'].map(lambda x: dict_val.get(x, [x]))
bed_df = bed_df.explode("mapped_group")
#filter reads that are not the original cds, i.e. incorrectly mapped
bed_df = bed_df[bed_df[0]==bed_df["mapped_group"]]

blocks = []

for group_name, group in bed_df.groupby(0):
    #sort by start position
    sorted_group = group.sort_values(1)
    new_cds = True
    block_index = 0
    for _, row in sorted_group[[1, 2]].iterrows():
        if new_cds:
            current_start = row[1]
            current_end = row[2]
            fragment_count = 1
            new_cds = False
        else:
            fragment_count += 1
            if row[1] <= current_end:
                current_end = row[2]
            else:
                blocks.append((group_name, block_index, current_start, current_end, fragment_count))
                current_start = row[1]
                current_end = row[2]
                block_index += 1
                fragment_count = 1
    #add the last block
    blocks.append((group_name, block_index, current_start, current_end, fragment_count))

blocks_df = pd.DataFrame(blocks, columns=["cds", "block_index", "start", "end", "fragment_count"])
blocks_df.to_csv(f"{prefix}_blocks.tsv", sep="\t", index=False)

cds_to_block_dict = blocks_df.groupby("cds").apply(lambda group: group[["block_index", "start", "end"]].values.tolist(), include_groups=False).to_dict()

in_data_set = 0
not_in_data_set = 0
all_cds = blocks_df["cds"].unique()
block_records = []

for record in SeqIO.parse(ref_transcriptome, "fasta"):
    if record.id in cds_to_block_dict:
        for block in cds_to_block_dict[record.id]:
            block_seq = record.seq[block[1]:(block[2]+1)]
            block_name = f"{record.id}_block{block[0]}"
            block_records.append(SeqIO.SeqRecord(block_seq, block_name, description=""))

SeqIO.write(block_records, sys.stdout, "fasta")
