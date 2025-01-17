import pandas
from Bio import SeqIO
import sys

bed_file = sys.argv[1]
ref_transcriptome = sys.argv[2]

bed_df = pandas.read_csv(bed_file, sep="\t", header=None)
bed_df = bed_df[bed_df[0] == bed_df[3].str.split("_", expand=True)[0] + "_" + bed_df[3].str.split("_", expand=True)[1]]

blocks = []

for group_name, group in bed_df.groupby(0):
    sorted_group = group.sort_values(1)
    new_cds = True
    block_index = 0
    for _, row in sorted_group[[1, 2]].iterrows():
        if new_cds:
            current_start = row[1]
            current_end = row[2]
            new_cds = False
        else:
            if row[1] <= current_end:
                current_end = row[2]
            else:
                blocks.append((group_name, block_index, current_start, current_end))
                current_start = row[1]
                current_end = row[2]
                block_index += 1
    blocks.append((group_name, block_index, current_start, current_end))

blocks_df = pandas.DataFrame(blocks, columns=["cds", "block_index", "start", "end"])

cds_to_block_dict = blocks_df.groupby("cds").apply(lambda group: group[["block_index", "start", "end"]].values.tolist(), include_groups=False).to_dict()

in_data_set = 0
not_in_data_set = 0
all_cds = blocks_df["cds"].unique()
block_records = []

for record in SeqIO.parse(ref_transcriptome, "fasta"):
    if record.id in cds_to_block_dict:
        for block in cds_to_block_dict[record.id]:
            block_seq = record.seq[block[1]:block[2]]
            block_name = f"{record.id}_block{block[0]}"
            block_records.append(SeqIO.SeqRecord(block_seq, block_name, description=""))

SeqIO.write(block_records, sys.stdout, "fasta")
