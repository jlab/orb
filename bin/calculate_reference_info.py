import pandas as pd
from Bio import SeqIO
import sys

blocks_file = sys.argv[1]
chim_block_file = sys.argv[2]
reference_fasta = sys.argv[3]
prefix = sys.argv[4]

blocks = pd.read_csv(blocks_file, sep="\t", header=None)
blocks.columns = ["cds", "start", "end", "block_name", "fragment_count"]
chimeric_blocks = pd.read_csv(chim_block_file, sep="\t")

blocks["length"] = blocks["end"] - blocks["start"]
chimeric_blocks["length"] = chimeric_blocks["overlap_block_end"] - chimeric_blocks["overlap_block_start"]

total_block_length = blocks["length"].sum()
total_chimeric_block_length = chimeric_blocks["length"].sum()

cds_lengths = []


for seq in SeqIO.parse(reference_fasta, "fasta"):
    cds_lengths.append(len(seq.seq))

total_cds_length = sum(cds_lengths)

pd.DataFrame({"Block Lengths": blocks["length"]}).to_csv(f"{prefix}_block_lengths.tsv", sep="\t", index=False)
pd.DataFrame({"Chimeric Blocks Lengths":chimeric_blocks["length"]}).to_csv(f"{prefix}_chimeric_block_lengths.tsv", sep="\t", index=False)
pd.DataFrame({"CDS Lengths": cds_lengths}).to_csv(f"{prefix}_cds_lengths.tsv", sep="\t", index=False)

cds_summary = pd.concat([pd.Series(cds_lengths).describe(), pd.Series({"total_length": total_cds_length})])
blocks_summary = pd.concat([blocks["length"].describe(), pd.Series({"total_length": total_block_length})])
chim_blocks_summary = pd.concat([chimeric_blocks["length"].describe(), pd.Series({"total_length": total_chimeric_block_length})])

pd.DataFrame({'CDS': cds_summary, 'Blocks': blocks_summary, 'Chimeric Blocks': chim_blocks_summary}).to_csv(f"{prefix}_reference_stats.tsv", sep="\t")
