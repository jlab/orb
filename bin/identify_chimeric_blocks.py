import pandas as pd
from Bio import SeqIO
import sys
import glob
import re


base_gtf_dir = sys.argv[1]
gene_summary = sys.argv[2]
translation_df = sys.argv[3]
all_blocks = sys.argv[4]
dup_ids_group = sys.argv[5]
prefix = sys.argv[6]

pd_gene = pd.read_csv(gene_summary)
pd_translation = pd.read_csv(translation_df, index_col=0)
pd_blocks = pd.read_csv(all_blocks, sep="\t")

translation_dict = pd_translation.set_index("species_long_name").to_dict()["accession"]

pd_gene["accession"] = pd_gene["origin_species"].apply(lambda x: translation_dict[x])

# for dedup mapping
group_dict = {}
for line in open(dup_ids_group):
    line = line.strip()
    ids = line.split("\t")[1]
    ids_list = ids.split(", ")
    group_dict[ids_list[0]] = ids_list

# in the gtf annotation the gene_id needs to extracted, use regex as it is more efficient and beautiful than string splits
pattern = r'gene_id\s+"([^"]+)"'
gene_to_genomic_location = {}


# i want to map the gene names to the genomic locations in order to infer chimeric blocks (i.e. blocks overlapping in the original genome)
for group_name, group in pd_gene.groupby("accession"):
    accession = re.sub(r"\.\d+$", "", group_name)
    ref_mapping = glob.glob(f"{base_gtf_dir}/{accession}.*/genomic.gtf")
    if len(ref_mapping) == 1:
        ref_mapping = ref_mapping[0]
    else:
        raise ValueError(f"Multiple or no genome files found for {accession}")

    pd_group_gtf = pd.read_csv(
        ref_mapping,
        sep="\t",
        header=None,
        skiprows=5,
        usecols=[0, 2, 3, 4, 8],
        names=["chromosome", "feature", "start", "end", "attributes"]
    )
    pd_group_gtf = pd_group_gtf[pd_group_gtf["feature"] == "CDS"]
    pd_group_gtf["gene_id"] = pd_group_gtf["attributes"].str.extract(pattern)

    # the same gene can occur in multiple locations, so it needs to be a list of lists
    gene_to_genomic_location = (
        pd_group_gtf.groupby("gene_id")[["chromosome", "start", "end"]]
        .apply(lambda x: list(zip(x["chromosome"], x["start"], x["end"])))
        .to_dict()
    )
    pd_gene.loc[group.index, "genomic_locations"] = group["gene_name"].map(gene_to_genomic_location).apply(lambda x: [[group_name] + list(i) for i in x] if isinstance(x, list) else pd.NA)

gene_location_dict = pd_gene.set_index("gene_name").to_dict()["genomic_locations"]

# we need to merge the genomic coordinates of genes in deduplicated groups
for k, v in group_dict.items():
    for id in v[1:]:
        if k not in gene_location_dict or id not in gene_location_dict:
            print(f"Warning: {k} or {id} not in gene_location_dict, make sure this is intended, skipping")
            continue
        gene_location_dict[k] += gene_location_dict.pop(id)


pd_blocks["genomic_location"] = pd_blocks["cds"].map(gene_location_dict)

pd_blocks_exploded = pd_blocks.explode("genomic_location")

# TODO: maybe remove the value error for publishing, others might have not a perfect datasource
if pd_blocks_exploded["genomic_location"].isna().sum() > 0:
    print(f"Warning: {pd_blocks_exploded['genomic_location'].isna().sum()} blocks have no genomic location")
    print("The dataframe:")
    print(pd_blocks_exploded[pd_blocks_exploded["genomic_location"].isna()])
    # raise ValueError("Some blocks have no genomic location, have a look at the dataframe")

pd_blocks_exploded = pd_blocks_exploded[~pd_blocks_exploded["genomic_location"].isna()]

pd_blocks_exploded['accession'], pd_blocks_exploded['chromosome'], pd_blocks_exploded['cds_genomic_start'], pd_blocks_exploded['cds_genomic_end'] = zip(*pd_blocks_exploded['genomic_location'])
pd_blocks_exploded.drop(columns=["genomic_location"], inplace=True)

pd_blocks_exploded["block_genomic_start"] = pd_blocks_exploded["start"] + pd_blocks_exploded["cds_genomic_start"]
pd_blocks_exploded["block_genomic_end"] = pd_blocks_exploded["end"] + pd_blocks_exploded["cds_genomic_start"]

chimeric_blocks = []
block_index = 0
for accession, blocks_per_accession in pd_blocks_exploded.groupby("accession"):
    for chromosome, blocks_per_chromosome in blocks_per_accession.groupby("chromosome"):
        sorted_group = blocks_per_chromosome.sort_values("block_genomic_start")
        new_chromosome = True
        for _, row in sorted_group.iterrows():
            if new_chromosome:
                current_start = row["block_genomic_start"]
                current_end = row["block_genomic_end"]
                origin_cds = [row["cds"]]
                origin_block_ids = [row["block_index"]]
                new_chromosome = False
            else:
                if row["block_genomic_start"] <= current_end:
                    current_end = row["block_genomic_end"]
                    origin_cds.append(row["cds"])
                    origin_block_ids.append(row["block_index"])
                else:
                    chimeric_blocks.append((f"chimblock_{block_index}", accession, origin_cds, origin_block_ids, chromosome, current_start, current_end))
                    current_start = row["block_genomic_start"]
                    current_end = row["block_genomic_end"]
                    origin_cds = [row["cds"]]
                    origin_block_ids = [row["block_index"]]
                    block_index += 1
    # add the last block
    chimeric_blocks.append((f"chimblock_{block_index}", group_name, origin_cds, origin_block_ids, chromosome, current_start, current_end))
    block_index += 1

blocks_df = pd.DataFrame(chimeric_blocks, columns=["block_name", "accession", "origin_cds", "origin_blocks", "chromosome", "genomic_start", "genomic_end"])
chim_blocks_df = blocks_df[blocks_df["origin_cds"].apply(len) > 1]
chim_blocks_df.reset_index(drop=True, inplace=True)
chim_blocks_records = []

for accession, blocks in chim_blocks_df.groupby("accession"):
    accession = re.sub(r"\.\d+$", "", accession)
    ref_transcriptome = glob.glob(f"{base_gtf_dir}/{accession}.*/{accession}.*_*_genomic.fna")
    if len(ref_transcriptome) != 1:
        print(ref_transcriptome)
        raise ValueError(f"Multiple or no genome files found for {accession}")
    for record in SeqIO.parse(ref_transcriptome[0], "fasta"):
        subset = blocks[blocks["chromosome"] == record.id]
        for index, row in subset.iterrows():
            chim_block_seq = record.seq[int(row["genomic_start"]):int(row["genomic_end"])]
            chim_block_name = f"{record.id}_chimeric_block{index}"
            chim_blocks_records.append(SeqIO.SeqRecord(chim_block_seq, chim_block_name, description=""))

chim_blocks_df.to_csv(f"{prefix}_chimeric_blocks.tsv", sep="\t", index=False)
SeqIO.write(chim_blocks_records, f"{prefix}_chimeric_blocks.fa", "fasta")
