import pandas as pd
import sys
import re
from Bio import SeqIO

all_contigs_ids = sys.argv[1]
mapped_ids = sys.argv[2]
chimeric_mapped_ids = sys.argv[3]
length_filtered_ids = sys.argv[4]
assembler_mapping = sys.argv[5]
dup_ids_group = sys.argv[6]
gene_summary = sys.argv[7]
contig_fasta = sys.argv[8]
prefix = sys.argv[9]

all_contigs_ids = pd.read_csv(all_contigs_ids, header=None).iloc[:, 0]
len_all_contigs_ids = len(all_contigs_ids)

try:
    mapped_ids = pd.read_csv(mapped_ids, header=None).iloc[:, 0]
except pd.errors.EmptyDataError:
    mapped_ids = pd.Series([])
try:
    chimeric_mapped_ids = pd.read_csv(chimeric_mapped_ids, header=None).iloc[:, 0]
except pd.errors.EmptyDataError:
    chimeric_mapped_ids = pd.Series([])
try:
    length_filtered_ids = pd.read_csv(length_filtered_ids, header=None).iloc[:, 0]
except pd.errors.EmptyDataError:
    length_filtered_ids = pd.Series([])

gene_summary = pd.read_csv(gene_summary, usecols=["gene_name", "orthogroup"]).set_index("gene_name")
gene_og_dict = gene_summary.to_dict()["orthogroup"]

# remove ids that are already mapped
chimeric_mapped_ids = chimeric_mapped_ids[~chimeric_mapped_ids.isin(mapped_ids)]

chimeric_mapped_ids_double_mapped = chimeric_mapped_ids[chimeric_mapped_ids.isin(mapped_ids)]

if len(chimeric_mapped_ids_double_mapped) > 0:
    print("Some contigs are double mapped, are removed from summary")
    print(chimeric_mapped_ids_double_mapped)

assembler_mapping = pd.read_csv(assembler_mapping, sep="\t", header=None)


duplicated_ids_dict = {}
for line in open(dup_ids_group):
    line = line.strip()
    ids = line.split("\t")[1]
    ids_list = ids.split(", ")
    dedupped_id = ids_list[0]
    for id in ids_list[1:]:
        duplicated_ids_dict[id] = dedupped_id


all_contigs_ids = all_contigs_ids[~all_contigs_ids.isin(pd.concat([mapped_ids, chimeric_mapped_ids, length_filtered_ids]))]

unmapped_contigs = all_contigs_ids[~all_contigs_ids.isin(assembler_mapping[0].unique())]

assembler_mapping = assembler_mapping[assembler_mapping[0].isin(all_contigs_ids)]

assembler_mapping["origin_gene"] = assembler_mapping[3].apply(lambda x:   re.sub(r"(.*?_.*?)_.*", r"\1", x))
assembler_mapping["origin_gene"] = assembler_mapping["origin_gene"].apply(lambda x: duplicated_ids_dict.get(x, x))

assembler_mapping["origin_orthogroup"] = assembler_mapping["origin_gene"].apply(lambda x: gene_og_dict.get(x))

contig_aggregation = assembler_mapping.groupby(0).agg(
    origin_gene_nunique=("origin_gene", "nunique"),
    origin_orthogroup_nunique=("origin_orthogroup", "nunique")
)

single_mapped_contigs = contig_aggregation[contig_aggregation["origin_gene_nunique"] == 1].index
multi_mapped_contigs = contig_aggregation[contig_aggregation["origin_gene_nunique"] != 1].index

number_ogs_per_multimapped = contig_aggregation[contig_aggregation["origin_gene_nunique"] != 1]["origin_orthogroup_nunique"].value_counts()

multi_og_contigs = number_ogs_per_multimapped[number_ogs_per_multimapped.index != 1].sum()
single_og_contigs = number_ogs_per_multimapped[number_ogs_per_multimapped.index == 1].sum()

mapped_contig_lengths = []
unmpapped_contig_lengths_single = []
unmpapped_contig_lengths_multi = []

single_mapped_contigs = set(single_mapped_contigs.to_series().astype(str))

mapped_ids = set(mapped_ids.astype(str))

for record in SeqIO.parse(contig_fasta, "fasta"):
    record_length = len(record.seq)
    if record.id in mapped_ids:
        mapped_contig_lengths.append(record_length)
    else:
        if record.id in single_mapped_contigs:
            unmpapped_contig_lengths_single.append(record_length)
        else:
            unmpapped_contig_lengths_multi.append(record_length)


result_dict = {"total_contigs": len_all_contigs_ids, "total_bases": sum(mapped_contig_lengths) + sum(unmpapped_contig_lengths_single) + sum(unmpapped_contig_lengths_multi),
               "mapped_contigs": len(mapped_ids), "chimeric_mapped_contigs": len(chimeric_mapped_ids), "length_filtered_contigs": len(length_filtered_ids),
               "unmapped_contigs": len(unmapped_contigs), "multi_mapped_contigs": len(multi_mapped_contigs), "multi_mapped_contigs_single_og":
               single_og_contigs, "multi_mapped_contigs_multi_og": multi_og_contigs, "single_mapped_contigs": len(single_mapped_contigs),
               "mapped_contig_bases": sum(mapped_contig_lengths), "unmapped_contig_bases": sum(unmpapped_contig_lengths_single) + sum(unmpapped_contig_lengths_multi),}

pd.DataFrame.from_dict(result_dict, orient="index", columns=[prefix]).to_csv(f"{prefix}_scores.tsv", sep="\t")

pd.DataFrame({prefix: mapped_contig_lengths + unmpapped_contig_lengths_single + unmpapped_contig_lengths_multi}).to_csv(f"{prefix}_all_contig_lengths.tsv", sep="\t", index=False)
pd.DataFrame({prefix: mapped_contig_lengths}).to_csv(f"{prefix}_mapped_contig_lengths.tsv", sep="\t", index=False)
pd.DataFrame({prefix: unmpapped_contig_lengths_single}).to_csv(f"{prefix}_unmapped_contig_lengths_single.tsv", sep="\t", index=False)
pd.DataFrame({prefix: unmpapped_contig_lengths_multi}).to_csv(f"{prefix}_unmapped_contig_lengths_multi.tsv", sep="\t", index=False)

pd.DataFrame({prefix: pd.concat([pd.Series(mapped_contig_lengths + unmpapped_contig_lengths_single + unmpapped_contig_lengths_multi).describe(), pd.Series({"total_length": sum(mapped_contig_lengths + unmpapped_contig_lengths_single + unmpapped_contig_lengths_multi)})])}).to_csv(f"{prefix}_all_contig_lengths_summary.tsv", sep="\t")
pd.DataFrame({prefix: pd.concat([pd.Series(mapped_contig_lengths).describe(), pd.Series({"total_length": sum(mapped_contig_lengths)})])}).to_csv(f"{prefix}_mapped_contig_lengths_summary.tsv", sep="\t")
pd.DataFrame({prefix: pd.concat([pd.Series(unmpapped_contig_lengths_multi).describe(), pd.Series({"total_length": sum(unmpapped_contig_lengths_multi)})])}).to_csv(f"{prefix}_unmapped_contig_lengths_multi_summary.tsv", sep="\t")
pd.DataFrame({prefix: pd.concat([pd.Series(unmpapped_contig_lengths_single).describe(), pd.Series({"total_length": sum(unmpapped_contig_lengths_single)})])}).to_csv(f"{prefix}_unmapped_contig_lengths_single_summary.tsv", sep="\t")