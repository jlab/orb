import pandas as pd
import json
import sys
from sklearn.metrics import confusion_matrix

assembler_dge_file = sys.argv[1]
map_file = sys.argv[2]
gene_summary_file = sys.argv[3]

prefix = sys.argv[4]

pval_name = sys.argv[5]
log_name = sys.argv[6]
row_prefix = sys.argv[7]

gene_summary = pd.read_csv(gene_summary_file, sep='\t')
assembler_dge = pd.read_csv(assembler_dge_file, sep='\t')


with open(map_file, "r") as file:
    contig_to_gene = json.load(file)

gene_summary["cds_de"] = (gene_summary[pval_name] < 0.05) & (gene_summary[log_name].abs() > 1)

assembler_dge["ContigName"] = assembler_dge["ContigName"].astype(str)

assembler_dge["original_cds"] = assembler_dge["ContigName"].map(contig_to_gene)

assembler_dge["gene_name"] = assembler_dge["original_cds"].str.replace(r"_block\d+$", "", regex=True)

assembler_dge["de"] = (assembler_dge[pval_name] < 0.05) & (assembler_dge[log_name].abs() > 1)

for cds, group in assembler_dge[assembler_dge["gene_name"].duplicated(keep=False)].groupby("gene_name"):
    if sum(group["de"]) > 0:
        assembler_dge.loc[group.index, "de"] = True

cds_contig_dict = assembler_dge.drop_duplicates(subset="gene_name").set_index("gene_name")[["ContigName", "de"]].apply(list, axis=1).to_dict()

contig_map, de_map = zip(*[cds_contig_dict.get(gene, ["", False]) for gene in gene_summary["ContigName"]])

gene_summary["assembler_contig"], gene_summary["assembler_de"] = contig_map, de_map

true_labels = (gene_summary["cds_de"]).astype(int)
pred_labels = (gene_summary["assembler_de"]).astype(int)

cm = confusion_matrix(true_labels, pred_labels, labels=[0, 1])

linear_cm = {f"{row_prefix}DE_TN": cm[0, 0], f"{row_prefix}DE_FP": cm[0, 1], f"{row_prefix}DE_FN": cm[1, 0], f"{row_prefix}DE_TP": cm[1, 1]}

pd.DataFrame.from_dict(linear_cm, orient="index", columns=[prefix]).to_csv(f"{row_prefix}{prefix}_linear_cm.tsv", sep="\t")
