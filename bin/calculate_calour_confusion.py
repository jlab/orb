import pandas as pd
import json
import sys
from sklearn.metrics import confusion_matrix

assembler_dge_file = sys.argv[1]
map_file = sys.argv[2]
gene_summary_file = sys.argv[3]

prefix = sys.argv[4]

row_prefix = "calour_"

gene_summary = pd.read_csv(gene_summary_file, sep='\t')
assembler_dge = pd.read_csv(assembler_dge_file, sep='\t')


with open(map_file, "r") as file:
    contig_to_gene = json.load(file)

assembler_dge["Name"] = assembler_dge["Name"].astype(str)
 
assembler_dge["original_cds"] = assembler_dge["Name"].map(contig_to_gene)

assembler_dge["gene_name"] = assembler_dge["original_cds"].str.replace(r"_block\d+$", "", regex=True)

assembler_dge = assembler_dge.drop_duplicates(subset="gene_name")["gene_name"].to_list()

reference_dge = gene_summary["gene_name"].to_list()

not_labels = [dge for dge in assembler_dge if dge not in reference_dge]

true_labels = [True for _ in reference_dge] + [False for _ in not_labels]
pred_labels = [dge in assembler_dge for dge in reference_dge] + [dge in assembler_dge for dge in not_labels]

cm = confusion_matrix(true_labels, pred_labels)

linear_cm = {f"{row_prefix}DE_TN": cm[0, 0], f"{row_prefix}DE_FP": cm[0, 1], f"{row_prefix}DE_FN": cm[1, 0], f"{row_prefix}DE_TP": cm[1, 1]}

pd.DataFrame.from_dict(linear_cm, orient="index", columns=[prefix]).to_csv(f"{row_prefix}{prefix}_linear_cm.tsv", sep="\t")
