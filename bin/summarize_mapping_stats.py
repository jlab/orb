import pandas as pd
import sys
import os

summary_df = sys.argv[1]

assembler_name = sys.argv[2]

gene_counts = gene_summary[[col for col in list(gene_summary.columns) if col.startswith("sample")]]

origin_gene_stats = summary_df["origin_gene_nunique"].describe()
origin_og_stats = summary_df["orthogroup_nunique"].describe()

exploded_genes = summary_df["origin_genes"].str.split(",").explode()

number_of_genes_with_one_read_mapped = exploded_genes.nunique()

number_of_reads_with_duplicated_mapping = exploded_genes[exploded_genes.duplicated()].nunique()

summary_genes = pd.DataFrame([summary_df["origin_gene_nunique"].describe(), summary_df["orthogroup_nunique"].describe()]).T

summary_genes.columns = [f"{assembler_name}_{col}" for col in summary_genes.columns]

summary_mapping_stats = pd.DataFrame([number_of_genes_with_one_read_mapped, number_of_reads_with_duplicated_mapping], columns=[f"{assembler_name}_mapping_stats"],
                                     index=["number_of_genes_with_one_read_mapped", "number_of_reads_with_duplicated_read_mapping"])

summary_mapping_stats.to_csv(f"{assembler_name}_mapping_stats.csv")
summary_genes.to_csv(f"{assembler_name}_contig_gene_stats.csv")