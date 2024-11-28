import pandas as pd
import sys
import re

gene_summary = sys.argv[1]
mapping_result = sys.argv[2]

gene_summary = pd.read_csv(gene_summary)
mapping_result = pd.read_csv(mapping_result, sep="\t", header=None)

gene_dict = gene_summary.set_index("gene_name")["orthogroup"].to_dict()

mapping_result["origin_gene"] = mapping_result[3].apply(lambda x:   re.sub(r"(.*?_.*?)_.*", r"\1", x))
mapping_result["orthogroup"] = mapping_result["origin_gene"].apply(lambda x: gene_dict[x])
mapping_result["contig_origin_gene"] = mapping_result[0].astype(str) + "_" + mapping_result["origin_gene"].astype(str)
filtered = mapping_result.drop_duplicates(subset="contig_origin_gene")


#TODO: maybe also add the actual names of the origin cds and the orthogroups
result = mapping_result.groupby(0).agg(
    origin_gene_count=("origin_gene", "count"),
    origin_gene_nunique=("origin_gene", "nunique"),
    orthogroup_nunique=("orthogroup", "nunique")
)
result["origin_genes"] = result.index.map(lambda x: ",".join(filtered[filtered[0] == x]["origin_gene"].to_list()))
result["origin_og"] = result.index.map(lambda x: ",".join(set(filtered[filtered[0] == x]["orthogroup"].to_list())))

result.to_csv(sys.stdout, sep="\t")
