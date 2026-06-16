process MINIMAP2CLASSIFICATION {
    tag "$meta.id"
    label "process_high"
    //conda "${moduleDir}/environment.yml"
    container "quay.io/tensulin/orb_toolchain:1.0"

    input:
    tuple val(meta), path(mapping), path(overlap_selection), path(contigs), path(gene_summary)

    output:
    tuple val(meta), path("${prefix}_recovered.json")   , emit: map
    tuple val(meta), path("${prefix}_minimap2_categories.tsv"), emit: categories
    tuple val(meta), path("${prefix}_minimap2_category_counts.tsv"), emit: category_counts
    tuple val(meta), path("${prefix}_contig_n_stats.tsv"), emit: contig_n_stats
    path  "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // seed used for random selection for cooptimal mappings
    def args = task.ext.args ?: '1234'
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    minimap2_classification.py ${mapping} ${overlap_selection} ${contigs} ${gene_summary} ${args} ${prefix} 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    python: "\$(python3 --version | sed 's/Python //')"
    pandas: "\$(python3 -c 'import pandas as pd; print(pd.__version__)')"
    END_VERSIONS
    """
}
