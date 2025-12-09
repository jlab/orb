process MERGE_ASSEMBLER_OG_COUNTS {
    tag "$meta.id"
    label "process_medium"
    //conda "${moduleDir}/environment.yml"
    //TODO: create a custom container with pandas and jq
    container "quay.io/tensulin/orb_toolchain:1.0"

    input:
    tuple val(meta), path(counts), path(assembler_mapping), path(gene_summary)

    output:
    tuple val(meta), path("${prefix}_merged_orthogroups.tsv")              , emit: merged_assembler_orthogroups
    path  "versions.yml"                                                   , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    merge_assembler_og_counts.py ${counts} ${assembler_mapping} ${gene_summary} ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    python: "\$(python3 --version | sed 's/Python //')"
    polars: "\$(python3 -c 'import polars as pl; print(pl.__version__)')"
    END_VERSIONS
    """
}
