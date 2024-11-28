process ASSEMBLYREADSTATS {
    tag "$meta.id"
    label "process_low"

    //conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/pandas:2.2.1"

    input:
    tuple val(meta), path(mapping)
    tuple val(meta2), path(summary)

    output:
    tuple val(meta), path("${prefix}.tsv")   , emit: counts
    tuple val(meta), path("${prefix}_merged.tsv")   , emit: merged_mapping

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    cat ${mapping} > ${prefix}_merged.tsv
    python /container/bin/count_genes_and_ogs_per_contig.py ${summary} ${prefix}_merged.tsv > ${prefix}.tsv
    """
}
