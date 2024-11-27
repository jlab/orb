process SCRIPT_MAPPING_STATS {
    tag "$meta.id"
    label "process_low"

    //conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/pandas:2.2.1"

    input:
    tuple val(meta), path(mapping)
    tuple val(meta2), path(summary)
    tuple val(meta3), path(logs)

    output:
    tuple val(meta), path("${prefix}.tsv")   , emit: counts
    path "${prefix}_log.tsv"                 , emit: logs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    python /container/bin/count_genes_and_ogs_per_contig.py ${summary} ${mapping} > ${prefix}.tsv
    python /container/bin/parse_bowtie2_logs.py ${logs} ${prefix} > ${prefix}_log.tsv
    """
}
