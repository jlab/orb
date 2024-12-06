process SUMMARIZEMAPPINGSTATS {
    tag "$meta.id"
    label "process_medium"

    //conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/pandas:2.2.1"

    input:
    tuple val(meta), path(mapping_stats)

    output:
    tuple val(meta), path("${prefix}_mapping_stats.csv")       , emit: mapping_stats
    tuple val(meta), path("${prefix}_contig_gene_stats")     , emit: contig_gene_stats

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    python /container/bin/summarize_mapping_stats.py ${mapping_stats} ${prefix}
    """
}
