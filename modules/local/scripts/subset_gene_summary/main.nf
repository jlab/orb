process SUBSETGENESUMMARY {
    tag "$meta.id"
    label "process_low"

    //conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/pandas:2.2.1"
    //container "ubuntu:plucky-20241111"
    //TODO: change container?   container='ubuntu:24.10'

    input:
    tuple val(meta), path(df)

    output:
    tuple val(meta), path("${prefix}_count.tsv")       , emit: df

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    python /container/bin/subset_gene_summary.py ${df} ${prefix}
    """
}
