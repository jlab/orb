process STACKDATAFRAMES {
    tag "$meta.id"
    label "process_low"

    //conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/pandas:2.2.1"

    input:
    tuple val(meta), path(dfs)

    output:
    tuple val(meta), path("${prefix}.tsv")       , emit: stacked_dfs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    python /container/bin/stack_dataframes.py ${dfs} > ${prefix}.tsv
    """
}
