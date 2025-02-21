process SORTDATAFRAME {
    tag "$meta.id"
    label "process_low"

    //conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/pandas:2.2.1"

    input:
    tuple val(meta), path(df)
    val axis

    output:
    tuple val(meta), path("${prefix}_sorted.tsv")       , emit: dataframe

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    python /container/bin/sort_dataframe.py ${df} ${axis} ${prefix}
    """
}
