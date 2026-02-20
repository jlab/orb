process MERGEMAPPINGLOGS {
    tag "$meta.id"
    label "process_low"

    //conda "${moduleDir}/environment.yml"
    container "quay.io/tensulin/orb_toolchain:1.0"

    input:
    tuple val(meta), path(logs)

    output:
    tuple val(meta), path("${prefix}_mean_logs.tsv")       , emit: log_mean
    tuple val(meta), path("${prefix}_merged_logs.tsv")     , emit: log_merged

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    merge_bowtie2_logs.py ${logs} ${prefix}
    """
}
