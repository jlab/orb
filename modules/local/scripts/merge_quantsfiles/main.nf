process MERGEQUANTSFFILES {
    tag "$meta.id"
    label "process_low"

    //conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/pandas:2.2.1"

    input:
    tuple val(meta), path(quant_files)

    output:
    tuple val(meta), path("${prefix}.tsv")       , emit: merged_quant_files

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    python /container/bin/merge_quant_files.py ${quant_files} > ${prefix}.tsv
    """
}
