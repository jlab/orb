process EDGER {
    tag "$meta.id"
    label "process_medium"
    
    //conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/bioconductor-edger:4.4.0--r44h3df3fcb_0"

    input:
    tuple val(meta), path(count_matrix)


    output:
    tuple val(meta), path("${prefix}_edgeR_full_table.tsv")       , emit: results
    path "versions.yml"                                           , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    Rscript /container/bin/edgeR_calculation.R ${count_matrix} ${prefix} ${args}
    """
}
