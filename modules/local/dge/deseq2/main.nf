process DESEQ2 {
    tag "$meta.id"
    label "process_medium"

    //conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/bioconductor-deseq2:1.46.0--r44he5774e6_0"

    input:
    tuple val(meta), path(count_matrix)


    output:
    tuple val(meta), path("${prefix}_DESeq2_full_table.tsv")       , emit: results
    path "versions.yml"                                            , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo "" > deletethisline.txt
    echo "Running edgeR calculation"
    Rscript /container/bin/deseq2_calculation.R ${count_matrix} ${prefix} ${args}
    """
}
