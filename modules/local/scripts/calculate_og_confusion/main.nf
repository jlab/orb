process CALCULATEOGCONFUSION {
    tag "$meta.id"
    label "process_medium"
    container "quay.io/biocontainers/pandas:2.2.1"

    input:
    tuple val(meta), path(dge_ass)
    tuple val(meta2), path(dge_ref)

    output:
    tuple val(meta), path("*_linear_cm.tsv")                          , emit: linear_cm

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    pip install scikit-learn
    calculate_og_confusion.py ${dge_ass} ${dge_ref} ${prefix} ${args}
    """
}
