process CALCULATEDGECONFUSION {
    tag "$meta.id"
    label "process_medium"
    container "quay.io/tensulin/orb_toolchain:1.0"

    input:
    tuple val(meta), path(dge_ass), path(mapping)
    tuple val(meta2), path(dge_ref)

    output:
    tuple val(meta), path("*_linear_cm.tsv")                          , emit: linear_cm

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    calculate_dge_confusion.py ${dge_ass} ${mapping} ${dge_ref} ${prefix} ${args}    
    """
}
