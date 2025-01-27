process EXTRACTMAPPEDIDS {
    tag "$meta.id"
    label "process_low"

    input:
    tuple val(meta), path(mapped_contigs), path(mapped_contigs_chim)

    output:
    tuple val(meta), path("${prefix}_mapped_ids.txt")       , emit: mapped_ids

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    jq -r 'keys[]' ${mapped_contigs} ${mapped_contigs_chim} > ${prefix}_mapped_ids.txt
    """
}
