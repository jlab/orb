process EXTRACTMAPPEDIDS {
    tag "$meta.id"
    label "process_low"

    container "quay.io/biocontainers/jq:1.5--4"

    input:
    tuple val(meta), path(mapped_contigs), path(mapped_contigs_chim)

    output:
    tuple val(meta), path("${prefix}_mapped_ids.txt")       , emit: mapped_ids
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    jq -r 'keys[]' ${mapped_contigs} ${mapped_contigs_chim} > ${prefix}_mapped_ids.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jq: \$( jq --version | sed 's/jq-//' )
    END_VERSIONS
    """
}
