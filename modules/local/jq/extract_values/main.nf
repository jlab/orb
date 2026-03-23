process EXTRACTMAPPEDVALUES {
    tag "$meta.id"
    label "process_low"

    container "quay.io/biocontainers/jq:1.5--4"

    input:
    tuple val(meta), path(json)

    output:
    tuple val(meta), path("${prefix}_values.txt")           , emit: values
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    jq -r 'keys[]' ${json} > ${prefix}_values.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jq: \$( jq --version | sed 's/jq-//' )
    END_VERSIONS
    """
}
