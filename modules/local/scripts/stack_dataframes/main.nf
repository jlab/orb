process STACKDATAFRAMES {
    tag "$meta.id"
    label "process_low"

    //conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/pandas:2.2.1"
    //container "ubuntu:plucky-20241111"
    //TODO: change container?   container='ubuntu:24.10'

    input:
    tuple val(meta), path(dfs)

    output:
    tuple val(meta), path("${prefix}.tsv")       , emit: stacked_dfs
    path  "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat ${dfs} > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    cat: "\$(cat --version | head -n1 | sed 's/cat (GNU coreutils) //')"
    END_VERSIONS
    """
}
