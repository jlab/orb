process SALMON_INDEX {
    tag "$meta.id:$transcript_fasta"
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/salmon:1.10.1--h7e5ed60_0' :
        'biocontainers/salmon:1.10.1--h7e5ed60_0' }"

    input:
    tuple val(meta), path(transcript_fasta)

    output:
    tuple val(meta), path("salmon")   , emit: index
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    salmon \\
        index \\
        --threads $task.cpus \\
        -t $transcript_fasta \\
        $args \\
        -i salmon

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
    END_VERSIONS
    """
}
