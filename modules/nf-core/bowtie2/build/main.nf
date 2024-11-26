process BOWTIE2_BUILD {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bowtie2:2.5.2--py39h6fed5c7_0' :
        'biocontainers/bowtie2:2.5.2--py39h6fed5c7_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("$meta.id")    , emit: index
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir $meta.id
    bowtie2-build $args --threads $task.cpus $fasta $meta.id/${fasta.baseName}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir $meta.id
    touch $meta.id/${fasta.baseName}.{1..4}.bt2
    touch $meta.id/${fasta.baseName}.rev.{1,2}.bt2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
    END_VERSIONS
    """
}
