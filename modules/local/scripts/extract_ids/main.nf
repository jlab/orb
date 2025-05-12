process EXTRACTIDS {
    tag "$meta.id"
    label "process_low"

    container "quay.io/biocontainers/jq:1.5--4"

    input:
    tuple val(meta), path(contigs), path(contigs_pre_length_filtering), path(contigs_post_length_filtering)

    output:
    tuple val(meta), path("${prefix}_contigs_ids.txt"), path("${prefix}_length_filtered_ids.txt") , emit: contig_ids
    path  "versions.yml"                                                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    grep "^>" ${contigs} | sed 's/>//' | awk '{print \$1}' > ${prefix}_contigs_ids.txt
    grep "^>" ${contigs_pre_length_filtering} > ids.txt
    grep "^>" ${contigs_post_length_filtering} >> ids.txt
    sort ids.txt | uniq -u | sed s/">"//g | awk '{print \$1}' > ${prefix}_length_filtered_ids.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jq: \$( jq --version | sed 's/jq-//' )
    END_VERSIONS
    """
}
