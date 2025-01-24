process GENERATEBLOCKS {
    tag "$meta.id"
    label "process_high"

    //TODO: build custom python container for steps with required packages
    container "quay.io/biocontainers/pandas:2.2.1"

    input:
    tuple val(meta), path(merged_beds), path(reference)
    tuple val(meta2), path(dup_seqs)

    output:
    tuple val(meta), path("${prefix}.fa")       , emit: blocks

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    pip install biopython==1.83
    python /container/bin/generate_blocks.py ${merged_beds} ${reference} ${dup_seqs} > ${prefix}.fa
    """
}
