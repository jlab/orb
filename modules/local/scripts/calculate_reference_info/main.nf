process CALCULATEREFERENCEINFO {
    tag "$meta.id"
    label "process_medium"
    //conda "${moduleDir}/environment.yml"
    //TODO: create a custom container with pandas and jq
    container "quay.io/biocontainers/pandas:2.2.1"

    input:
    tuple val(meta), path(ref_fasta)
    tuple val(meta2), path(blocks_df)
    tuple val(meta3) path(chim_blocks_df)

    output:
    tuple val(meta), path("${prefix}_reference_stats.tsv")                             , emit: ref_stats
    tuple val(["id": "CDS"]), path("${prefix}_cds_lengths.tsv")                        , emit: cds_lenghts
    tuple val(["id": "Blocks"]), path("${prefix}_block_lengths.tsv")                   , emit: block_lengths
    tuple val(["id": "Chimeric Blocks"]), path("${prefix}_chimeric_block_lengths.tsv") , emit: chimeric_block_lengths

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    python /container/bin/calculate_reference_info.py  ${blocks_df} ${chim_blocks_df} ${ref_fasta} ${prefix}
    """
}
