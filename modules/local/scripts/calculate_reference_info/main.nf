process CALCULATEREFERENCEINFO {
    tag "$meta.id"
    label "process_medium"
    //conda "${moduleDir}/environment.yml"
    //TODO: create a custom container with pandas and jq
    container "quay.io/biocontainers/pandas:2.2.1"

    input:
    tuple val(meta), path(ref_fasta)
    tuple val(meta2), path(blocks_df)
    tuple val(meta3), path(chim_blocks_df)

    output:
    path("${prefix}_reference_stats.tsv")   , emit: ref_stats
    path("${prefix}_cds_lengths.tsv")                        , emit: cds_lenghts
    path("${prefix}_block_lengths.tsv")                      , emit: block_lengths
    path("${prefix}_chimeric_block_lengths.tsv")             , emit: chimeric_block_lengths
    path  "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    pip install biopython==1.83
    python /container/bin/calculate_reference_info.py  ${blocks_df} ${chim_blocks_df} ${ref_fasta} ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    python: "\$(python3 --version | sed 's/Python //')"
    pandas: "\$(python3 -c 'import pandas as pd; print(pd.__version__)')"
    END_VERSIONS
    """
}
