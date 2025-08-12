process IDENTIFYCHIMERICBLOCKS {
    tag "$meta.id"
    label "process_high"

    //TODO: build custom python container for steps with required packages
    container "quay.io/biocontainers/pandas:2.2.1"

    input:
    tuple val(meta), path(blocks_tsv)
    tuple val(meta2), path(dup_seqs)
    path base_gtf_dir
    path translation_df
    path gene_summary
    val min_chimeric_overlap

    output:
    tuple val(meta), path("${prefix}_chimeric_blocks.fa")    , emit: chim_blocks
    tuple val(meta), path("${prefix}_chimeric_blocks.tsv")   , emit: chim_blocks_tsv
    path "versions.yml"                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    pip install biopython==1.83
    python /container/bin/identify_chimeric_blocks.py ${base_gtf_dir} ${gene_summary} ${translation_df} ${blocks_tsv} ${dup_seqs} ${min_chimeric_overlap} ${prefix}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    python: "\$(python3 --version | sed 's/Python //')"
    biopython: "\$(python3 -c 'import Bio; print(Bio.__version__)')"
    pandas: "\$(python3 -c 'import pandas as pd; print(pd.__version__)')"
    END_VERSIONS
    """
}
