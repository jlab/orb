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
    tuple val(meta), path("${prefix}_blocks.tsv")   , emit: blocks_tsv
    path  "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    pip install biopython==1.83
    python /container/bin/generate_blocks.py ${merged_beds} ${reference} ${dup_seqs} ${prefix} > ${prefix}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    python: "\$(python3 --version | sed 's/Python //')"
    biopython: "\$(python3 -c 'import Bio; print(Bio.__version__)')"
    pandas: "\$(python3 -c 'import pandas as pd; print(pd.__version__)')"
    END_VERSIONS
    """
}
