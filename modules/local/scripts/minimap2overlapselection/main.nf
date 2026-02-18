process MINIMAP2OVERLAPSELECTION {
    tag "$meta.id"
    label "process_high"
    //conda "${moduleDir}/environment.yml"
    container "quay.io/tensulin/orb_toolchain:1.0"

    input:
    tuple val(meta), path(mapping), path(contigs)

    output:
    tuple val(meta), path("${prefix}_overlap_recovered.json")   , emit: map
    tuple val(meta), path("${prefix}_minimap2_overlap_blocks.tsv"), emit: categories
    tuple val(meta), path("${prefix}_overlap_n_counts.tsv"), emit: n_counts
    path  "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // seed used for random selection for cooptimal mappings
    def args = task.ext.args ?: '1234'
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    minimap2_overlap_selection.py ${mapping} ${contigs} ${args} ${prefix} 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    python: "\$(python3 --version | sed 's/Python //')"
    pandas: "\$(python3 -c 'import pandas as pd; print(pd.__version__)')"
    END_VERSIONS
    """
}
