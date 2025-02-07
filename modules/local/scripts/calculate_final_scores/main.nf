process CALCULATEFINALSCORES {
    tag "$meta.id"
    label "process_medium"
    //conda "${moduleDir}/environment.yml"
    //TODO: create a custom container with pandas and jq
    container "quay.io/biocontainers/pandas:2.2.1"

    input:
    tuple val(meta), path(merged_scores), path(blocks_df), path(chim_blocks_df)

    output:
    tuple val(meta), path("${prefix}_orb_scores.tsv")       , emit: final_scores

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    python /container/bin/calculate_ratios.py  ${merged_scores} ${blocks_df} ${chim_blocks_df} > ${prefix}_orb_scores.tsv 
    """
}
