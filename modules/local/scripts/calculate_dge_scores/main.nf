process CALCULATEDGESCORES {
    tag "$meta.id"
    label "process_medium"
    container "quay.io/biocontainers/pandas:2.2.1"

    input:
    tuple val(meta), path(contig_ids), path(length_filtered_contig_ids)

    output:
    tuple val(meta), path("${prefix}_scores.tsv")                          , emit: scores
    tuple val(meta), path("${prefix}_all_contig_lengths.tsv")              , emit: all_contig_lengths

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    python /container/bin/calculate_dge_scores.py ${contig_ids} ${mapped_scores} ${mapped_chim_scores} ${length_filtered_contig_ids} 
    """
}
