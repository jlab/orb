process GATHERRESULTS {
    tag "$meta.id"
    label "process_medium"
    //conda "${moduleDir}/environment.yml"
    //TODO: create a custom container with pandas and jq
    container "quay.io/biocontainers/pandas:2.2.1"

    input:
    tuple val(meta), path(contig_ids), path(length_filtered_contig_ids), path(mapped_scores), path(mapped_chim_scores), path(assembler_mapping), path(dup_seqs), path(gene_summary)

    output:
    tuple val(meta), path("${prefix}_scores.tsv")       , emit: scores

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    python /container/bin/gather_results.py ${contig_ids} ${mapped_scores} ${mapped_chim_scores} ${length_filtered_contig_ids} \\
                                            ${assembler_mapping} ${dup_seqs} ${gene_summary} ${prefix} 
    """
}
