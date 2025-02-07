process GATHERRESULTS {
    tag "$meta.id"
    label "process_medium"
    //conda "${moduleDir}/environment.yml"
    //TODO: create a custom container with pandas and jq
    container "quay.io/biocontainers/pandas:2.2.1"
    //container "quay.io/biocontainers/biopython:1.68--py35_0"

    input:
    tuple val(meta), path(contig_ids), path(length_filtered_contig_ids), path(mapped_scores), path(mapped_chim_scores), path(assembler_mapping), path(dup_seqs), path(gene_summary), path(contigs_fasta)

    output:
    tuple val(meta), path("${prefix}_scores.tsv")                          , emit: scores
    tuple val(meta), path("${prefix}_all_contig_lengths.tsv")              , emit: all_contig_lengths
    tuple val(meta), path("${prefix}_unmapped_contig_lengths_single.tsv")  , emit: unmapped_contig_lengths_single
    tuple val(meta), path("${prefix}_unmapped_contig_lengths_multi.tsv")   , emit: unmapped_contig_lengths_multi
    tuple val(meta), path("${prefix}_mapped_contig_lengths.tsv")           , emit: mapped_contig_lengths

    tuple val(meta), path("${prefix}_all_contig_lengths_summary.tsv")             , emit: all_contig_lengths_summary
    tuple val(meta), path("${prefix}_unmapped_contig_lengths_single_summary.tsv") , emit: unmapped_contig_lengths_single_summary
    tuple val(meta), path("${prefix}_unmapped_contig_lengths_multi_summary.tsv")  , emit: unmapped_contig_lengths_multi_summary
    tuple val(meta), path("${prefix}_mapped_contig_lengths_summary.tsv")          , emit: mapped_contig_lengths_summary

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    echo "asdas" > test.txt
    pip install Biopython
    python /container/bin/gather_results.py ${contig_ids} ${mapped_scores} ${mapped_chim_scores} ${length_filtered_contig_ids} \\
                                            ${assembler_mapping} ${dup_seqs} ${gene_summary} ${contigs_fasta} ${prefix} 
    """
}
