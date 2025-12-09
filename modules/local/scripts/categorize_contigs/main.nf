process CATEGORIZECONTIGS {
    tag "$meta.id"
    label "process_medium"
    //conda "${moduleDir}/environment.yml"
    container "quay.io/tensulin/orb_toolchain:1.0"

    input:
    tuple val(meta), path(contig_ids), path(length_filtered_contig_ids), path(mapped_scores), path(mapped_chim_scores), path(assembler_mapping), path(gene_summary), path(contigs_fasta)

    output:
    tuple val(meta), path("${prefix}_contigs_categorised.tsv")           , emit: contig_categorisation
    path  "versions.yml"                                                   , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    categorize_contigs.py ${contig_ids} ${mapped_scores} ${mapped_chim_scores} ${length_filtered_contig_ids} \\
                                                ${assembler_mapping} ${gene_summary} ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    python: "\$(python3 --version | sed 's/Python //')"
    polars: "\$(python3 -c 'import polars as pl; print(pl.__version__)')"
    END_VERSIONS
    """
}
