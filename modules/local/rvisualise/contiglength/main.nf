process LENGTHVISUALIZE_ASSEMBLY_CONTIG {
    input:
    tuple val(meta), path(contigs)

    output:
    path "*${meta.id}_contig_lengths_histogram.png"

    script:
    """
    contig_length_visualiser.R ${contigs} ${meta.id}_contig_lengths_histogram.png
    """
}
