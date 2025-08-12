process STACKDATAFRAMESWITHHEADER {
    tag "$meta.id"
    label "process_low"

    container "quay.io/biocontainers/pandas:2.2.1"

    input:
    tuple val(meta), path(dfs)

    output:
    tuple val(meta), path("${prefix}.tsv"), emit: stacked_dfs

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    def first_file = dfs[0]
    def remaining_files = dfs.size() > 1 ? dfs[1..-1].join(" ") : ""

    """
    cp ${first_file} ${prefix}.tsv

    for file in ${remaining_files}; do
        tail -n +2 "\$file" >> ${prefix}.tsv
    done
    """
}
