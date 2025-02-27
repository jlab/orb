process CALOUR {
    tag "$meta.id"
    label "process_medium"
        
    environment = [ 'PATH' : '/bin:/usr/bin:/usr/local/bin' ]
    container "quay.io/tensulin/calour:latest"

    input:
    tuple val(meta), path(count_matrix), val(min_abundance)


    output:
    tuple val(meta), path("${prefix}_calour_full_table.tsv")       , emit: results
    path "versions.yml"                                            , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    python3 -c "import calour; print(f'Calour: \"{calour.__version__}\"')" > versions.yml
    python3 /container/bin/calour_calculation.py ${count_matrix} ${min_abundance} ${prefix} ${args}
    """
}
