process MAKEHISTOGRAMS {
    tag "$meta.id"
    label "process_low"

    //conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/pandas:2.2.1"

    input:
    tuple val(meta), path(dfs)

    output:
    tuple val(meta), path("${prefix}_histograms.png")        , emit: hists
    tuple val(meta), path("${prefix}_log_histograms.png")    , emit: hists_log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    pip install matplotlib
    python /container/bin/make_histogramms.py ${prefix} ${dfs}
    """
}
