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
    path  "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    pip install matplotlib
    python /container/bin/make_histogramms.py ${prefix} ${dfs}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    python: "\$(python3 --version | sed 's/Python //')"
    matplotlib: "\$(python3 -c 'import matplotlib; print(matplotlib.__version__)')"
    pandas: "\$(python3 -c 'import pandas as pd; print(pd.__version__)')"
    END_VERSIONS
    """
}
