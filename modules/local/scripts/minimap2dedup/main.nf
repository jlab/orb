process MINIMAP2FILTER {
    tag "$meta.id"
    label "process_low"

    //conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/pandas:0.23.4--py36hf8a1672_0"

    input:
    tuple val(meta), path(mapping)

    output:
    tuple val(meta), path("*.json")   , emit: map
    path  "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    python /container/bin/minimap2_filter.py ${mapping}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    python: "\$(python3 --version | sed 's/Python //')"
    pandas: "\$(python3 -c 'import pandas as pd; print(pd.__version__)')"
    END_VERSIONS
    """
}
