process MINIMAP2_DEDUP {
    tag "$meta.id"
    label "process_low"

    //conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/pandas:0.23.4--py36hf8a1672_0"

    input:
    tuple val(meta), path(mapping)

    output:
    tuple val(meta), path("*.json")   , emit: map

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    python /container/bin/minimap2_deduplicate.py ${mapping} 
    """
}
