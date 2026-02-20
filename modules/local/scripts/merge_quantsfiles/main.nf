process MERGEQUANTSFFILES {
    tag "$meta.id"
    label "process_low"
 //   cache false

    //conda "${moduleDir}/environment.yml"
    container "quay.io/tensulin/orb_toolchain:1.0"

    input:
    tuple val(meta), path(quant_files)

    output:
    tuple val(meta), path("${prefix}.tsv")       , emit: merged_quant_files
    path  "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    merge_quant_files.py ${quant_files} > ${prefix}.tsv


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    python: "\$(python3 --version | sed 's/Python //')"
    pandas: "\$(python3 -c 'import pandas as pd; print(pd.__version__)')"
    END_VERSIONS
    
    """
}
