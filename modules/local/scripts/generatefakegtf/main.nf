process GENERATEFAKEGTF {
    tag "$meta.id"
    label "process_low"

    //conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/pandas:2.2.1"
    //container "ubuntu:plucky-20241111"
    //TODO: change container?   container='ubuntu:24.10'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}_fake_annotation.gtf")    ,  emit: gtf

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    awk '                    
    /^>/ {
        if (seq) {
            print header "\\tFAKE\\tgene\\t1\\t" length(seq) "\\t.\\t+\\t.\\tgene_id \\"" header "\\"; transcript_id \\"" header "\\";";
        }
        header = substr(\$1, 2);  
        seq = "";
        next;
    }
    { seq = seq \$0; }
    END {
        if (seq) {
            print header "\\tFAKE\\tgene\\t1\\t" length(seq) "\\t.\\t+\\t.\\tgene_id \\"" header "\\"; transcript_id \\"" header "\\";";
        }
    }
    ' ${fasta} > ${prefix}_fake_annotation.gtf
    """
}
