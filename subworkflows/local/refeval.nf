include { BOWTIE2_BUILD                } from '../../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN                } from '../../modules/nf-core/bowtie2/align/main'
include { SEQKIT_GREP as SEQKIT_GREP2  } from '../../modules/nf-core/seqkit/grep/main'
include { GENERATEFAKEGTF              } from '../../modules/local/scripts/generatefakegtf/main'
include { SUBREAD_FEATURECOUNTS        } from '../../modules/nf-core/subread/featurecounts/main'
include { EDGER                        } from '../../modules/local/dge/edger/main'
include { DESEQ2                       } from '../../modules/local/dge/deseq2/main'
include { SALMON_INDEX                 } from '../../modules/local/salmon/index/main'
include { SALMON_QUANT                 } from '../../modules/local/salmon/quant/main'
include { SALMON_QUANT as SALMON_QUANT_ALIGN } from '../../modules/nf-core/salmon/quant/main'
include { MERGEQUANTSFFILES            } from '../../modules/local/scripts/merge_quantsfiles/main'

workflow REFEVAL {
    take:
    reference
    reads
    bowtie2_alignment

    main:
    def prefix = params.run_prefix ?: "prefix"

    GENERATEFAKEGTF(
        reference
    )

    bowtie2_alignment.map {
        [["id": prefix  + "_reference"], it[1]]
    }.groupTuple()
    .set { bowtie2_alignment_with_ref }
    
    bowtie2_alignment_with_ref.join (
        GENERATEFAKEGTF.out.gtf, by: 0
    ).set { ref_with_fake_gtf }

    SUBREAD_FEATURECOUNTS(
        ref_with_fake_gtf
    )
    
    bowtie2_alignment.combine(
        Channel.fromPath("/homes/tlin/Projects/refbasedassemblereval/empty_placeholder.txt")
    ).combine(
        GENERATEFAKEGTF.out.gtf.map { it[1] }
    ).combine(
        reference.map { it[1] }
    ).combine(
        Channel.of(true)
    ).combine(
        Channel.of("A")
    ).multiMap { it ->
        bam: tuple(it[0], it[1])
        index: it[2]
        gtf: it[3]
        transcript_fasta: it[4]
        alignment_mode: it[5]
        lib_type: it[6] 
    }.set { bam_with_gtf }

    SALMON_QUANT_ALIGN(
        bam_with_gtf.bam,
        bam_with_gtf.index,
        bam_with_gtf.gtf,
        bam_with_gtf.transcript_fasta,
        bam_with_gtf.alignment_mode,
        bam_with_gtf.lib_type
    )
    
    SALMON_INDEX (
        reference
    )

    reads.combine(
        SALMON_INDEX.out.index.map {
            [it[1]]
        }
    ).set{ reads_index_combined }

    SALMON_QUANT(
        reads_index_combined
    )

    SALMON_QUANT_ALIGN.out.results.map {
        [["id": prefix + "_quant_alignment"], it[1]]
    }.groupTuple().concat(
        SALMON_QUANT.out.results.map {
            [["id": prefix + "_quant"], it[1]]
        }.groupTuple()
    ).set { quant_files }
    
    MERGEQUANTSFFILES(
        quant_files
    )
}