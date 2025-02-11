include { BOWTIE2_BUILD                } from '../../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN                } from '../../modules/nf-core/bowtie2/align/main'
include { SEQKIT_GREP as SEQKIT_GREP2  } from '../../modules/nf-core/seqkit/grep/main'
include { GENERATEFAKEGTF              } from '../../modules/local/scripts/generatefakegtf/main'
include { SUBREAD_FEATURECOUNTS        } from '../../modules/nf-core/subread/featurecounts/main'
include { EDGER                        } from '../../modules/local/dge/edger/main'
include { DESEQ2                       } from '../../modules/local/dge/deseq2/main'

workflow DGEEVAL {
    take:
    contigs
    mapped_ids
    reads

    main:
    contigs.join(
        mapped_ids
    ).multiMap { it ->
        contigs: tuple(it[0], it[1])
        mapped_ids: tuple(it[2])
    }.set { contigs_mapped_ids }

    SEQKIT_GREP2(
        contigs_mapped_ids.contigs,
        contigs_mapped_ids.mapped_ids
    )

    BOWTIE2_BUILD(
        SEQKIT_GREP2.out.filter
    )

    reads
    .combine(
        BOWTIE2_BUILD.out.index
    )
    .combine(Channel.of(true))
    .combine(Channel.of(false))
    .multiMap { it ->
        reads: tuple(
             [
                "id"           : it[0]["id"] + "_" + it[2]["id"],
                "single_end"   : it[0]["single_end"],
                "read_id"      : it[0]["id"],
                "assembler_id" : it[2]["id"]
            ],
            it[1]
        )
        index: tuple(it[2], it[3])
        save_unaligned: it[4]
        sort_bam: it[5]
    }.set { reads_index_ch }

    BOWTIE2_ALIGN(
        reads_index_ch.reads,
        reads_index_ch.index,
        reads_index_ch.save_unaligned,
        reads_index_ch.sort_bam
    )

    GENERATEFAKEGTF(
        SEQKIT_GREP2.out.filter
    )
    
    BOWTIE2_ALIGN.out.aligned.map { 
        tuple(["id": it[0]["assembler_id"]], it[1])
    }
    .groupTuple()
    .combine(
        GENERATEFAKEGTF.out.gtf, by: 0
    )
    .set { aligned_gtf_ch }

    SUBREAD_FEATURECOUNTS(
        aligned_gtf_ch
    )

    EDGER(
        SUBREAD_FEATURECOUNTS.out.counts
    )

    DESEQ2(
        SUBREAD_FEATURECOUNTS.out.counts
    )
    
    emit:
    counts_ch = SUBREAD_FEATURECOUNTS.out.counts
}