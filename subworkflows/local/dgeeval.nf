include { BOWTIE2_BUILD                } from '../../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN                } from '../../modules/nf-core/bowtie2/align/main'
include { SEQKIT_GREP as SEQKIT_GREP2  } from '../../modules/nf-core/seqkit/grep/main'
include { GENERATEFAKEGTF              } from '../../modules/local/scripts/generatefakegtf/main'
include { EDGER; EDGER as REFEDGER     } from '../../modules/local/dge/edger/main'
include { DESEQ2; DESEQ2 as REFDESEQ2  } from '../../modules/local/dge/deseq2/main'
include { SALMON_INDEX                 } from '../../modules/local/salmon/index/main'
include { SALMON_QUANT                 } from '../../modules/local/salmon/quant/main'
include { MERGEQUANTSFFILES            } from '../../modules/local/scripts/merge_quantsfiles/main'
include { CALCULATEDGECONFUSION as CALCULATEDGECONFUSIONDESEQ; CALCULATEDGECONFUSION as CALCULATEDGECONFUSIONEDGER } from '../../modules/local/scripts/calculate_dge_confusion/main'
include { STACKDATAFRAMESWITHHEADER    } from '../../modules/local/scripts/stack_dataframes_with_header/main'
include { MERGEDATAFRAMES              } from '../../modules/local/scripts/merge_dataframes/main'
include { SORTDATAFRAME                } from '../../modules/local/scripts/sort_dataframe/main'
include { CALOUR; CALOUR as REFCALOUR  } from '../../modules/local/dge/calour/main'
include { SUBSETGENESUMMARY            } from '../../modules/local/scripts/subset_gene_summary/main'
include { CALCULATECALOURCONFUSION     } from '../../modules/local/scripts/calculate_calour_confusion/main'

workflow DGEEVAL {
    take:
    contigs
    mapped_ids
    reads
    reference_summary
    contig_cds_map
    ch_versions

    main:
    def prefix = params.run_prefix ?: "prefix"

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

    ch_versions = ch_versions.mix(SEQKIT_GREP2.out.versions)

    SALMON_INDEX(
        SEQKIT_GREP2.out.filter
    )

    reads.combine(
        SALMON_INDEX.out.index
    ).map { it ->
        reads: tuple(
             [
                "id"           : it[2]["id"] + "_" + it[0]["id"],
                "single_end"   : it[0]["single_end"],
                "assembler_id" : it[2]["id"]
            ],
            it[1],
            it[3]
        )
    }
    .set{ reads_index_combined }

    SALMON_QUANT(
        reads_index_combined
    )

    ch_versions = ch_versions.mix(SALMON_QUANT.out.versions)

    SALMON_QUANT.out.results.map {
        tuple(["id": it[0]["assembler_id"]], it[1])
    }.groupTuple().set { quant_files_grouped }
    
    MERGEQUANTSFFILES (
        quant_files_grouped
    )

    ch_versions = ch_versions.mix(MERGEQUANTSFFILES.out.versions)

    EDGER(
        MERGEQUANTSFFILES.out.merged_quant_files
    )

    ch_versions = ch_versions.mix(EDGER.out.versions)

    DESEQ2(
        MERGEQUANTSFFILES.out.merged_quant_files
    )

    ch_versions = ch_versions.mix(DESEQ2.out.versions)

    CALOUR(
        MERGEQUANTSFFILES.out.merged_quant_files.combine(
            Channel.of(params.calour_min_reads)
        )
    )

    ch_versions = ch_versions.mix(CALOUR.out.versions)

    REFEDGER(
        reference_summary
    )

    ch_versions = ch_versions.mix(REFEDGER.out.versions)

    REFDESEQ2(
        reference_summary
    )

    ch_versions = ch_versions.mix(REFDESEQ2.out.versions)

    SUBSETGENESUMMARY(
        reference_summary
    )

    ch_versions = ch_versions.mix(SUBSETGENESUMMARY.out.versions)

    REFCALOUR(
        SUBSETGENESUMMARY.out.df.combine(
            Channel.of(params.calour_min_reads)
        )
    )

    ch_versions = ch_versions.mix(REFCALOUR.out.versions)    

    DESEQ2.out.results.join(contig_cds_map).combine(
        REFDESEQ2.out.results
    ).multiMap { it ->
        ass: tuple(it[0], it[1], it[2])
        ref: tuple(it[3], it[4])
    }.set { deseq2_dge }

    CALOUR.out.results.join(contig_cds_map).combine(
        REFCALOUR.out.results
    ).multiMap { it ->
        ass: tuple(it[0], it[1], it[2])
        ref: tuple(it[3], it[4])
    }.set { calour_dge }

    CALCULATEDGECONFUSIONDESEQ(
        deseq2_dge.ass,
        deseq2_dge.ref
    )

    CALCULATECALOURCONFUSION(
        calour_dge.ass,
        calour_dge.ref
    )

    EDGER.out.results.join(contig_cds_map).combine(
        REFEDGER.out.results
    ).multiMap { it ->
        ass: tuple(it[0], it[1], it[2])
        ref: tuple(it[3], it[4])
    }.set { edger_dge }

    CALCULATEDGECONFUSIONEDGER(
        edger_dge.ass,
        edger_dge.ref
    )

    CALCULATEDGECONFUSIONDESEQ.out.linear_cm.concat(
        CALCULATEDGECONFUSIONEDGER.out.linear_cm
    ).concat(
        CALCULATECALOURCONFUSION.out.linear_cm
    ).groupTuple()
    .set{ cms }

    STACKDATAFRAMESWITHHEADER(
        cms
    )
    
    STACKDATAFRAMESWITHHEADER.out.stacked_dfs.map {
        [["id": prefix], it[1]]
    }.groupTuple()
    .set { stacked_cms } 

    MERGEDATAFRAMES(
        stacked_cms
    )

    SORTDATAFRAME(
        MERGEDATAFRAMES.out.merged_dfs,
        Channel.of(1)
    )

    emit:
    eval = SORTDATAFRAME.out.dataframe
    ch_versions
}