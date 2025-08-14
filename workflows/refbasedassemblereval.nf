/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowAssemblereval.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK_CONTIGS } from '../subworkflows/local/input_check_contig'
include { INPUT_CHECK_READS   } from '../subworkflows/local/input_check_reads'
include { DGEEVAL             } from '../subworkflows/local/dgeeval'
include { REFEVAL             } from '../subworkflows/local/refeval'
include { LENGTHVISUALIZE_ASSEMBLY_CONTIG as ASSEMLYLENGTHVISUALIZE_ASSEMBLY_CONTIG} from '../modules/local/rvisualise/contiglength'
include { LENGTHVISUALIZE_ASSEMBLY_CONTIG as SAMPLELENGTHVISUALIZE_ASSEMBLY_CONTIG } from '../modules/local/rvisualise/contiglength'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { QUAST as QUAST_SAMPLE       } from '../modules/nf-core/quast/main'
include { QUAST as QUAST_ASSEMBLER    } from '../modules/nf-core/quast/main'
include { QUAST as QUAST              } from '../modules/nf-core/quast/main'
include { SALMON_INDEX                } from '../modules/local/salmon/index/main'
include { SALMON_QUANT                } from '../modules/local/salmon/quant/main'
include { MINIMAP2_MAP; MINIMAP2_MAP as MINIMAP2_MAP_CHIMERIC } from '../modules/local/minimap2/map/main'
include { MINIMAP2FILTER; MINIMAP2FILTER as MINIMAP2_FILTER_CHIMERIC } from '../modules/local/scripts/minimap2dedup/main'
include { BOWTIE2_BUILD; BOWTIE2_BUILD as BOWTIE2_BUILD_REFERENCE } from '../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN; BOWTIE2_ALIGN as BOWTIE2_ALIGN_REFERENCE } from '../modules/nf-core/bowtie2/align/main'
include { BEDTOOLS_BAMTOBED; BEDTOOLS_BAMTOBED as BEDTOOLS_BAMTOBED_REFERENCE } from '../modules/nf-core/bedtools/bamtobed/main'
include { SCRIPT_MAPPING_STATS        } from '../modules/local/scripts/mapping_stats/main'
include { ASSEMBLYREADSTATS           } from '../modules/local/scripts/assembly_read_stats/main'
include { MERGEMAPPINGLOGS            } from '../modules/local/scripts/merge_mapping_logs/main'
include { MERGEALLLOGS                } from '../modules/local/scripts/merge_all_logs/main'
include { SUMMARIZEMAPPINGSTATS       } from '../modules/local/scripts/summarizemapping/main'
include { MERGEDATAFRAMES as MERGEDATAFRAMESMAPPING; MERGEDATAFRAMES as MERGEDATAFRAMECONTIG;
MERGEDATAFRAMES as MERGEDATAFRAMESSUMMARIES } from '../modules/local/scripts/merge_dataframes/main'
include { GENERATEBLOCKS              } from '../modules/local/scripts/generate_blocks/main'
include { STACKDATAFRAMES; STACKDATAFRAMES as STACKDATAFRAMES_REFERENCE } from '../modules/local/scripts/stack_dataframes/main'
include { SEQKIT_RMDUP                } from '../modules/nf-core/seqkit/rmdup/main'
include { IDENTIFYCHIMERICBLOCKS      } from '../modules/local/scripts/identify_chimeric_blocks/main'
include { EXTRACTMAPPEDIDS            } from '../modules/local/scripts/extract_mapped_ids/main'
include { SEQKIT_SEQ                  } from '../modules/nf-core/seqkit/seq/main'
include { SEQKIT_GREP                 } from '../modules/nf-core/seqkit/grep/main'
include { GATHERRESULTS               } from '../modules/local/scripts/gather_results/main'
include { EXTRACTMAPPEDVALUES; EXTRACTMAPPEDVALUES as EXTRACTMAPPEDVALUES_CHIMERIC } from '../modules/local/jq/extract_values/main'
include { EXTRACTIDS                  } from '../modules/local/scripts/extract_ids/main'
include { CALCULATEFINALSCORES        } from '../modules/local/scripts/calculate_final_scores/main'
include { CALCULATEREFERENCEINFO      } from '../modules/local/scripts/calculate_reference_info/main'
include { MAKEHISTOGRAMS              } from '../modules/local/scripts/make_histograms/main'
include { SORTDATAFRAME                } from '../modules/local/scripts/sort_dataframe/main'
include { STACKDATAFRAMESWITHHEADER as MERGEDRESULTS } from '../modules/local/scripts/stack_dataframes_with_header/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow ASSEMBLEREVAL {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //

    INPUT_CHECK_READS (
        file(params.reads)
    )

    ch_versions = ch_versions.mix(INPUT_CHECK_READS.out.versions)

    def prefix = params.run_prefix ?: "prefix"

    sample_assemblers = Channel.fromPath(params.contigs)
    reads = INPUT_CHECK_READS.out.reads
    reference_cds = Channel.fromPath(params.reference_cds)
    gene_summary = Channel.fromPath(params.gene_summary)
    Channel.of(["id": prefix]).combine(gene_summary).set{ gene_summary_prefix }
    
    reference_cds.map {
        assembler ->
        [["id": prefix + "_reference"], assembler]
    }.set{ reference_cds_val }

    sample_assemblers.map {
        assembler ->
        [["id": assembler.baseName.replace('_contigs', '')], assembler]
    }.set{ assembler_contigs }

    

    SEQKIT_RMDUP(
        reference_cds_val
    )

    ch_versions = ch_versions.mix(SEQKIT_RMDUP.out.versions)

    Channel.of(["id": prefix]).combine(
        Channel.fromPath(params.blocks)
    ).set{ blocks_ch }

    Channel.of(["id": prefix]).combine(
        Channel.fromPath(params.blocks_tsv)
    ).set{ blocks_tsv_ch }
    
    Channel.of(["id": prefix]).combine(
        Channel.fromPath(params.overlap_blocks)
    ).set{ overlap_blocks_ch }
    
    Channel.of(["id": prefix]).combine(
        Channel.fromPath(params.overlap_blocks_tsv)
    ).set{ overlap_blocks_tsv_ch }


    assembler_contigs
    .combine(blocks_ch)
    .multiMap { it ->
        contigs: tuple(it[0], it[1])
        blocks: tuple(it[2], it[3])
    }.set{ contigs_blocks }

    MINIMAP2_MAP(
        contigs_blocks.contigs,
        contigs_blocks.blocks
    )

    ch_versions = ch_versions.mix(MINIMAP2_MAP.out.versions)

    MINIMAP2FILTER (
        MINIMAP2_MAP.out.mapping
    )

    ch_versions = ch_versions.mix(MINIMAP2_MAP.out.versions)

    sample_assemblers.map {
        assembler ->
        [["id":assembler.baseName.replace('_contigs', '_chimeric')], assembler]
    }
    .combine(overlap_blocks_ch)
    .multiMap { it ->
        contigs: tuple(it[0], it[1])
        chimeric_blocks: tuple(it[2], it[3])
    }.set{ contigs_chimeric_blocks }

    
    MINIMAP2_MAP_CHIMERIC(
        contigs_chimeric_blocks.contigs,
        contigs_chimeric_blocks.chimeric_blocks
    )

    ch_versions = ch_versions.mix(MINIMAP2_MAP_CHIMERIC.out.versions)

    MINIMAP2_FILTER_CHIMERIC (
        MINIMAP2_MAP_CHIMERIC.out.mapping
    )

    ch_versions = ch_versions.mix(MINIMAP2_FILTER_CHIMERIC.out.versions)

    MINIMAP2_FILTER_CHIMERIC.out.map.map {
        [["id":it[0]["id"].replace('_chimeric', '')], it[1]]
    }.set{ chimeric_mapping }

    MINIMAP2FILTER.out.map.join(
        chimeric_mapping
    ).set { all_mapping }

    EXTRACTMAPPEDIDS(
        all_mapping
    )

    ch_versions = ch_versions.mix(EXTRACTMAPPEDIDS.out.versions)


    sample_assemblers.map {
        assembler ->
        [["id":assembler.baseName.replace('_contigs', '')], assembler]
    }.join(
        EXTRACTMAPPEDIDS.out.mapped_ids
    ).multiMap { it ->
        contigs: tuple(it[0], it[1])
        mapped_ids: tuple(it[2])
    }.set{ contigs_mapped_ids }
    
    SEQKIT_GREP(
        contigs_mapped_ids.contigs,
        contigs_mapped_ids.mapped_ids
    )

    ch_versions = ch_versions.mix(SEQKIT_GREP.out.versions)

    SEQKIT_SEQ(
        SEQKIT_GREP.out.filter
    )

    ch_versions = ch_versions.mix(SEQKIT_SEQ.out.versions)

    BOWTIE2_BUILD(
        SEQKIT_SEQ.out.fastx
    )

    ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)

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

    ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions)

    BEDTOOLS_BAMTOBED(
        BOWTIE2_ALIGN.out.aligned
    )

    ch_versions = ch_versions.mix(BEDTOOLS_BAMTOBED.out.versions)

    BEDTOOLS_BAMTOBED.out.bed
    .map {
        meta = ["id": it[0]["assembler_id"]]
        [meta, it[1]]
    }
    .groupTuple()
    .set{ contigs_bed }

    STACKDATAFRAMES(
        contigs_bed
    )

    ch_versions = ch_versions.mix(STACKDATAFRAMES.out.versions)

    EXTRACTMAPPEDVALUES(
        MINIMAP2FILTER.out.map
    )

    ch_versions = ch_versions.mix(EXTRACTMAPPEDVALUES.out.versions)

    EXTRACTMAPPEDVALUES_CHIMERIC (
        MINIMAP2_FILTER_CHIMERIC.out.map
    )

    ch_versions = ch_versions.mix(EXTRACTMAPPEDVALUES_CHIMERIC.out.versions)

    assembler_contigs.join(
        SEQKIT_GREP.out.filter
    ).join(
        SEQKIT_SEQ.out.fastx
    ). set { contigs_filtered }

    EXTRACTIDS(
        contigs_filtered
    )

    ch_versions = ch_versions.mix(EXTRACTIDS.out.versions)

    EXTRACTIDS.out.contig_ids.join(
        EXTRACTMAPPEDVALUES.out.values
    ).join(
        EXTRACTMAPPEDVALUES_CHIMERIC.out.values.map {
            [["id":it[0]["id"].replace('_chimeric', '')], it[1]]
        }
    ).join(
        STACKDATAFRAMES.out.stacked_dfs
    ).combine(
        SEQKIT_RMDUP.out.dup_seqs.map {
            it[1]
        }
    ).combine(
        Channel.fromPath(params.gene_summary)
    ).join(
        assembler_contigs
    ).set { combined_results }

    GATHERRESULTS(
        combined_results
    )

    ch_versions = ch_versions.mix(GATHERRESULTS.out.versions)

    GATHERRESULTS.out.scores.map {
        [it[1]]
    }.collect()
    .map {
        [["id": prefix + "_all_scores"], it]
    }.set { all_scores }

    MERGEDATAFRAMESMAPPING(
        all_scores
    )

    ch_versions = ch_versions.mix(MERGEDATAFRAMESMAPPING.out.versions)

    MERGEDATAFRAMESMAPPING.out.merged_dfs.combine(
        blocks_tsv_ch.map {
            it[1]
        }
    )
    .combine(
        overlap_blocks_tsv_ch.map {
            it[1]
        }
    ).set { all_scores_with_block }

    CALCULATEFINALSCORES(
        all_scores_with_block
    )

    ch_versions = ch_versions.mix(CALCULATEFINALSCORES.out.versions)
    
    CALCULATEREFERENCEINFO(
        reference_cds_val,
        blocks_tsv_ch,
        overlap_blocks_tsv_ch
    )

    ch_versions = ch_versions.mix(CALCULATEREFERENCEINFO.out.versions)

    Channel.of(["id": prefix + "_all_contigs"])
    .concat(
        Channel.of(["id": prefix + "_unmapped_multi"]).combine(
            GATHERRESULTS.out.unmapped_contig_lengths_multi_summary.map {
                it[1]
            }
            .collect()
        )
    ).concat(
        Channel.of(["id": prefix + "_mapped_contigs"]).combine(
            GATHERRESULTS.out.mapped_contig_lengths_summary.map {
                it[1]
            }
            .collect()
        )
    ).combine(
        CALCULATEREFERENCEINFO.out.ref_stats
    ).map{
        [it[0], it[1..-1]]
    }. set { contig_summaries }

    MERGEDATAFRAMESSUMMARIES(
        contig_summaries
    )

    Channel.of(["id": prefix + "_all_contigs"]).combine(
        GATHERRESULTS.out.all_contig_lengths.map {
                it[1]
        }
        .collect()
    )
    .concat(
        Channel.of(["id": prefix + "_unmapped_single"]).combine(
            GATHERRESULTS.out.unmapped_contig_lengths_single.map {
                it[1]
            }
            .collect()
        )
    )
    .concat(
        Channel.of(["id": prefix + "_unmapped_multi"]).combine(
            GATHERRESULTS.out.unmapped_contig_lengths_multi.map {
                it[1]
            }
            .collect()
        )
    ).concat(
        Channel.of(["id": prefix + "_mapped_contigs"]).combine(
            GATHERRESULTS.out.mapped_contig_lengths.map {
                it[1]
            }
            .collect()
        )
    ).combine(
        CALCULATEREFERENCEINFO.out.cds_lenghts
    ).combine(
        CALCULATEREFERENCEINFO.out.block_lengths
    ).combine(
        CALCULATEREFERENCEINFO.out.chimeric_block_lengths
    ).map{
        [it[0], it[1..-1]]
    }.set { contig_lengths }
    
    MAKEHISTOGRAMS(
        contig_lengths
    )

    ch_versions = ch_versions.mix(MAKEHISTOGRAMS.out.versions)

    SORTDATAFRAME(
        CALCULATEFINALSCORES.out.final_scores,
        Channel.of(1)
    )

    ch_versions = ch_versions.mix(SORTDATAFRAME.out.versions)
    
    DGEEVAL(
        assembler_contigs,
        EXTRACTMAPPEDIDS.out.mapped_ids,
        reads,
        gene_summary_prefix,
        MINIMAP2FILTER.out.map,
        ch_versions
    )
    /*
    REFEVAL(
        reference_cds_val,
        reads,
        BOWTIE2_ALIGN_REFERENCE.out.aligned
    )
    */

    SORTDATAFRAME.out.dataframe.combine(
        DGEEVAL.out.eval
    ).map {
        [it[0], [it[1], it[3]]]
    }.set { all_scores_dge }

    MERGEDRESULTS(
        all_scores_dge
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

process RENAMECONTIGS {
    input:
    tuple val(meta), path(contigs)

    output:
    path "${sample}_*.fa"

    script:

    sample = meta.id
    """
    for contig in ${contigs}; do
        cp \$contig ${sample}__\${contig}
    done
    """
}

process COVERAGEVISUALIZE_ASSEMBLY_COUNTS {
    input:
    tuple val(meta), path(counts)

    output:
    path "${meta.id}_counts_histogram.png"

    script:
    //TODO: find out how to make the Rscript path relative
    """
    contig_visualiser.R ${counts}/quant.sf ${meta.id}_counts_histogram.png
    """
}

process MULTICOVERAGEVISUALISER_ASSEMBLY_COUNTS {
    input:
    tuple val(meta), path(counts)

    output:
    path "*_counts_histogram.png"

    script:
    //TODO: find out how to make the Rscript path relative

    def quant_files = counts.toList().collect { it.plus("/quant.sf")}
    def quant_names = counts.toList().collect { it.toString().split(":").last().split("_")[0]}

    def renameCommands = quant_files.indices.collect { i ->
        def file = quant_files[i]
        def name = quant_names[i]
        "cp ${file} ${name}"
    }.join("\n")

    def input_file_names = quant_names.join(" ")

    """
    ${renameCommands}

    contig_coverage_visualiser.R ${input_file_names} ${meta.id}_counts_histogram.png
    """
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
