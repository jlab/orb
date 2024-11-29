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
include { INPUT_CHECK_READS } from '../subworkflows/local/input_check_reads'
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
include { MINIMAP2_MAP                } from '../modules/local/minimap2/map/main'
include { MINIMAP2_DEDUP              } from '../modules/local/scripts/minimap2dedup/main'
include { BOWTIE2_BUILD               } from '../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN               } from '../modules/nf-core/bowtie2/align/main'
include { BEDTOOLS_BAMTOBED           } from '../modules/nf-core/bedtools/bamtobed/main'
include { SCRIPT_MAPPING_STATS        } from '../modules/local/scripts/mapping_stats/main'
include { ASSEMBLYREADSTATS           } from '../modules/local/scripts/assembly_read_stats/main'
include { MERGEMAPPINGLOGS            } from '../modules/local/scripts/merge_mapping_logs/main'
include { MERGEALLLOGS                } from '../modules/local/scripts/merge_all_logs/main'

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

    //ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    // TODO: OPTIONAL, you can use nf-validation plugin to create an input channel from the samplesheet with Channel.fromSamplesheet("input")
    // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
    // ! There is currently no tooling to help you write a sample sheet schema

    sample_assemblers = Channel.fromPath(params.contigs)
    reads = INPUT_CHECK_READS.out.reads
    reference_cds = Channel.fromPath(params.reference_cds)
    gene_summary = Channel.fromPath(params.gene_summary)

    sample_assemblers.map {
        assembler ->
        [["id":assembler.baseName.replace('_contigs', '')], assembler]
    }.set{ id_sample_assemblers }

    //sample_assemblers.view()
    //reads.view()
    /*

    sample_assemblers.map {
        assembler ->
        [["id":assembler.baseName.replace('_contigs', '')], assembler]
    }
    .combine(reference_cds)
    .multiMap { it ->
        contigs: tuple(it[0], it[1])
        reference: tuple(["id": "mt_reference"], it[2])
    }.set{ contigs_ref }
    
    MINIMAP2_MAP(
        contigs_ref.contigs,
        contigs_ref.reference
    )

    MINIMAP2_DEDUP (
        MINIMAP2_MAP.out.mapping
    )

    SALMON_INDEX(
        contigs_ref.contigs
    )
    reads.combine(SALMON_INDEX.out.index).map {
        meta = ["id": it[0]["id"] + "_" + it[2]["id"]] + ["sample_id": it[0]["id"]] + ["single_end": it[0]["single_end"]] + ["assembler": it[2]["id"]]
        [meta, it[1], it[3]]
    }.set{ reads_combined }

    SALMON_QUANT (
        reads_combined
    )
    */

    BOWTIE2_BUILD(
        id_sample_assemblers
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

    BEDTOOLS_BAMTOBED(
        BOWTIE2_ALIGN.out.aligned
    )

    //TODO call the script and cat the bed files

    BEDTOOLS_BAMTOBED.out.bed.combine(
        BOWTIE2_ALIGN.out.log, by: 0
    ).combine(
        gene_summary
    ).multiMap { it ->
        bed: tuple(it[0], it[1])
        gene_summary: tuple([:], it[3])
        log: tuple(it[0], it[2])
    }.set{ bed_gene_summary }

    SCRIPT_MAPPING_STATS(
        bed_gene_summary.bed,
        bed_gene_summary.gene_summary,
        bed_gene_summary.log
    )

    BEDTOOLS_BAMTOBED.out.bed.map {
        meta = ["id": it[0]["assembler_id"]]
        [meta, it[1]]
    }.groupTuple()
    .combine(
        gene_summary
    ).multiMap { it ->
        bed: tuple(it[0], it[1])
        gene_summary: tuple([:], it[2])
    }.set{ bed_gene_summary }


    ASSEMBLYREADSTATS(
        bed_gene_summary.bed,
        bed_gene_summary.gene_summary
    )

    SCRIPT_MAPPING_STATS.out.logs.map {
        meta = ["id": it[0]["assembler_id"]]
        [meta, it[1]]
    }.groupTuple()
    .set{ logs }

    MERGEMAPPINGLOGS (
        logs
    )

    MERGEMAPPINGLOGS.out.log_mean.map {
        meta = ["id": "all_assembler"]
        [meta, it[1]]
    }.groupTuple()
    .set { all_logs }

    MERGEALLLOGS (
        all_logs
    )

    /*
    ass_contigs = sample_assemblers.map {
        sample_id, assembler ->
        contigs = assembler.values()
        [sample_id, contigs]
    }


    renamed_contigs = RENAMECONTIGS(
        ass_contigs
    )
    .flatten()
    .map {
        contig ->
        def assembler = contig.getName().split("__")[1].minus("_contigs.fa")
        [["id": assembler], contig]
    }
    .groupTuple()

    QUAST_ASSEMBLER (
        renamed_contigs,
        [[:],[]],
        [[:],[]]
    )

    ASSEMLYLENGTHVISUALIZE_ASSEMBLY_CONTIG(renamed_contigs)

    QUAST_SAMPLE (
        ass_contigs,
        [[:],[]],
        [[:],[]]
    )

    SAMPLELENGTHVISUALIZE_ASSEMBLY_CONTIG(ass_contigs)

    all_contigs = Channel.of(["id":"all_assemblies"])
    .concat(renamed_contigs
        .map {
            it[1]
        }.collect()
    )
    .toList()

    QUAST (
        all_contigs,
        [[:],[]],
        [[:],[]]
    )

    contigs_per_sample_and_assembler = sample_assemblers
    .map {
        sample, assemblers ->
        def list = []
        sample["single_end"] = false
        for (assembler in assemblers.keySet()) {
            def sampleCopy = sample.clone()
            sampleCopy["assembler"] = assembler
            list.add([sampleCopy, assemblers[assembler]])
        }
        [list]
    }
    .flatten()
    .collate(2)

    SALMON_INDEX(
        contigs_per_sample_and_assembler
    )

    index_with_assemblername = SALMON_INDEX.out.index
        .map {
            meta, index ->
            meta_copy = meta.clone()
            ass_copy = meta.clone()
            assembler = ass_copy["assembler"]
            meta_copy.remove("assembler")
            meta_copy.remove("sample")
            [meta_copy, index, assembler]
        }

    reads_and_contigs_per_sample_and_assembler = reads.map { it ->
        def id = it[0]
        id_copy = id.clone()
        id_copy["id"] = id_copy["id"].replace("_T1", "")
        [id_copy, it[1]]
    }

    reads_combined = reads_and_contigs_per_sample_and_assembler//.view()
        .combine(index_with_assemblername, by:0)
        .map {
            meta, reads, contigs, assembler ->
            meta_copy = meta.clone()
            meta_copy["sample"] = meta_copy["id"]
            meta_copy["assembler"] = assembler
            meta_copy["id"] = meta_copy["id"].plus(":").plus(assembler)
            [meta_copy, reads, contigs]
        }


    SALMON_QUANT (
        reads_combined
    )

    ch_read_counts_grouped_by_sample = SALMON_QUANT.out.results.map {
        meta, results ->
        meta_copy = meta.clone()
        meta_copy.remove("assembler")
        meta_copy["id"] = meta["sample"]
        meta_copy.remove("sample")
        [meta_copy, results]
    }.groupTuple()

    ch_read_counts_grouped_by_sample

    MULTICOVERAGEVISUALISER_ASSEMBLY_COUNTS(ch_read_counts_grouped_by_sample)

    ch_versions = ch_versions.mix(QUAST_SAMPLE.out.versions)
    ch_versions = ch_versions.mix(SALMON_INDEX.out.versions)

    //COVERAGEVISUALIZE_ASSEMBLY_COUNTS(
    //    SALMON_QUANT.out.results
    //)


    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
    */
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
