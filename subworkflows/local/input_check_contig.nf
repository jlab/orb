//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK_CONTIGS {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    //SAMPLESHEET_CHECK ( samplesheet )
    Channel.fromPath(samplesheet)
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channel(it) }
        .set { assembler }

    emit:
    assembler                                     // channel: [ val(meta), [ reads ] ]
  //  versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample

    // add path(s) of the fastq file(s) to the meta map
    def assembler = [:]
    for (key in row.keySet()) {
        if ( key.contains("assembler:") ) {

            if (!file(row.get(key)).exists()) {
                exit 1, "ERROR: Please check input samplesheet -> in column ${key} this file does not exist! 🥲 \n${row.get(key)}"
            }
            assembler[key.minus("assembler:")] = file(row.get(key))
        }
    }

   // exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
 //   println assembler_meta
    assembler_meta = [ meta, assembler ]
    return assembler_meta
}
