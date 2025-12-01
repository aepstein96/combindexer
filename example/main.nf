include { BARCODE_READS as BARCODE_READS_GEX } from "./modules.nf"
include { BARCODE_READS as BARCODE_READS_TARG } from "./modules.nf"

/*
include { ALIGN_READS_AND_COUNT } from "./modules.nf"
include { ALIGN_READS } from "./modules.nf"
include { MERGE_CSV } from "./modules.nf"
include { MERGE_SAM } from "./modules.nf"
include { COUNT_TARGETED_READS_FROM_CSV } from "./modules.nf"
include { COUNT_TARGETED_READS_FROM_SAM } from "./modules.nf"
include { FILTER_TARGETED_READS } from "./modules.nf"

include { CREATE_ADATAS } from "./modules.nf"
include { MERGE_ADATAS } from "./modules.nf"
include { ADD_TARGETED_READS_TO_ADATAS } from "./modules.nf"
*/


workflow {
    main:
    samples_ch = Channel
        .fromPath(params.samples_file)
        .splitText()        .map { it.trim() }
        .filter { it }
        .filter { it != "Sample Name"}

    samples_ch_targ = samples_ch.map { "Targeted_" + it }
    samples_ch_gex = samples_ch.map { "GEX_" + it }

    barcoded_reads_gex_ch = BARCODE_READS_GEX(samples_ch_gex, Channel.fromPath(params.input_folder), Channel.fromPath(params.barcode_folder_gex))
    barcoded_reads_targ_ch = BARCODE_READS_TARG(samples_ch_targ, Channel.fromPath(params.input_folder), Channel.fromPath(params.barcode_folder_targ))

    // gene_count_ch = ALIGN_READS_AND_COUNT(barcoded_reads_gex_ch)
    // target_csv_count_ch = COUNT_TARGETED_READS_FROM_CSV(barcoded_reads_targ_ch)

    publish:
    barcode_reads_gex = barcoded_reads_gex_ch
    barcode_reads_targ = barcoded_reads_targ_ch
}

output {

    barcode_reads_gex {
        path { sample, files -> "results_gex/$sample" }
    }

    barcode_reads_targ {
        path { sample, files -> "results_targ/$sample" }
    }
}

