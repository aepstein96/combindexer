process BARCODE_READS {
    tag "$sample"

    input:
    val(sample)
    path(input_folder)
    path(barcode_folder)

    output:
    tuple val(sample), path(output_folder)

    script:
    output_folder = "reads_barcoded"
    
    """
    python /ru-auth/local/home/aepstein/cursor_projects/combindexer/combindexer/readBarcoding_single.py --sample "${sample}" \
        --input_folder "${input_folder}" \
        --output_folder "${output_folder}" \
        --barcode_folder "${barcode_folder}"
    """
}

process ALIGN_READS {
    tag "$sample"

    input:
    tuple val(sample), val(cb_len), val(umi_len), path(barcoded_fastq_folder)

    output:
    tuple val(sample), path(output_folder)

    script:
    output_folder = "reads_aligned"

    """
    STAR \
    --runThreadN "${THREADS}" \
    --genomeDir "${reference_folder}/hg38_star_index" \
    --readFilesIn "${barcoded_fastq_folder}/all_samples_cdna.fastq.gz" \
                    "${barcoded_fastq_folder}/all_samples_barcodes.fastq.gz" \
    --readFilesCommand zcat \
    --outFileNamePrefix "${output_folder_genefull}/" \
    --soloType CB_UMI_Simple \
    --soloBarcodeMate 0 \
    --soloCBstart 1 --soloCBlen "${cb_len}" \
    --soloUMIstart "${umi_start}" --soloUMIlen "${umi_len}" \
    --soloFeatures GeneFull \
    --soloStrand Unstranded \
    --outSAMstrandField intronMotif \
    --soloUMIdedup 1MM_Directional \
    --soloCellFilter EmptyDrops_CR \
    --soloCBwhitelist None
    """

}

process CREATE_ADATAS {

    input:
    tuple val(sample), path(output_folder)

    output:
    tuple val(sample), path(output_folder)

    script:
    output_folder = "reads_merged"

    """
    
    """

}

process MERGE_ADATAS {

    input:
    path(input_folder)

    output:
    path(output_folder)

    script:
    output_folder = "adatas_merged"

    """
    
    """
}