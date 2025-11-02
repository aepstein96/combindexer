#!/bin/bash
set -e 



# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

input_folder="data/fastq_files"
barcoded_fastq_folder="data/barcoded_fastqs"
output_folder="data/output"
barcode_folder="barcodes"
reference_folder="references"
cb_len=30
umi_len=8

# Data download sources
data_source="" # PLACEHOLDER: once final fastq files are uploaded, this will be the url.
reference_source="" # PLACEHOLDER

# Change to script directory so all relative paths are resolved correctly
cd "${SCRIPT_DIR}"

# --- Detect number of threads automatically ---
THREADS=$(nproc || sysctl -n hw.ncpu || echo 8)
export THREADS


# Download data if not already present
if [[ ! -d "${input_folder}" ]]; then
  mkdir -p "${input_folder}"
  # TODO: download data from data_source
fi

# Download references if not already present
if [[ ! -d "${reference_folder}" ]]; then
  mkdir -p "${reference_folder}"
  # TODO: download references
fi

# Barcode reads
echo "Barcoding reads..."
mkdir -p ${barcoded_fastq_folder}
python "${PROJECT_ROOT}/combindexer/readBarcoding.py" \
    --input_folder "${input_folder}" \
    --output_folder "${barcoded_fastq_folder}" \
    --barcode_folder "${barcode_folder}" \
    --cores "${THREADS}" \
    #--debug_mode

# Merge all cdna and barcodes fastqs into a single file
echo "Merging all *_cdna.fastq.gz files into all_samples_cdna.fastq.gz ..."
cat "${barcoded_fastq_folder}"/GEX*_cdna.fastq.gz > "${barcoded_fastq_folder}/all_samples_cdna.fastq.gz"

echo "Merging all *_barcodes.fastq.gz files into all_samples_barcodes.fastq.gz ..."
cat "${barcoded_fastq_folder}"/GEX*_barcodes.fastq.gz > "${barcoded_fastq_folder}/all_samples_barcodes.fastq.gz"
  

# Run STARsolo
echo "Running STARsolo..."
mkdir -p "${output_folder}"
umi_start=$((cb_len + 1))
STAR \
   --runThreadN "${THREADS}" \
  --genomeDir "${reference_folder}/hg38_star_index" \
  --readFilesIn "${barcoded_fastq_folder}/all_samples_barcodes.fastq.gz" \
                "${barcoded_fastq_folder}/all_samples_cdna.fastq.gz" \
  --readFilesCommand zcat \
  --outFileNamePrefix "${output_folder}/" \
  --soloType CB_UMI_Simple \
  --soloBarcodeMate 0 \
  --soloCBstart 1 --soloCBlen "${cb_len}" \
  --soloUMIstart "${umi_start}" --soloUMIlen "${umi_len}" \
  --soloFeatures Gene \
  --soloStrand Reverse \
  --soloUMIdedup 1MM_Directional \
  --soloCellFilter EmptyDrops_CR \
  --soloCBwhitelist None