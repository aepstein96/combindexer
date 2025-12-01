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
#python "${PROJECT_ROOT}/combindexer/readBarcoding.py" \
#    --input_folder "${input_folder}" \
#    --output_folder "${barcoded_fastq_folder}" \
#    --barcode_folder "${barcode_folder}" \
#    --cores "${THREADS}" \
    #--debug_mode

# Merge all cdna and barcodes fastqs into a single file
echo "Merging all *_cdna.fastq.gz files into all_samples_cdna.fastq.gz ..."
#cat "${barcoded_fastq_folder}"/GEX*_cdna.fastq.gz > "${barcoded_fastq_folder}/all_samples_cdna.fastq.gz"

echo "Merging all *_barcodes.fastq.gz files into all_samples_barcodes.fastq.gz ..."
cat "${barcoded_fastq_folder}"/GEX*_barcodes.fastq.gz > "${barcoded_fastq_folder}/all_samples_barcodes.fastq.gz"
  

# Run STARsolo - Test 1: GeneFull alone
echo "=========================================="
echo "Running STARsolo Test 1: GeneFull alone"
echo "=========================================="
mkdir -p "${output_folder_genefull}"
umi_start=$((cb_len + 1))
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

echo "Test 1 (GeneFull) completed successfully"

# Run STARsolo - Test 2: GeneFull + winFlankNbins
echo "=========================================="
echo "Running STARsolo Test 2: GeneFull + winFlankNbins 6"
echo "=========================================="
mkdir -p "${output_folder_genefull_flank}"
umi_start=$((cb_len + 1))
STAR \
   --runThreadN "${THREADS}" \
  --genomeDir "${reference_folder}/hg38_star_index" \
  --readFilesIn "${barcoded_fastq_folder}/all_samples_cdna.fastq.gz" \
                "${barcoded_fastq_folder}/all_samples_barcodes.fastq.gz" \
  --readFilesCommand zcat \
  --outFileNamePrefix "${output_folder_genefull_flank}/" \
  --soloType CB_UMI_Simple \
  --soloBarcodeMate 0 \
  --soloCBstart 1 --soloCBlen "${cb_len}" \
  --soloUMIstart "${umi_start}" --soloUMIlen "${umi_len}" \
  --soloFeatures GeneFull \
  --soloStrand Unstranded \
  --outSAMstrandField intronMotif \
  --winFlankNbins 6 \
  --soloUMIdedup 1MM_Directional \
  --soloCellFilter EmptyDrops_CR \
  --soloCBwhitelist None

echo "Test 2 (GeneFull + winFlankNbins) completed successfully"

# Comparison analysis
echo ""
echo "=========================================="
echo "COMPARING RESULTS: GeneFull vs GeneFull + winFlankNbins"
echo "=========================================="

python3 << 'PYTHON_EOF'
import sys
import os

def read_summary(filepath):
    """Read STARsolo Summary.csv and extract key metrics"""
    metrics = {}
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if ',' in line:
                    key, value = line.strip().split(',', 1)
                    metrics[key] = value
    except Exception as e:
        print(f"Error reading {filepath}: {e}", file=sys.stderr)
        return None
    return metrics

def read_features_stats(filepath):
    """Read Features.stats and extract key counts"""
    stats = {}
    try:
        with open(filepath, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    key = parts[0].strip()
                    try:
                        value = int(parts[-1])
                        stats[key] = value
                    except:
                        pass
    except Exception as e:
        print(f"Error reading {filepath}: {e}", file=sys.stderr)
        return None
    return stats

# Paths
genefull_summary = "data/output_genefull/Solo.out/GeneFull/Summary.csv"
genefull_features = "data/output_genefull/Solo.out/GeneFull/Features.stats"
genefull_flank_summary = "data/output_genefull_flank/Solo.out/GeneFull/Summary.csv"
genefull_flank_features = "data/output_genefull_flank/Solo.out/GeneFull/Features.stats"

print("\n" + "="*70)
print("COMPARISON: GeneFull vs GeneFull + winFlankNbins 6")
print("="*70 + "\n")

# Read data
gf_summary = read_summary(genefull_summary)
gf_features = read_features_stats(genefull_features)
gff_summary = read_summary(genefull_flank_summary)
gff_features = read_features_stats(genefull_flank_features)

if not all([gf_summary, gf_features, gff_summary, gff_features]):
    print("ERROR: Could not read all output files. Check paths.")
    sys.exit(1)

# Compare key metrics
comparisons = [
    ("Gene mapping rate", "Reads Mapped to Gene: Unique Gene", "%"),
    ("noNoFeature reads", "noNoFeature", "reads"),
    ("yesWLmatch reads", "yesWLmatch", "reads"),
    ("Mean reads/cell", "Mean Reads per Cell", "reads"),
    ("Mean UMI/cell", "Mean UMI per Cell", "UMIs"),
    ("Mean genes/cell", "Mean Gene per Cell", "genes"),
    ("Total genes detected", "Total Gene Detected", "genes"),
    ("Estimated cells", "Estimated Number of Cells", "cells"),
]

print(f"{'Metric':<30} {'GeneFull':<20} {'GeneFull+Flank':<20} {'Improvement':<15}")
print("-"*85)

for label, key, unit in comparisons:
    if key in gf_summary:
        try:
            gf_val = float(gf_summary[key])
            gff_val = float(gff_summary[key])
            improvement = f"+{((gff_val-gf_val)/gf_val*100):.1f}%"
            print(f"{label:<30} {gf_val:<20.2f} {gff_val:<20.2f} {improvement:<15}")
        except:
            pass
    elif key in gf_features:
        gf_val = gf_features[key]
        gff_val = gff_features[key]
        if "noNoFeature" in label:
            improvement = f"-{((gf_val-gff_val)/gf_val*100):.1f}%"
        else:
            improvement = f"+{((gff_val-gf_val)/gf_val*100):.1f}%"
        print(f"{label:<30} {gf_val:<20,} {gff_val:<20,} {improvement:<15}")

print("\n" + "="*70)
print("KEY INSIGHTS:")
print("="*70)

# Calculate key improvements
try:
    gf_map_rate = float(gf_summary.get("Reads Mapped to Gene: Unique Gene", 0))
    gff_map_rate = float(gff_summary.get("Reads Mapped to Gene: Unique Gene", 0))
    map_improvement = ((gff_map_rate - gf_map_rate) / gf_map_rate * 100) if gf_map_rate > 0 else 0
    
    gf_nofeature = gf_features.get("noNoFeature", 0)
    gff_nofeature = gff_features.get("noNoFeature", 0)
    nofeature_reduction = ((gf_nofeature - gff_nofeature) / gf_nofeature * 100) if gf_nofeature > 0 else 0
    
    print(f"✓ Gene mapping rate: {gf_map_rate:.2f}% → {gff_map_rate:.2f}% ({map_improvement:+.1f}%)")
    print(f"✓ noNoFeature reads: {gf_nofeature:,} → {gff_nofeature:,} ({nofeature_reduction:+.1f}% reduction)")
    print(f"✓ Gene-matched reads: {gf_features.get('yesWLmatch', 0):,} → {gff_features.get('yesWLmatch', 0):,}")
    
    if map_improvement > 0:
        print(f"\n✓ winFlankNbins 6 IMPROVED gene mapping by {map_improvement:.1f}%")
    else:
        print(f"\n⚠ winFlankNbins 6 had minimal impact ({map_improvement:+.1f}%)")
        
except Exception as e:
    print(f"Error calculating improvements: {e}")

print("="*70 + "\n")
PYTHON_EOF

echo ""
echo "=========================================="
echo "Comparison complete! Check results above."
echo "=========================================="