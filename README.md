# Combindexer

Universal processing pipeline for combinatorial indexing data with flexible barcode correction and UMI assignment.

## Overview

Combindexer is a flexible barcoding pipeline for custom single-cell genomics methods, with a special focus on combinatorial indexing. This pipeline offers several key features not found in other, similar packages:

- **Flexible Barcode Configuration**: Define multiple barcode types with custom positions, correction distances, and orientations. Uniquely define cells by combinations of barcodes. 
- **Individual barcode whitelists**: Packages such as STARsolo allow multiple barcodes, but require a whitelist that individually contains each combination of barcodes. With combinatorial indexing methods, this can be millions of combinations! Combindexer allows individual barcode whitelists to be specified for each barcode.
- **Flexible barcode correction**: 

### Key Features

- **Error Correction**: Automatic correction of sequencing errors (SNPs and indels) using Levenshtein distance
- **UMI Extraction**: Supports UMI extraction from any position in the read
- **Multi-Read Processing**: Processes multiple reads simultaneously (e.g., R1, R2, I1, I2)
- **Conditional Barcode Matching**: Supports hierarchical barcode matching with conditional next steps
- **Read Trimming**: Automatic trimming based on reverse complement matching between paired reads
- **Parallel Processing**: Multiprocessing support for processing multiple samples simultaneously
- **Flexible Output**: Can output to CSV files, FASTQ files, or both based on barcode configuration

## Installation

### Requirements

- Python 3.7+
- Required packages:
  - `fire` (for CLI interface)
  - Standard library: `collections`, `contextlib`, `multiprocessing`, `functools`, `gzip`, `glob`, `sys`, `os`

### Installation Steps

```bash
# Clone or navigate to the repository
cd combindexer

# Install dependencies (if using a virtual environment)
pip install fire
```

## Directory Structure

```
combindexer/
├── combindexer/
│   ├── readBarcoding.py      # Main barcoding pipeline
│   └── countTargetedReads.py # Read counting and statistics
├── barcodes/                 # Configuration directory
│   ├── Barcode-*.csv        # Barcode definition files
│   ├── Reads.csv            # Read configuration
│   ├── Samples.csv          # Sample list
│   ├── Settings.csv         # Pipeline settings (optional)
│   └── sheet_names.txt      # Barcode processing order (optional)
└── README.md
```

## Configuration

### 1. Barcode Definition Files (`Barcode-*.csv`)

Each barcode type requires a CSV file named `Barcode-<NAME>.csv` in the `barcodes/` directory.

**Header Format:**
```csv
Length to use,<value>,,,,,
Plates to use,<value>,,,,,
Correction distance,<value>,,,,,
Orientation,<forward|reverse>,,,,,
Correct SNPs,<True|False>,,,,,
Correct Indels,<True|False>,,,,,
,,,,,,
Barcode,Plate,Well,Condition,Read Type,Gene,Next Commands,Output File Number
```

**Fields:**
- **Length to use**: Number of bases to use for matching (`all` to use full length)
- **Plates to use**: Comma-separated plate names or `all`
- **Correction distance**: Maximum Levenshtein distance for error correction (0 = no correction)
- **Orientation**: `forward` or `reverse` (reverse complement matching)
- **Correct SNPs**: Enable/disable single nucleotide polymorphism correction
- **Correct Indels**: Enable/disable insertion/deletion correction
- **Barcode**: Full barcode sequence
- **Plate/Well/Condition**: Metadata for sample identification
- **Read Type**: Type identifier for the barcode
- **Gene**: Gene identifier (if applicable)
- **Next Commands**: Semicolon-separated commands to execute after matching this barcode
- **Output File Number**: Which output file(s) to write to (semicolon-separated numbers or `csv`)

**Special Barcodes:**
- `No match`: Default for unmatched barcodes
- `Ambiguous`: For barcodes matching multiple targets

**Example:**
```csv
Length to use,10,,,,,
Plates to use,all,,,,,
Correction distance,1,,,,,
Orientation,forward,,,,,
Correct SNPs,True,,,,,
Correct Indels,True,,,,,
,,,,,,
Barcode,Plate,Well,Condition,Read Type,Gene,Next Commands,Output File Number
GAGTTCTAAG,,A1,core,TargRT,,,
GCTCTTAGAA,,A2,core,TargRT,,,
No match,,,,Unknown,Unknown,,
```

### 2. Read Configuration (`Reads.csv`)

Defines how to process each sequencing read file.

**Format:**
```csv
Read,Input File,Commands,Trim Partner,Min Length,Output File 1,Output File 2,Output File 3
```

**Fields:**
- **Read**: Read identifier (e.g., R1, R2, I1, I2)
- **Input File**: Glob pattern for input FASTQ files (use `{SAMPLE}` placeholder)
- **Commands**: Semicolon-separated commands defining barcode extraction:
  - `UMI:<length>`: Extract UMI of specified length
  - `Bar:<name>`: Match barcode from `Barcode-<name>.csv`
  - `Goto:<position>`: Jump to absolute position in read
  - `Gap:<length>`: Skip forward by specified length
  - `Trim`: Mark position for trimming (uses position after barcode extraction)
- **Trim Partner**: Read name to use for reverse complement trimming
- **Min Length**: Minimum read length after trimming
- **Output File N**: Output file patterns (use `{SAMPLE}` placeholder)

**Example:**
```csv
Read,Input File,Commands,Trim Partner,Min Length,Output File 1,Output File 2,Output File 3
R1,{SAMPLE}_*R1*.fastq.gz,UMI:8;Bar:RT;Gap:30;Bar:RTPrimer,,,
R2,{SAMPLE}_*R3*.fastq.gz,Bar:SSSprimer,R1,20,{SAMPLE}.R2.fastq.gz,,
I5,{SAMPLE}_*R2*.fastq.gz,Bar:P5,,,,
I7,{SAMPLE}_*I1*.fastq.gz,Bar:P7,,,,
```

### 3. Sample List (`Samples.csv`)

Simple CSV file listing sample names:

```csv
Sample Name
Targeted_E2
Targeted_F2
```

### 4. Optional Configuration

**`sheet_names.txt`**: Define custom barcode processing order:
```
Barcode-RT,Barcode-P5,Barcode-P7,Barcode-Gene
```

**`Settings.csv`**: Optional settings file (currently not directly used by the pipeline but may be referenced for other tools):
```csv
Date,2023-04-02 00:00:00
Experiment Name,Joint P53 sequencing
Input Folder,/path/to/input
Output Folder,/path/to/output
```

## Usage

### Basic Usage

```bash
python combindexer/readBarcoding.py \
    --input_folder /path/to/input/fastq \
    --output_folder /path/to/output \
    --barcode_folder /path/to/barcodes \
    --cores 24 \
    --default_write_setting "1;csv"
```

### Parameters

- `input_folder`: Directory containing input FASTQ files
- `output_folder`: Directory for output files
- `barcode_folder`: Directory containing configuration files (`Reads.csv`, `Samples.csv`, `Barcode-*.csv`)
- `cores`: Number of CPU cores for parallel processing
- `default_write_setting`: Default output targets (semicolon-separated: file numbers or `csv`)
- `debug_mode`: Enable debug output (optional, defaults to False)

### Output Files

For each sample, the pipeline generates:

1. **CSV Output** (`{SAMPLE}_output.csv`): Contains metadata and read sequences
   - Columns: `wells_conditions,UMI,read_types,<read1>,<read2>,...`
   - `wells_conditions`: Plate-well pairs and conditions (semicolon-separated, joined by `&`)
   - `UMI`: Concatenated UMI sequences from all reads
   - `read_types`: Read type identifiers from barcodes
   - Subsequent columns: Sequences for each read

2. **FASTQ Files**: If specified in `Reads.csv`, individual FASTQ files with demultiplexed reads
   - Header format: `@<wells_conditions>,<UMI>,<read_types>,<original_read_id>`

3. **`reads_to_align.txt`**: List of output file patterns for downstream alignment

### Read Counting

Use `countTargetedReads.py` to generate read statistics:

```bash
python combindexer/countTargetedReads.py \
    <sample_name> \
    <input_folder> \
    <barcoded_folder> \
    <umis_folder> \
    <output_file>
```

**Parameters:**
- `sample_name`: Sample identifier
- `input_folder`: Original FASTQ input folder
- `barcoded_folder`: Folder containing `{SAMPLE}_output.csv`
- `umis_folder`: Folder containing UMI count files (`{SAMPLE}_output_umi_counts.csv`)
- `output_file`: Output CSV file for statistics

**Statistics Generated:**
- Raw reads
- Barcoded reads
- Assigned UMI reads
- Gene identified reads
- Gene identified with unique UMI counts (>=1, >=2, >=3, >=4, >=5)

## Command Reference

### Barcode Matching Commands

- **`Bar:<name>`**: Match barcode from `Barcode-<name>.csv`
- **`UMI:<length>`**: Extract UMI sequence of specified length
- **`Goto:<position>`**: Jump to absolute base position in read
- **`Gap:<length>`**: Skip forward by specified number of bases
- **`Trim`**: Mark current position for trimming

### Barcode Error Correction

The pipeline uses Levenshtein distance for error correction:

- **SNPs**: Single nucleotide substitutions
- **Indels**: Insertions and deletions (with post-context compensation)
- **Distance**: Configurable maximum edit distance (typically 0-2)
- **Ambiguity Handling**: Ambiguous matches (matching multiple barcodes) can be marked or excluded

### Conditional Barcode Matching

Barcodes can specify `Next Commands` to modify subsequent processing:
- After matching a barcode, pipeline can jump to different positions
- Supports hierarchical matching (e.g., match primer, then gene-specific barcode)

**Example Flow:**
1. Match RT primer barcode
2. Skip gap (e.g., `Gap:30`)
3. Match gene-specific RTPrimer barcode
4. Conditional commands may modify further processing

### Read Trimming

- Reads can be trimmed based on reverse complement of partner read
- Uses `Trim Partner` field in `Reads.csv`
- Minimum length after trimming can be specified
- Trimming searches for reverse complement sequence in the read

## Workflow Example

1. **Setup Configuration:**
   - Define barcode files for each barcode type (P5, P7, RT, Gene, etc.)
   - Configure `Reads.csv` with read processing commands
   - List samples in `Samples.csv`

2. **Run Pipeline:**
   ```bash
   python combindexer/readBarcoding.py \
       --input_folder ./fastq_files \
       --output_folder ./processed \
       --barcode_folder ./barcodes \
       --cores 24
   ```

3. **Process Output:**
   - CSV files contain barcode assignments and sequences
   - FASTQ files (if configured) are ready for alignment
   - Use `countTargetedReads.py` for quality metrics

4. **Downstream Analysis:**
   - Align reads using `reads_to_align.txt` as reference
   - Use UMI counts for deduplication
   - Analyze barcode combinations for sample identification

## Advanced Features

### Output File Routing

Barcodes can specify which output files to write to:
- Include `Output File Number` in barcode CSV
- Use semicolon-separated list: `1;2;csv`
- Pipeline writes to specified files based on barcode match
- Default: `default_write_setting` parameter (typically `"1;csv"`)

### Multi-Plate Support

- Filter barcodes by plate using `Plates to use` in barcode CSV
- Combine barcode information across multiple plates
- Plate information included in output metadata

### Debug Mode

Enable detailed logging:
```bash
python combindexer/readBarcoding.py ... --debug_mode True
```

Shows:
- Barcode matching decisions
- Sequence extraction details
- File paths and patterns

## Troubleshooting

### Common Issues

1. **No fastq files found**: Check `Input File` patterns in `Reads.csv` match your file naming
2. **Barcode matching fails**: Verify `Length to use` and `Correction distance` settings
3. **Ambiguous matches**: Review barcode sequences for similarity, adjust correction distance
4. **Trimming issues**: Check `Trim Partner` configuration and minimum length requirements

### File Format Requirements

- Input FASTQ files must be valid (4 lines per read)
- Supports both compressed (`.gz`) and uncompressed FASTQ
- CSV files must use comma delimiters
- Barcode sequences should not contain ambiguous bases (N) unless intentional

## Citation

If you use Combindexer in your research, please cite appropriately.

## License

[Add license information if applicable]
