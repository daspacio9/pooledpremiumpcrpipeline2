## Configuration Overview

The **Pooled Premium PCR Pipeline 2** is a Snakemake workflow designed for demultiplexing and consensus sequence generation from pooled barcoded samples submitted to the Plasmidsaurus Premium PCR service. Configuration is managed through two main files:

- **`config.yaml`** - Main workflow parameters and settings
- **`ref/*`**

The workflow performs the following key steps:

1. Prepare reference barcode-pairs file from supplied primer and barcode data (when prepare_barcode_reference with --force flag is provided)
2. Demultiplex pooled sequencing data using barcode pairs
3. Generate consensus sequences for each barcode group
4. Produce synthetic AB1 trace files and alignment reports for each barcode group
5. Generate quality statistics and coverage plots

## Configuration Parameters

### Main Pipeline Settings (`config.yaml`)

| Parameter | Type | Description | Example |
|-----------|------|-------------|---------|
| **min_depth** | integer | Minimum read depth for consensus generation | `10` |
| **input_fastq** | string | Path to input FASTQ file in `sequences/` folder | `26_03_05_DA487-492_S2R3LF_1_1.fastq.gz` |
| **primer_consensus** | string | FASTA file with primer sequences in `ref/` folder | `cole1-primer-consensus.fasta` |
| **barcode_groups** | string | CSV file defining barcode groups in `ref/` folder | `cole1-barcode-groups-48.csv` |
| **barcodes** | string | FASTA file with barcode sequences in `ref/` folder | `cole1-barcodes-48.fasta` |
| **filter_size** | integer | Minimum read length (bp) after adapter trimming | `1500` |
| **cutadapt_min_overlap** | integer | Minimum matching bases between read and adapter | `15` |
| **cutadapt_error_rate** | float | Maximum allowed error rate for adapter matching (0-1) | `0.1` |
| **medaka_model** | string | Medaka consensus polishing model | `r1041_e82_400bps_hac_v5.0.0` |
| **medaka_spoa_threads** | integer | Number of parallel threads for consensus generation | `8` |
| **debug** | boolean | Enable debug mode for verbose logging | `False` |

### Input Data Structure

#### Reference Files (in `ref/` folder)

The pipeline requires three reference files for barcode-based demultiplexing:

- **`cole1-primer-consensus.fasta`** - Consensus sequences for primer sets
- **`cole1-barcode-groups-48.csv`** - Mapping of barcodes to sample groups
- **`cole1-barcodes-48.fasta`** - FASTA sequences of all barcodes

For different barcode systems (e.g., non-colE1 plasmids), create a new working directory and prepare different reference files using the `prepare_barcode_reference` rule.

#### Sequencing Data (in `sequences/` folder)

Place your demultiplexed or pooled FASTQ files (gzipped):

```bash
sequences/
├── sample_batch_1.fastq.gz
├── sample_batch_2.fastq.gz
└── ...
```

### Advanced Parameters

| Parameter | Type | Default | Notes |
|-----------|------|---------|-------|
| **cutadapt_min_overlap** | integer | `15` | Affects barcode matching stringency. Lower values are more permissive. |
| **cutadapt_error_rate** | float | `0.1` | Allows ~1.2 mismatches in 6bp stretch, ~3 mismatches in 15bp |
| **filter_size** | integer | `1500` | Filters out adapter dimers and very short reads |
| **medaka_model** | string | `r1041_e82_400bps_hac_v5.0.0` | Must match your sequencing device; this is optimized for oxford nanopore r1041 |
| **medaka_spoa_threads** | integer | `8` | Adjust based on available CPU resources |

## Output Structure

The workflow generates outputs in the following directory structure:

```
├── ab1/                           # Synthetic AB1 trace files (ABIF format, max 5kb chunks)
│   ├── group_1_0.ab1
│   └── ...
├── consensus/                     # Consensus sequences
│   ├── group_1_consensus.fastq
│   └── ...
├── aln/                           # Alignment files
│   └── *.bam
├── consensus_split/               # Split consensus sequences by length
├── demux/                         # Demultiplexed reads
├── logs/                          # Processing logs
├── report/                        # Analysis reports
│   ├── coverage.pdf
│   ├── mismatch_freq.pdf
│   └── ...
├── sequences/                     # Processed sequences
├── demux_stats.csv                # Summary statistics table
└── consensus_summary.csv          # Consensus sequence metadata
```

### Key Output Files

| File | Description | Usage |
|------|-------------|-------|
| **\*.ab1** | Synthetic chromatogram trace files | Open in Benchling, SnapGene, APE for alignment and confidence checking |
| **group_\*_consensus.fastq** | Consensus sequences in FASTQ format | Alignment to reference plasmid |
| **consensus_summary.csv** | Metadata for all consensus sequences | Summary of results and confidence metrics |
| **coverage.pdf** | Coverage depth plot | Visualize read mapping across consensus |
| **mismatch_freq.pdf** | Mismatch frequency plot | Identify problematic regions |
| **demux_stats.csv** | Barcode demultiplexing statistics | Track read counts per barcode group |

## Usage Examples

### Running the Workflow

First, prepare the reference files (if using a new barcode system):

```bash
snakemake -s ../workflow/Snakefile -j 4 prepare_barcode_reference --use-conda --force
```

Then run the main demultiplexing and consensus workflow:

```bash
# Dry run to check for errors
snakemake -s ../workflow/Snakefile -j 4 --use-conda -np

# Execute the workflow
snakemake -s ../workflow/Snakefile -j 4 --use-conda -p
```

### Reusing a Working Directory

To restart with new input data while keeping the directory structure:

```bash
snakemake -s ../workflow/Snakefile -j 4 --use-conda clean
```

This archives the previous run with a timestamp and prepares for a new run.

## Reference Files Preparation

To use different barcode systems or primer sets, you can create custom reference files:

1. Create a new working directory:
   ```bash
   mkdir working_directory_custom_locus
   cd working_directory_custom_locus
   ```

2. Place your reference files in `ref/` folder and update `config.yaml` with new filenames

3. Run the prepare_barcode_reference step with your new files

See the main [README.md](../README.md) for more detailed workflow information.
