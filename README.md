# Snakemake workflow: `pooledpremiumpcrpipeline2`


[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/<owner>/<repo>/workflows/Tests/badge.svg?branch=main)](https://github.com/<owner>/<repo>/actions?query=branch%3Amain+workflow%3ATests)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![workflow catalog](https://img.shields.io/badge/Snakemake%20workflow%20catalog-darkgreen)](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/<owner>/<repo>)

A Snakemake workflow for `demultiplex and consensus sequence generation after pooling barcoded samples for the plasmidsaurus premium pcr service`

- [Snakemake workflow: `pooledpremiumpcrpipeline2`](#snakemake-workflow-name)
  - [Usage](#usage)
  - [Deployment options](#deployment-options)
  - [Authors](#authors)
  - [References](#references)
  - [TODO](#todo)

## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/<owner>/<repo>).

Detailed information about input data and workflow configuration can also be found in the [`config/README.md`](config/README.md). 

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this repository or its DOI.

## Deployment options

To run the workflow from command line, change the working directory.

```bash
cd path/to/snakemake-workflow-name
git clone https://github.com/daspacio9/pooledpremiumpcrpipeline2.git

```
Enter your git credentials

Create the pooled premium pcr pipeline 2 conda environment (pppp) with all the necessary dependencies:

```bash
cd pooledpremiumpcrpipeline2

```
Build the conda environment, activate, change to the example
```bash
conda env create --file workflow/envs/environment.yaml

conda activate pppp

cd .test/

```
Adjust options in the default config file `config/config.yaml`. Always change to your sequences name.
Before running the complete workflow, you can perform a dry run using:

## Running the Example: prepare_reference
Run the pipeline with a dryrun (-np), then run without to execute

The pipeline runs with a separate workflow to prepare the barcode pairs reference file from the supplied files. This allows other indexing primers to be used in the pipeline. By providing a different barcode groups csv, barcode sequence fasta file, and primer consensus file
```bash
snakemake -s ../workflow/Snakefile -j -4 prepare_reference
```
## Running the Example: demultiplex and consensus
The barcode-pairs reference file should now exist. Start by using the touch command to change the system timestamp of the sequence files you downloaded from a potentially different timezone.

```bash
touch sequences/your_sequences.fastq

gzip sequences/your_sequences.fastq

snakemake -s ../workflow/Snakefile -j 4 -np

snakemake -s ../workflow/Snakefile -j 4 -p
```

Replace "../workflow/Snakefile" with the actual path to the Snakefile where it is installed. When running the Example this is one folder up.

## Multiple Runs in the Same Directory: clean
To reuse a directory for another sequencing batch, use the clean rule:
```bash
snakemake -s ../workflow/Snakefile -j 4 clean
```
This neatly archives a run with the date and time in the archive name and prepares for a new run. 

## Key Outputs for Analysis
The key outputs which are listed below. The report plots are useful for investigating the confidence of a consensus sequence which can help in interpreting the sequencing result. Each example_consensus.fastq can be aligned to the expected plasmid sequence. consensus_summary.csv is a table of all consensus sequences generated. example coverage.pdf shows the number of reads mapped to each base present in the constructed consensus sequence (coverage) and the number of basecalls in the subreads that match the chosen consensus basecall in the consensus sequence. example_mismatch_freq.pdf plots the only the freq of mismatches to the consensus basecall. 





## Working Directory (needs to be updated for pppp2)
The working directory is organized as follows. The required files are the inputs. If barcode-pairs.txt is not present, then prepare_reference needs to be run first.

```bash
│──────────────────────────────────────────────────────────
│                                               KEY OUTPUTS
├── demux_stats.csv
├── consensus
│   ├── bulk_consensus
│   │   ├── consensus.fastq
│   │   ├── poa.fasta
│   │   ├── poa.fasta.fai
│   │   ├── subreads_to_spoa.bam
│   │   └── subreads_to_spoa.bam.bai
│   └── consensus_summary.csv
├── consensus_split
│   ├── example_consensus.fastq
│   ├── example_consensus.fastq.fai
│   └── .workflow_flag.done
├── report
│   ├── example_coverage.pdf
│   └── example_mismatch_freq.pdf
```

```bash
├── config.yaml
├── ref
│   ├── cole1-barcodes-48.fasta
│   ├── cole1-primer-consensus.fasta
│   └── cole1-barcode-groups-48.csv
├── sequences
│   └── R8L89G-1-pDA260-269.fastq.gz
│                                                    INPUTS
│──────────────────────────────────────────────────────────
│                                                   OUTPUTS
├── ref
│   └── barcode-pairs.txt
├── aln
│   ├── example.sam
│   ├── example.bam
│   └── example.bam.bai
├── demux_stats.csv
├── demux
│   ├── low_depth
│   │   └── example_discarded.fastq.gz
│   ├── example.fastq.gz
│   └── unknown.fastq.gz
├── consensus
│   ├── example_consensus
│   │   ├── subreads_to_spoa.bam
│   │   ├── subreads_to_spoa.bam.bai
│   │   ├── consensus.hdf
│   │   ├── poa.fasta.fai
│   │   └── poa.fasta
│   ├── example_consensus.fastq
│   ├── example_consensus.fastq.fai
│   └── consensus_summary.csv
├── report
│   ├── example_coverage.pdf
│   ├── example_coverage.txt
│   ├── example_mismatch_freq.pdf
│   ├── example_mpileup.txt 
│   └── example_pypileup.txt
├── logs
│   ├── aln
│   ├── cutadapt
│   └── consensus
```







## Authors

- Derek Aspacio
  - UCI
  - ORCID profile
  - home page

## References

> Köster, J., Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., & Nahnsen, S. _Sustainable data analysis with Snakemake_. F1000Research, 10:33, 10, 33, **2021**. https://doi.org/10.12688/f1000research.29032.2.

## TODO

- Replace `<owner>` and `<repo>` everywhere in the template with the correct user name/organization, and the repository name. The workflow will be automatically added to the [snakemake workflow catalog](https://snakemake.github.io/snakemake-workflow-catalog/index.html) once it is publicly available on Github.
- Replace `<name>` with the workflow name (can be the same as `<repo>`).
- Replace `<description>` with a description of what the workflow does.
- Update the [deployment](#deployment-options), [authors](#authors) and [references](#references) sections.
- Update the `README.md` badges. Add or remove badges for `conda`/`singularity`/`apptainer` usage depending on the workflow's [deployment](#deployment-options) options.
- Do not forget to also adjust the configuration-specific `config/README.md` file.
