# Snakemake workflow: `pooledpremiumpcrpipeline2 2.0.0`


[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/daspacio9/pooledpremiumpcrpipeline2/workflows/Tests/badge.svg?branch=main)](https://github.com/daspacio9/pooledpremiumpcrpipeline2/actions?query=branch%3Amain+workflow%3ATests)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![workflow catalog](https://img.shields.io/badge/Snakemake%20workflow%20catalog-darkgreen)](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/daspacio9/pooledpremiumpcrpipeline2)

A Snakemake workflow for demultiplex and consensus sequence generation after pooling barcoded samples for the plasmidsaurus premium pcr service

- [Snakemake workflow: `pooledpremiumpcrpipeline2`](#snakemake-workflow-name)
  - [Usage](#usage)
  - [Deployment options: 'conda'](#deployment-options)
  - [Authors: 'Derek Sprayberry Aspacio'](#authors)
  - [References](#references)

## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/daspacio9/pooledpremiumpcrpipeline2).
Detailed information about input data and workflow configuration can also be found in the [`config/README.md`](config/README.md). 
If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this repository or its DOI.

## Deployment options
Use conda
To run the workflow from command line, navigate to the folder where you want to clone the repository.
```bash
cd path/to/location/on/your/desktop
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
Run the pipeline with a dryrun (-np), then run without to execute.
The pipeline runs with a separate workflow to prepare the barcode pairs reference file from the supplied files. This allows other indexing primers to be used in the pipeline. By providing a different barcode groups csv, barcode sequence fasta file, and primer consensus file
```bash
snakemake -s ../workflow/Snakefile -j -4 prepare_barcode_reference --use-conda -np 
snakemake -s ../workflow/Snakefile -j -4 prepare_barcode_reference  --use-conda --force
```
## Running the Example: demultiplex and consensus
The barcode-pairs reference file should always exist. Change it by using the prepare_reference module above with --force to overwrite the existing colE1 version of the reference file. Start by using the touch command to change the system timestamp of the sequence files you downloaded from a potentially different timezone.

```bash
touch sequences/your_sequences.fastq
gzip sequences/your_sequences.fastq
snakemake -s ../workflow/Snakefile -j 4 --use-conda -np 
snakemake -s ../workflow/Snakefile -j 4 --use-conda -p
```
Replace "../workflow/Snakefile" with the actual path to the Snakefile where it is installed. When running the Example this is one directory up.

There is a template working directory named "working_directory/". I typically run all samples that use the same primer set in this directory. Use the clean function below to prepare the directory for a new run. When I create a new primer set for non-colE1 plasmid sequencing, I use a new working directory with the same structure but different files in the ref/ folder. For example, "working_directory_ura3locus" for multiplex sequencing of yeast genomic integrations at the ura 3 locus. 

## Multiple Runs in the Same Directory: clean
To reuse a directory for another sequencing batch, use the clean rule:
```bash
snakemake -s ../workflow/Snakefile -j 4 --use-conda clean
```
This neatly archives a run with the date and time in the archive name and prepares for a new run.

## Key Outputs for Analysis
The key outputs which are listed below. The report plots are useful for investigating the confidence of a consensus sequence which can help in interpreting the sequencing result. 

  example_group.ab1 files can be aligned to the expected sequence in Benchling, SnapGene, APE, etc. These are generated synthetically based on the number of raw base counts for each base present in the consensus sequence. File size limits of ~8 kb is why the ab1 sequences are limited to 5 kb "chunks". Each base count is normalized to the most frequent base count and then used to generate a synthetic chromatogram trace for each base with height proportional to frequence of base counts. Useful for checking confidence when aligning to the expected reference sequence. 
  
  Each example_group_consensus.fastq can be aligned to the expected plasmid sequence. 
  
  consensus_summary.csv is a table of all consensus sequences generated. 
  
  example coverage.pdf shows the number of reads mapped to each base present in the constructed consensus sequence (coverage) and the number of basecalls in the subreads that match the chosen consensus basecall in the consensus sequence. 
  
  example_mismatch_freq.pdf plots the only the freq of mismatches to the consensus basecall. 

  demux_stats.csv lists the number of subreads used to generate each consensus.

## Working Directory
The working directory is organized as follows. The required files are the inputs. If barcode-pairs.txt is not present, then prepare_reference needs to be run first.

```bash
│──────────────────────────────────────────────────────────
│                                               KEY OUTPUTS
├── demux_stats.csv
├── ab1
│   └── example_group_0.ab1
├── consensus
│   ├── bulk_consensus
│   │   ├── consensus.fastq
│   │   ├── poa.fasta
│   │   ├── poa.fasta.fai
│   │   ├── subreads_to_spoa.bam
│   │   └── subreads_to_spoa.bam.bai
│   └── consensus_summary.csv
├── consensus_split
│   ├── example_group_consensus.fastq
│   ├── example_group_consensus.fastq.fai
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
  - University of California, Irvine
  - ORCID: 0000-0002-8210-1811

## References

> Köster, J., Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., & Nahnsen, S. _Sustainable data analysis with Snakemake_. F1000Research, 10:33, 10, 33, **2021**. https://doi.org/10.12688/f1000research.29032.2.

> Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34:3094-3100. doi:10.1093/bioinformatics/bty191

> Medaka: Consensus generation for long reads. Oxford Nanopore Technologies. https://github.com/nanoporetech/medaka

> Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. 
EMBnet Journal, 17(1):10-12.

> Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., ... & Durbin, R. (2009). The sequence alignment/map format and SAMtools. Bioinformatics, 25(16):2078-2079.

> Cock, P. J., Chang, J. T., Chang, K. M., et al. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25(3):422-423.

> Shen, W., Le, S., Li, Y., & Hu, F. (2016). SeqKit: A cross-platform and ultrafast toolkit for FASTA/Q file manipulation. PLoS ONE, 11(10):e0163962.
