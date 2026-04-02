# \HEADER\-------------------------------------------------------------------------
#
#  CONTENTS      : Snakemake nanopore data pipeline
#
#  DESCRIPTION   : demultiplexing of nanopore reads by end barcodes using cutadapt package
#
#  RESTRICTIONS  : none
#
#  REQUIRES      : none
#
# ---------------------------------------------------------------------------------


### Subread classification and consensus generation rules
# Classify reads in fastq to subreads following the required formatting by medaka smolecule package
# -----------------------------------------------------
rule fastq_to_fastq_subreads:
    input:
        inf="demux/{sample}.fastq.gz",
    output:
        outf="demux/temp/{sample}_subreads_batch.fastq.gz",
    log:
        logf="logs/consensus/{sample}_fastq_to_fastq_subreads.log",
    conda:
        "../envs/consensus.yaml"
    script:
        "../scripts/fastq_to_fastq_subreads.py"


# Rule to generate consensus sequences using Medaka from subreads
# Helper to get all samples from your data folder
# -----------------------------------------------------
rule medaka_consensus_from_subreads:
    input:
        # This collects ALL fastq.gz files into one list for parallel threading on Medaka smolecule
        fastas=expand(
            "demux/temp/{sample}_subreads_batch.fastq.gz", sample=get_passed_samples
        ),
        check=directory("demux/filtered_list")  # Require checkpoint output
    output:
        outDir=(directory("consensus/bulk_consensus")),
        consensus="consensus/bulk_consensus/consensus.fastq",
    log:
        "logs/consensus/bulk_consensus.log",
    conda:
        "../envs/consensus.yaml"
    threads: config["medaka_spoa_threads"]
    params:
        min_depth=(config["min_depth"] - 1),
        medaka_model=config["medaka_model"],
    shell:
        r"""
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting medaka consensus generation" >> {log}
        rm -rf {output.outDir} >> {log} 2>&1
        
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running medaka smolecule" >> {log}
        medaka smolecule \
            {output.outDir} \
            {input.fastas} \
            --method spoa \
            --threads {threads} \
            --depth {params.min_depth} \
            --qualities \
            --batch_size 20 \
            --model {params.medaka_model} \
            --spoa_min_coverage {params.min_depth} \
        >> {log} 2>&1
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Medaka consensus generation complete" >> {log}
        """


split_dir = "consensus_split"


# Rule to split the parallel process fastq back into individual files
# -----------------------------------------------------
checkpoint split_fastq:
    input:
        batchconsensus="consensus/bulk_consensus/consensus.fastq",
    output:
        out_dir=directory(split_dir),
    log:
        logf="logs/consensus/split_fastq.log",
    conda:
        "../envs/consensus.yaml"
    script:
        "../scripts/split_fastq.py"


# Rule to create a summary CSV with sample, sequence, and quality score columns from all consensus FASTA files
# -----------------------------------------------------
rule consensus_summary_csv:
    input:
        # consensuses = get_split_files
        consensuses=lambda wildcards: get_splits(
            wildcards, prefix="", suffix="_consensus.fastq"
        ),
    output:
        outf=report("consensus/consensus_summary.csv", category="consensus"),
    log:
        logf="logs/consensus/consensus_summary_csv.log",
    conda:
        "../envs/consensus.yaml"
    script:
        "../scripts/consensus_summary_csv.py"
