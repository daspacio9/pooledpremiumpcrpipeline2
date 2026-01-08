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



# Function to get passed samples after demux checkpoint is complete
# -----------------------------------------------------
def get_passed_samples(wildcards): #wildcards is required argument for all snakemake dynamic input/output objects to be accessed in the function
# This function reads the demux_stats.csv file to get samples that passed the depth filter
    passed_samples = []
    demux_stats_file = checkpoints.demux_stats.get().output[0] # get the output file of the checkpoint
    with open(demux_stats_file) as f:
        next(f)  # skip header
        for line in f:
            sample, n_reads = line.strip().split(",")
            if sample != "unknown" and int(n_reads) >= config["min_depth"]:
                passed_samples.append(sample)
    return passed_samples


### Subread classification and consensus generation rules
# Classify reads in fastq to subreads following the required formatting by medaka smolecule package
# -----------------------------------------------------
rule fastq_to_fastq_subreads:
    input:
        "demux/{sample}.fastq.gz"
    output:
        "demux/temp/{sample}_subreads_batch.fastq.gz"
    run:
        from Bio import SeqIO
        import gzip
        records = []
        with gzip.open(input[0], "rt") as fin:
            for i, record in enumerate(SeqIO.parse(fin, "fastq")):
                rec_id = f"{wildcards.sample}_subread{i+1}"
                record.id = rec_id
                record.description = ""
                records.append(record)
        with open(output[0], "w") as fout:
            SeqIO.write(records, fout, "fastq")
        


# Rule to generate consensus sequences using Medaka from subreads
# Helper to get all samples from your data folder
# -----------------------------------------------------
rule medaka_consensus_from_subreads:
    input:
        # This collects ALL fastq.gz files into one list for parallel threading on Medaka smolecule
        fastas = expand("demux/temp/{sample}_subreads_batch.fastq.gz", sample=get_passed_samples)
    output:
        outDir = (directory("consensus/bulk_consensus")),
        consensus = "consensus/bulk_consensus/consensus.fastq"
    log:
        "logs/consensus/bulk_consensus.log"
    threads: config["medaka_spoa_threads"]
    shell:
        r"""
        rm -rf {output.outDir}
        medaka smolecule \
            {output.outDir} \
            {input.fastas} \
            --method spoa \
            --threads {threads} \
            --depth {config[min_depth]} \
            --qualities \
            --batch_size 20 \
            --model {config[medaka_model]} \
            --spoa_min_coverage {config[min_depth]} \
        > {log} 2>&1
        """

split_dir = "consensus_split"
# Rule to split the parallel process fastq back into individual files
# -----------------------------------------------------
checkpoint split_fastq:
    input:
        "consensus/bulk_consensus/consensus.fastq"
    output:
        out_dir = directory(split_dir)
    run:
        import os
        from pathlib import Path
        from Bio import SeqIO # Requires biopython
        Path(output.out_dir).mkdir(parents=True, exist_ok=True)
        
        with open(input[0], "r") as handle:
            for record in SeqIO.parse(handle, "fastq"):
                # Clean up the ID: Benchling doesn't like special chars in names
                # We strip common sequencer chars like ':' or '/'
                clean_name = record.id.split("_")[0] + "_consensus"
                print(clean_name)
                file_path = os.path.join(output.out_dir, f"{clean_name}.fastq")
                with open(file_path, "w") as out_f:
                    SeqIO.write(record, out_f, "fastq")


# Helper function that determines which files are produced by the split_fastq checkpoint
# -----------------------------------------------------
def get_split_files(wildcards):
    # 1. Ensure the checkpoint has finished
    checkpoint_output = checkpoints.split_fastq.get(**wildcards).output[0]
    # 2. List the files that were actually created
    # This glob_wildcards looks inside the directory produced by the checkpoint
    filenames = glob_wildcards(os.path.join(checkpoint_output, "{sample}_consensus.fastq")).sample
    print(filenames)
    # 3. Return the full paths to the next rule
    return expand(os.path.join(checkpoint_output, "{sample}_consensus.fastq"), sample=filenames)


# Rule to create a summary CSV with sample, sequence, and quality score columns from all consensus FASTA files
# -----------------------------------------------------
rule consensus_summary_csv:
    input:
        consensuses = get_split_files
    output:
        report("consensus/consensus_summary.csv", category = "consensus")
    run:
        import csv
        from Bio import SeqIO
        with open(output[0], "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["sample", "sequence", "qscore"])
            print(input)
            for fastq_path in input.consensuses:
                sample = os.path.basename(fastq_path).replace("_consensus.fastq", "")
                for record in SeqIO.parse(fastq_path, "fastq"):
                    qscore = sum(record.letter_annotations["phred_quality"]) / len(record)
                    writer.writerow([sample, str(record.seq), qscore])

