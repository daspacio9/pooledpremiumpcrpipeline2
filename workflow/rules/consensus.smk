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
        inf = "demux/{sample}.fastq.gz"
    output:
        outf = "demux/temp/{sample}_subreads_batch.fastq.gz"
    log:
        logf = "logs/consensus/{sample}_fastq_to_fastq_subreads.log"
    run:
        from Bio import SeqIO
        import gzip
        records = []
        with gzip.open(input.inf, "rt") as fin:
            for i, record in enumerate(SeqIO.parse(fin, "fastq")):
                rec_id = f"{wildcards.sample}_subread{i+1}"
                record.id = rec_id
                record.description = ""
                records.append(record)
        with open(output.outf, "w") as fout:
            SeqIO.write(records, fout, "fastq")
        with open(log.logf, "w") as logf:
            logf.write("Demux stats prepared at {output.outf}\n")


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
    params:
        min_depth = config["min_depth"],
        medaka_model = config["medaka_model"]
    threads: config["medaka_spoa_threads"]
    shell:
        r"""
        rm -rf {output.outDir} >> {log} 2>&1
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
        """

split_dir = "consensus_split"
# Rule to split the parallel process fastq back into individual files
# -----------------------------------------------------
checkpoint split_fastq:
    input:
        batchconsensus = "consensus/bulk_consensus/consensus.fastq"
    output:
        out_dir = directory(split_dir)
    log: 
        logf = "logs/consensus/split_fastq.log"
    run:
        import os
        from pathlib import Path
        from Bio import SeqIO # Requires biopython
        Path(output.out_dir).mkdir(parents=True, exist_ok=True)
        with open(input.batchconsensus, "r") as handle:
            for record in SeqIO.parse(handle, "fastq"):
                # Clean up the ID: Benchling doesn't like special chars in names
                # We strip common sequencer chars like ':' or '/'
                clean_name = record.id.split("_")[0] + "_consensus"
                file_path = os.path.join(output.out_dir, f"{clean_name}.fastq")

                
                with open(log.logf, "w") as logf:
                    logf.write(f"Writing consensus for sample {clean_name}\n")
                    logf.write(f"to file {file_path}\n")
                with open(file_path, "w") as out_f:
                    SeqIO.write(record, out_f, "fastq")

# Rule to create a summary CSV with sample, sequence, and quality score columns from all consensus FASTA files
# -----------------------------------------------------
rule consensus_summary_csv:
    input:
        # consensuses = get_split_files
        consensuses = lambda wildcards: get_splits(wildcards, prefix="", suffix="_consensus.fastq")
    output:
        outf = report("consensus/consensus_summary.csv", category = "consensus")
    log: 
        logf = "logs/consensus/consensus_summary_csv.log"
    run:
        import csv
        from Bio import SeqIO
        with open(output.outf, "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["sample", "sequence", "qscore"])
            with open(log.logf, "w") as logf:
                logf.write("Generating consensus summary CSV for {}\n".format(input.consensuses))
            for fastq_path in input.consensuses:
                sample = os.path.basename(fastq_path).replace("_consensus.fastq", "")
                for record in SeqIO.parse(fastq_path, "fastq"):
                    qscore = sum(record.letter_annotations["phred_quality"]) / len(record)
                    writer.writerow([sample, str(record.seq), qscore])

