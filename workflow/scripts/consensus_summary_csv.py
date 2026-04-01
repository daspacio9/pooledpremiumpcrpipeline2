"""
Create a summary CSV with sample, sequence, and quality score columns from all consensus FASTA files.
"""
import csv
import os
from Bio import SeqIO

with open(snakemake.output.outf, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["sample", "sequence", "qscore"])
    
    with open(snakemake.log.logf, "w") as logf:
        logf.write(f"Generating consensus summary CSV for {snakemake.input.consensuses}\n")
    
    for fastq_path in snakemake.input.consensuses:
        sample = os.path.basename(fastq_path).replace("_consensus.fastq", "")
        for record in SeqIO.parse(fastq_path, "fastq"):
            qscore = sum(record.letter_annotations["phred_quality"]) / len(record)
            writer.writerow([sample, str(record.seq), qscore])