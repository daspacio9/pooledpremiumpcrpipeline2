"""
Create a summary CSV with sample, sequence, and quality score columns from all consensus FASTA files.
"""
import csv
import os
from Bio import SeqIO
from common import log_msg

with open(snakemake.log.logf, "w") as logf:
    log_msg(logf, "Starting consensus summary CSV generation")
    log_msg(logf, f"Processing {len(snakemake.input.consensuses)} consensus files")
    
    with open(snakemake.output.outf, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["sample", "sequence", "qscore"])
        log_msg(logf, f"CSV header written to {snakemake.output.outf}")
        
        total_records = 0
        for idx, fastq_path in enumerate(snakemake.input.consensuses, 1):
            sample = os.path.basename(fastq_path).replace("_consensus.fastq", "")
            log_msg(logf, f"Processing file {idx}/{len(snakemake.input.consensuses)}: {sample}")
            
            record_count = 0
            for record in SeqIO.parse(fastq_path, "fastq"):
                qscore = sum(record.letter_annotations["phred_quality"]) / len(record)
                writer.writerow([sample, str(record.seq), qscore])
                record_count += 1
            
            total_records += record_count
            log_msg(logf, f"Finished {sample}: wrote {record_count} records")
        
        log_msg(logf, f"CSV generation complete. Total records written: {total_records}")