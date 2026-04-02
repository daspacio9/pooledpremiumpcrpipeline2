"""
Classify reads in fastq to subreads following the required formatting by medaka smolecule package.
Renames reads with sample name and subread index.
"""
from Bio import SeqIO
import gzip
from common import log_msg

with open(snakemake.log.logf, "w") as logf:
    log_msg(logf, f"Starting fastq_to_fastq_subreads for sample: {snakemake.wildcards.sample}")
    log_msg(logf, f"Input file: {snakemake.input.inf}")
    
    log_msg(logf, "Opening input fastq file")
    records = []
    record_count = 0
    with gzip.open(snakemake.input.inf, "rt") as fin:
        for i, record in enumerate(SeqIO.parse(fin, "fastq")):
            rec_id = f"{snakemake.wildcards.sample}_subread{i+1}"
            record.id = rec_id
            record.description = ""
            records.append(record)
            record_count += 1
    
    log_msg(logf, f"Loaded {record_count} records from input file")
    log_msg(logf, f"Writing {record_count} renamed records to output file: {snakemake.output.outf}")
    
    with gzip.open(snakemake.output.outf, "wt") as fout:
        SeqIO.write(records, fout, "fastq")
    
    log_msg(logf, "Fastq to fastq subreads conversion completed successfully")