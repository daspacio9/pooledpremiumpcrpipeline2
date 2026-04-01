"""
Classify reads in fastq to subreads following the required formatting by medaka smolecule package.
Renames reads with sample name and subread index.
"""
from Bio import SeqIO
import gzip

records = []
with gzip.open(snakemake.input.inf, "rt") as fin:
    for i, record in enumerate(SeqIO.parse(fin, "fastq")):
        rec_id = f"{snakemake.wildcards.sample}_subread{i+1}"
        record.id = rec_id
        record.description = ""
        records.append(record)

with open(snakemake.output.outf, "w") as fout:
    SeqIO.write(records, fout, "fastq")

with open(snakemake.log.logf, "w") as logf:
    logf.write(f"Demux stats prepared at {snakemake.output.outf}\n")