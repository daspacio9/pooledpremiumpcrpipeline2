"""
Split the parallel process fastq back into individual files by sample.
Extracts sample name from record ID and writes each sample's consensus to separate file.
"""
import os
from pathlib import Path
from Bio import SeqIO

Path(snakemake.output.out_dir).mkdir(parents=True, exist_ok=True)

with open(snakemake.input.batchconsensus, "r") as handle:
    for record in SeqIO.parse(handle, "fastq"):
        # Extract sample name from record ID (format: {sample}_subread{i})
        # Remove the '_subread' part to get back the original sample name
        parts = record.id.rsplit("_subread", 1)
        sample_name = parts[0]
        clean_name = sample_name + "_consensus"
        file_path = os.path.join(snakemake.output.out_dir, f"{clean_name}.fastq")

        with open(snakemake.log.logf, "w") as logf:
            logf.write(f"Writing consensus for sample {clean_name}\n")
            logf.write(f"to file {file_path}\n")
        
        with open(file_path, "w") as out_f:
            SeqIO.write(record, out_f, "fastq")