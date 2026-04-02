"""
Split the parallel process fastq back into individual files by sample.
Extracts sample name from record ID and writes each sample's consensus to separate file.
"""
import os
from pathlib import Path
from Bio import SeqIO
from collections import defaultdict
from common import log_msg

with open(snakemake.log.logf, "w") as logf:
    log_msg(logf, f"Starting fastq splitting from: {snakemake.input.batchconsensus}")
    log_msg(logf, f"Output directory: {snakemake.output.out_dir}")
    
    log_msg(logf, "Creating output directory")
    Path(snakemake.output.out_dir).mkdir(parents=True, exist_ok=True)
    
    log_msg(logf, "Reading batch consensus file and splitting by sample")
    sample_counts = defaultdict(int)
    
    with open(snakemake.input.batchconsensus, "r") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            # Extract sample name from record ID (format: {sample}_subread{i})
            parts = record.id.rsplit("_subread", 1)
            sample_name = parts[0]
            clean_name = sample_name + "_consensus"
            file_path = os.path.join(snakemake.output.out_dir, f"{clean_name}.fastq")
            
            log_msg(logf, f"Writing consensus for sample: {clean_name}")
            
            with open(file_path, "w") as out_f:
                SeqIO.write(record, out_f, "fastq")
            sample_counts[clean_name] += 1
    
    log_msg(logf, "Split summary:")
    for sample, count in sorted(sample_counts.items()):
        log_msg(logf, f"  - {sample}: {count} records")
    
    log_msg(logf, f"Fastq splitting completed. Created {len(sample_counts)} sample files")