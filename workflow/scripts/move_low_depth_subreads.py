"""
Identify samples below depth threshold and create placeholder files for downstream moving.
"""
import pandas as pd
import os
from common import log_msg

with open(snakemake.log.logf, "w") as logf:
    log_msg(logf, "Starting low depth subread identification")
    log_msg(logf, f"Depth threshold: {snakemake.params.threshold}")
    
    log_msg(logf, "Creating output directories")
    os.makedirs(snakemake.output.out, exist_ok=True)
    os.makedirs("demux/.low_depth", exist_ok=True)
    
    log_msg(logf, f"Reading demux statistics from: {snakemake.input.csv}")
    df = pd.read_csv(snakemake.input.csv)
    log_msg(logf, f"Total samples in statistics: {len(df)}")
    
    # Identify files below threshold
    log_msg(logf, f"Identifying samples with reads < {snakemake.params.threshold}")
    low_read_files = df[df['n_reads'] < snakemake.params.threshold]['sample'].tolist()
    log_msg(logf, f"Found {len(low_read_files)} samples below threshold")
    
    # Create empty "flag" files to tell Snakemake what to move
    log_msg(logf, "Creating placeholder files for low depth samples")
    for f in low_read_files:
        flag_path = os.path.join(snakemake.output.out, f"{f}.fastq.gz")
        with open(flag_path, 'w') as out:
            out.write("")
        log_msg(logf, f"  - Created placeholder: {flag_path}")
    
    log_msg(logf, f"Low depth identification complete. Marked {len(low_read_files)} samples for moving")