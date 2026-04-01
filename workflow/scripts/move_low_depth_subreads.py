"""
Identify samples below depth threshold and create placeholder files for downstream moving.
"""
import pandas as pd
import os

os.makedirs(snakemake.output.out, exist_ok=True)
os.makedirs("demux/.low_depth", exist_ok=True)

df = pd.read_csv(snakemake.input.csv)

# Identify files below threshold
low_read_files = df[df['n_reads'] < snakemake.params.threshold]['sample'].tolist()

with open(snakemake.log.logf, "w") as logf:
    logf.write("Demux stats log\n")

# Create empty "flag" files to tell Snakemake what to move
for f in low_read_files:
    with open(os.path.join(snakemake.output.out, f"{f}.fastq.gz"), 'w') as out:
        out.write("")