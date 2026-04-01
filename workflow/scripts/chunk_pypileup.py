"""
Chunk the pypileup file for the ab1 writer at 5 kb intervals.
AB1 format has a max length of ~8 kb, but 5kb chunks are used as a standard size.
"""
import pandas as pd
import os
from common import log_msg

with open(snakemake.log[0], 'w') as logf:
    log_msg(logf, f"Starting pypileup chunking for sample: {snakemake.wildcards.sample}")
    log_msg(logf, f"Input file: {snakemake.input[0]}")
    
    log_msg(logf, "Reading pypileup file")
    df = pd.read_csv(snakemake.input[0], sep='\t')
    log_msg(logf, f"Total rows in pypileup: {len(df)}")
    
    chunk_size = snakemake.params.chunk_size
    log_msg(logf, f"Chunk size: {chunk_size} rows per chunk")
    
    chunks_dir = "report/pileup_chunks"
    log_msg(logf, f"Creating chunks directory: {chunks_dir}")
    os.makedirs(chunks_dir, exist_ok=True)
    
    log_msg(logf, "Creating chunks")
    num_chunks = 0
    for i, chunk_start in enumerate(range(0, len(df), chunk_size)):
        chunk_end = min(chunk_start + chunk_size, len(df))
        df_chunk = df.iloc[chunk_start:chunk_end]
        chunk_file = f"{chunks_dir}/{snakemake.wildcards.sample}_pypileup_chunk_{i}.tsv"
        df_chunk.to_csv(chunk_file, sep='\t', index=False)
        log_msg(logf, f"  Chunk {i}: rows {chunk_start}-{chunk_end} ({len(df_chunk)} rows) -> {chunk_file}")
        num_chunks = i + 1
    
    log_msg(logf, f"Chunking complete. Created {num_chunks} chunks for {snakemake.wildcards.sample}")