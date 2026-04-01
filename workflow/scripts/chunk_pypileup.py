"""
Chunk the pypileup file for the ab1 writer at 5 kb intervals.
AB1 format has a max length of ~8 kb, but 5kb chunks are used as a standard size.
"""
import pandas as pd
import os

df = pd.read_csv(snakemake.input[0], sep='\t')
chunk_size = snakemake.params.chunk_size
chunks_dir = "report/pileup_chunks"
os.makedirs(chunks_dir, exist_ok=True)

num_chunks = 0
for i, chunk_start in enumerate(range(0, len(df), chunk_size)):
    chunk_end = min(chunk_start + chunk_size, len(df))
    df_chunk = df.iloc[chunk_start:chunk_end]
    chunk_file = f"{chunks_dir}/{snakemake.wildcards.sample}_pypileup_chunk_{i}.tsv"
    df_chunk.to_csv(chunk_file, sep='\t', index=False)
    num_chunks = i + 1

with open(snakemake.log[0], 'w') as f:
    f.write(f"Chunked {snakemake.wildcards.sample} pypileup into {num_chunks} chunks\n")