# \HEADER\-------------------------------------------------------------------------
#
#  CONTENTS      : Snakemake nanopore data pipeline
#
#  DESCRIPTION   : Workflow rules for writing the ab1 files from the pypileup chunks. Ab1 writer from Vincent Hu.
#  RESTRICTIONS  : none
#
#  REQUIRES      : none
#
# ---------------------------------------------------------------------------------

### Chunking the pypileup file for the ab1 writer at 5 kb. ab1 format has a max length of ~8 kb according to plasmidsaurus, but 5kb is a nice round number.
# -----------------------------------------------------
checkpoint chunk_pypileup:
    input:
        "report/{sample}_pypileup.tsv"
    output:
        chunk_flag = touch("report/.chunks_{sample}.done")
    log:
        "logs/aln/{sample}_chunk_pypileup.log"
    run:
        import pandas as pd
        import os
        
        # Read the pypileup file
        df = pd.read_csv(input[0], sep='\t')
        
        # Create chunks directory if it doesn't exist
        chunks_dir = os.path.join(os.path.dirname(input[0]), "chunks")
        os.makedirs(chunks_dir, exist_ok=True)
        
        # Split into chunks of 5000 rows
        chunk_size = 5000
        for i, chunk_start in enumerate(range(0, len(df), chunk_size)):
            chunk_end = min(chunk_start + chunk_size, len(df))
            df_chunk = df.iloc[chunk_start:chunk_end]
            
            # Output chunk file
            chunk_file = f"{chunks_dir}/{wildcards.sample}_pypileup_chunk_{i}.tsv"
            df_chunk.to_csv(chunk_file, sep='\t', index=False)
        
        with open(log[0], 'w') as f:
            f.write(f"Chunked {wildcards.sample} pypileup into {i+1} chunks of {chunk_size} rows\n")


### Chunking the pypileup file for the ab1 writer
# -----------------------------------------------------
rule write_ab1:
    input:
        chunk_flag = "report/.chunks_{sample}.done",
        chunks = "report/pileup_chunks/{sample}_pypileup_chunk_{chunk}.tsv"
    output:
        flag = touch("ab1/.ab1_{sample}_{chunk}.done"), #flag file to indicate this chunk is done,
        ab1_file = "ab1/{sample}_{chunk}.ab1"
    log:
        "logs/ab1writer/{sample}_{chunk}_write_ab1.log"
    run:
        import pandas as pd
        from common import write_ab1
        from common import generate_trace
        
        # Read each chunk
        for chunk_file in input.chunks:
            if os.path.exists(chunk_file):
                base_counts, seq = generate_trace(chunk_file)
                write_ab1(output.ab1_file, seq, base_counts['G'], base_counts['A'], base_counts['T'], base_counts['C'])
                with open(log[0], 'w') as f:
                    f.write(f"Wrote AB1 file for {wildcards.sample}_{wildcards.chunk} with {len(seq)} bases\n")


rule ab1_done:
    input:
        lambda wildcards: get_splits(wildcards, prefix=".ab1_", suffix=".done")
    output:
        touch(".ab1.done")
    log:
        logf = "logs/aln/ab1_done.log"
    run:
        with open(output[0], "w") as f:
            f.write("AB1 writing completed.\n")