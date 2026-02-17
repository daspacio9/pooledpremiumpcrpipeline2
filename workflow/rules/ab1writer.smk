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
    params:
        chunk_size = 5000
    run:
        import pandas as pd
        import os
        
        df = pd.read_csv(input[0], sep='\t')
        chunk_size = params.chunk_size
        chunks_dir = "report/pileup_chunks"
        os.makedirs(chunks_dir, exist_ok=True)
        
        num_chunks = 0
        for i, chunk_start in enumerate(range(0, len(df), chunk_size)):
            chunk_end = min(chunk_start + chunk_size, len(df))
            df_chunk = df.iloc[chunk_start:chunk_end]
            chunk_file = f"{chunks_dir}/{wildcards.sample}_pypileup_chunk_{i}.tsv"
            df_chunk.to_csv(chunk_file, sep='\t', index=False)
            num_chunks = i + 1
        
        with open(log[0], 'w') as f:
            f.write(f"Chunked {wildcards.sample} pypileup into {num_chunks} chunks\n")


rule write_ab1:
    input:
        chunk_flag = "report/.chunks_{sample}.done",
        chunk_file = "report/pileup_chunks/{sample}_pypileup_chunk_{chunk}.tsv"
    output:
        flag = touch("ab1/.ab1_{sample}_{chunk}.done"),
        ab1_file = "ab1/{sample}_{chunk}.ab1"
    log:
        "logs/ab1writer/{sample}_{chunk}_write_ab1.log"
    run:
        from common import write_ab1, generate_trace, parse_pypileup
        
        df_norm, cons_seq = parse_pypileup(input.chunk_file)
        base_counts = generate_trace(df_norm, cons_seq)
        write_ab1(output.ab1_file, cons_seq, base_counts['G'], base_counts['A'], base_counts['T'], base_counts['C'])

rule ab1_done:
    input:
        ab1_files = lambda wildcards: get_all_ab1_files(wildcards)
    output:
        touch(".ab1_done.done")
    log:
        logf = "logs/aln/ab1_done.log"
    run:
        with open(output[0], "w") as f:
            f.write("AB1 writing completed.\n")