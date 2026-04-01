"""
Generate AB1 file from pypileup chunk data.
Parses pypileup data, generates base trace, and writes to AB1 format.
"""
from common import write_ab1, generate_trace, parse_pypileup, log_msg
with open(snakemake.log.logf, "w") as logf:
    log_msg(logf, f"Starting AB1 file generation for {snakemake.wildcards.sample} chunk {snakemake.wildcards.chunk}")
    log_msg(logf, f"Input pypileup chunk: {snakemake.input.chunk_file}")

    log_msg(logf, "Parsing pypileup data")
    df_norm, cons_seq = parse_pypileup(snakemake.input.chunk_file)
    log_msg(logf, f"Parsed consensus sequence length: {len(cons_seq)} bp")

    log_msg(logf, "Generating base trace")
    base_counts = generate_trace(df_norm, cons_seq)
    log_msg(logf, f"Generated trace for bases: G, A, T, C")

    log_msg(logf, f"Writing AB1 file to: {snakemake.output.ab1_file}")
    write_ab1(snakemake.output.ab1_file, cons_seq, base_counts['G'], base_counts['A'], base_counts['T'], base_counts['C'])

    log_msg(logf, "AB1 file written successfully")
