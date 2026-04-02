"""
Parse mpileup file and generate pypileup TSV with base counts.
"""
import pandas as pd
from common import log_msg, parse_mpileup

with open(snakemake.log.logf, "w") as logf:
    log_msg(logf, f"Starting mpileup parsing for sample: {snakemake.wildcards.sample}")
    log_msg(logf, f"Input mpileup file: {snakemake.input[0]}")
    
    log_msg(logf, "Parsing mpileup file")
    df_pileup = parse_mpileup(snakemake.input[0])
    
    if df_pileup is None:
        error_msg = f"ERROR: Failed to parse mpileup file: {snakemake.input[0]}"
        log_msg(logf, error_msg)
        raise ValueError(error_msg)
    
    log_msg(logf, f"Parsed {len(df_pileup)} positions from mpileup")
    log_msg(logf, f"Writing pypileup TSV to {snakemake.output[0]}")
    
    df_pileup.to_csv(snakemake.output[0], sep="\t", index=False)
    
    log_msg(logf, "Mpileup parsing completed successfully")
