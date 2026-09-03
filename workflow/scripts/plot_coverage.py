"""
Generate coverage and mismatch frequency plots from pypileup data.
Creates two PDF files: coverage plot and substitution frequency plot.
"""
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import log_msg, log_warning, print_colored, RESET, PURPLE
import sys

with open(snakemake.log.logf, "w") as logf:
    log_msg(logf, f"Starting coverage plot generation for sample: {snakemake.wildcards.sample}")
    
    log_msg(logf, f"Reading pypileup data from: {snakemake.input[0]}")
    df_ref = pd.read_csv(snakemake.input[0], sep="\t")
    log_msg(logf, f"Loaded {len(df_ref)} coverage data points")
    log_msg(logf, f"Coverage statistics - Min: {df_ref['coverage'].min()}, Max: {df_ref['coverage'].max()}, Mean: {df_ref['coverage'].mean():.2f}")
    
    log_msg(logf, "Calculating mismatch frequencies")
    df_ref['freq_mismatch'] = (df_ref['coverage'] - df_ref['ref_match']) / df_ref['coverage']
    df_ref['freq_mismatch'] = df_ref['freq_mismatch'].fillna(1)
    log_msg(logf, f"Mismatch statistics - Min: {df_ref['freq_mismatch'].min():.4f}, Max: {df_ref['freq_mismatch'].max():.4f}, Mean: {df_ref['freq_mismatch'].mean():.4f}")
    
    # Plot coverage
    log_msg(logf, f"Generating coverage plot to {snakemake.output[0]}")
    plt.figure(figsize=(10, 4))
    plt.plot(df_ref["pos"], df_ref["coverage"], color='blue', linewidth=0.2)
    plt.plot(df_ref["pos"], df_ref["ref_match"], color='green', linewidth=0.4, alpha=0.5)
    plt.title(f"Coverage Plot for {snakemake.wildcards.sample}")
    plt.xlabel("Position")
    plt.ylabel("Depth")
    plt.yscale("linear")
    plt.ylim(0, round(df_ref["coverage"].max()*1.1))
    plt.legend(["Coverage", "Ref Match"])
    plt.grid(False)
    plt.savefig(snakemake.output[0])
    plt.close()
    log_msg(logf, "Coverage plot generated successfully")

    # Plot frequency of mismatches
    log_msg(logf, f"Generating mismatch frequency plot to {snakemake.output[1]}")
    plt.figure(figsize=(10, 4))
    plt.plot(df_ref["pos"], df_ref["freq_mismatch"], color='red', linewidth=0.5)
    plt.fill_between(df_ref["pos"], df_ref["freq_mismatch"], color='red', alpha=0.3)
    plt.title(f"Substitution Plot for {snakemake.wildcards.sample}")
    plt.xlabel("Position")
    plt.ylabel("freq")
    plt.yscale("linear")
    plt.ylim(0, (df_ref["freq_mismatch"]).max()*1.1)
    plt.legend(["freq_mismatch"])
    plt.grid(False)
    plt.savefig(snakemake.output[1])
    plt.close()
    log_msg(logf, "Mismatch frequency plot generated successfully")

        
    # Identify quality issues
    low_coverage_threshold = df_ref['coverage'].max() * 0.2  # or use fixed value like 50
    high_mismatch_threshold = 0.30

    low_coverage_positions = df_ref[df_ref['coverage'] < low_coverage_threshold]
    high_mismatch_positions = df_ref[df_ref['freq_mismatch'] > high_mismatch_threshold]

    num_low_coverage = len(low_coverage_positions)
    num_high_mismatch = len(high_mismatch_positions)

    # Extract barcode group from sample name (adjust pattern based on your naming)
    sample_name = snakemake.wildcards.sample
    # Example: if sample is "group-C1" extract "C1" or map from config

    # Create warning messages
    if num_low_coverage > 0 or num_high_mismatch > 0:
        warning_msg = f"⚠ QC WARNINGS for {sample_name}:"
        log_warning(f"{PURPLE}⚠ QC WARNINGS for {sample_name}:{RESET}")
        print_colored(warning_msg)
        log_msg(logf, warning_msg)
        
        if num_low_coverage > 0:
            msg = f"  • Low coverage (<20%): {num_low_coverage} bases"
            print_colored(msg)
            log_msg(logf, msg)
            
        if num_high_mismatch > 0:
            msg = f"  • High mismatch (>30%): {num_high_mismatch} bases"
            print_colored(msg)
            log_msg(logf, msg)

    with open(snakemake.output.plot_flag, "w") as _f:
        _f.write("Plotting completed.\n")
    
    log_msg(logf, "All plotting steps completed successfully")