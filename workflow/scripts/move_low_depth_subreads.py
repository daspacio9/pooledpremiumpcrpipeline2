"""
Identify samples below depth threshold and create placeholder files for downstream moving.
Includes error handling for when no barcode groups pass the min_depth threshold.
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
    
    # Identify files below threshold, excluding 'unknown' which is handled by finalize_demux
    log_msg(logf, f"Identifying samples with reads < {snakemake.params.threshold}")
    low_read_files = df[
        (df['n_reads'] < snakemake.params.threshold) & (df['sample'] != 'unknown')
    ]['sample'].tolist()
    log_msg(logf, f"Found {len(low_read_files)} samples below threshold (excluding 'unknown')")
    
    # Get passed samples (not including 'unknown')
    passed_samples = df[
        (df['n_reads'] >= snakemake.params.threshold) & (df['sample'] != 'unknown')
    ]['sample'].tolist()
    
    # Check if NO barcode groups pass the threshold (common critical error)
    if len(passed_samples) == 0:
        log_msg(logf, "ERROR: No barcode groups have reads >= min_depth threshold!")
        
        # Try to read filter statistics to report filtered sequences
        filtered_count = "unknown"
        filter_threshold = "unknown"
        try:
            if os.path.exists("logs/filter/filter_stats.csv"):
                stats_df = pd.read_csv("logs/filter/filter_stats.csv")
                stats_dict = dict(zip(stats_df['metric'], stats_df['value']))
                filtered_count = int(stats_dict.get('filtered_reads', 'unknown'))
                filter_threshold = int(stats_dict.get('filter_threshold_bp', 'unknown'))
                log_msg(logf, f"Filter statistics read: {filtered_count} reads filtered (below {filter_threshold} bp)")
        except Exception as e:
            log_msg(logf, f"Warning: Could not read filter statistics: {e}")
        
        # Build detailed error message
        total_barcode_groups = len(df) - (1 if 'unknown' in df['sample'].values else 0)
        error_msg = (
            f"\n{'='*80}\n"
            f"ERROR: PIPELINE FAILURE - All barcode groups below minimum depth threshold\n"
            f"{'='*80}\n"
            f"\n"
            f"DEPTH FILTER SUMMARY:\n"
            f"  - Min depth threshold: {snakemake.params.threshold} reads\n"
            f"  - Barcode groups with reads >= threshold: {len(passed_samples)}\n"
            f"  - Barcode groups with reads < threshold: {total_barcode_groups}\n"
            f"  - Total barcode groups evaluated: {total_barcode_groups}\n"
            f"\n"
            f"LENGTH FILTER SUMMARY:\n"
            f"  - Minimum read length filter: {filter_threshold} bp\n"
            f"  - Sequences filtered out: {filtered_count}\n"
            f"\n"
            f"AFFECTED BARCODE GROUPS:\n"
        )
        
        # List barcode groups with their read counts
        for sample in low_read_files:
            sample_reads = df[df['sample'] == sample]['n_reads'].values[0]
            error_msg += f"  - {sample}: {sample_reads} reads\n"
        
        error_msg += (
            f"\n"
            f"RECOMMENDATIONS:\n"
            f"  1. Check input data quality and quantity\n"
            f"  2. Verify barcode sequences and demultiplexing accuracy\n"
            f"  3. Consider lowering the 'min_depth' threshold in config if appropriate\n"
            f"  4. Review filter-size parameter - may be too restrictive\n"
            f"  5. Check for sequencing or library preparation issues\n"
            f"{'='*80}\n"
        )
        
        log_msg(logf, error_msg)
        raise RuntimeError(error_msg)
    
    log_msg(logf, f"Samples passing depth threshold: {len(passed_samples)}")
    for sample in passed_samples:
        sample_reads = df[df['sample'] == sample]['n_reads'].values[0]
        log_msg(logf, f"  - {sample}: {sample_reads} reads")
    
    # Create empty "flag" files to tell Snakemake what to move
    log_msg(logf, "Creating placeholder files for low depth samples")
    for f in low_read_files:
        flag_path = os.path.join(snakemake.output.out, f"{f}.fastq.gz")
        with open(flag_path, 'w') as out:
            out.write("")
        log_msg(logf, f"  - Created placeholder: {flag_path}")
    
    log_msg(logf, f"Low depth identification complete. Marked {len(low_read_files)} samples for moving")