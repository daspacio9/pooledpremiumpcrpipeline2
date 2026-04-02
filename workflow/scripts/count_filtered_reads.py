"""
Count reads filtered out by length threshold and report statistics.
Compares input reads vs output reads to calculate filtered count.
"""
import gzip
from common import log_msg

def count_lines_in_file(filepath):
    """Count lines in a file, handling both gzipped and uncompressed formats."""
    try:
        # Try to open as gzipped first
        with gzip.open(filepath, "rt") as f:
            return sum(1 for _ in f)
    except (OSError, gzip.BadGzipFile):
        # Fall back to regular file if not gzipped
        with open(filepath, "r") as f:
            return sum(1 for _ in f)

with open(snakemake.log.logf, "w") as logf:
    log_msg(logf, "Starting filtered read count analysis")
    log_msg(logf, f"Filter threshold: {snakemake.params.min_length} bp")
    
    # Count input reads
    log_msg(logf, f"Counting reads in input: {snakemake.input.source}")
    input_lines = count_lines_in_file(snakemake.input.source)
    input_reads = input_lines // 4
    log_msg(logf, f"Input reads: {input_reads}")
    
    # Count output (filtered) reads
    log_msg(logf, f"Counting reads in output: {snakemake.input.filtered}")
    output_lines = count_lines_in_file(snakemake.input.filtered)
    output_reads = output_lines // 4
    log_msg(logf, f"Output reads (after filtering): {output_reads}")
    
    # Calculate filtered count
    filtered_count = input_reads - output_reads
    filter_rate = (filtered_count / input_reads * 100) if input_reads > 0 else 0
    
    log_msg(logf, f"Reads filtered out: {filtered_count} ({filter_rate:.2f}%)")
    log_msg(logf, f"Retention rate: {(100 - filter_rate):.2f}%")
    
    # Write statistics file
    with open(snakemake.output.stats, "w") as f:
        f.write("metric,value\n")
        f.write(f"input_reads,{input_reads}\n")
        f.write(f"output_reads,{output_reads}\n")
        f.write(f"filtered_reads,{filtered_count}\n")
        f.write(f"filter_threshold_bp,{snakemake.params.min_length}\n")
    
    log_msg(logf, f"Filter statistics written to: {snakemake.output.stats}")
