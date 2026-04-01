"""
Count reads in each demuxed fastq.gz file and write statistics to demux_stats.csv.
"""
import gzip
import os
from common import log_msg

with open(snakemake.log.logf, "w") as logf:
    log_msg(logf, "Starting demux statistics generation")
    
    all_files = snakemake.input.samplefiles + [snakemake.input.unkfiles]
    log_msg(logf, f"Processing {len(all_files)} demux sample files")
    
    with open(snakemake.output[0], "w") as outf:
        outf.write("sample,n_reads\n")
        log_msg(logf, "CSV header written")
        
        total_reads = 0
        for idx, samplefile in enumerate(all_files, 1):
            sample = os.path.basename(samplefile).replace(".fastq.gz", "")
            log_msg(logf, f"Counting reads in file {idx}/{len(all_files)}: {sample}")
            
            with gzip.open(samplefile, "rt") as f:
                cntr1 = 0
                for line in f:
                    cntr1 += 1
            
            n_reads = cntr1 // 4
            total_reads += n_reads
            outf.write(f"{sample},{n_reads}\n")
            log_msg(logf, f"Sample {sample}: {n_reads} reads (file: {samplefile})")
        
        log_msg(logf, f"Demux statistics complete. Total reads across all samples: {total_reads}")