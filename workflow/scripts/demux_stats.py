"""
Count reads in each demuxed fastq.gz file and write statistics to demux_stats.csv.
"""
import gzip
import os

with open(snakemake.output[0], "w") as outf:
    outf.write("sample,n_reads\n")
    
    with open(snakemake.log.logf, "w") as logf:
        logf.write("Demux stats log\n")
        
        for samplefile in snakemake.input.samplefiles + [snakemake.input.unkfiles]:
            sample = os.path.basename(samplefile).replace(".fastq.gz", "")
            with gzip.open(samplefile, "rt") as f:
                cntr1 = 0
                for line in f:
                    cntr1 += 1
            n_reads = cntr1 // 4
            outf.write(f"{sample},{n_reads}\n")
            logf.write(f"Sample {sample}: {n_reads} reads")