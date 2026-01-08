# \HEADER\-------------------------------------------------------------------------
#
#  CONTENTS      : Snakemake nanopore data pipeline
#
#  DESCRIPTION   : demultiplexing of nanopore reads by end barcodes using cutadapt package
#
#  RESTRICTIONS  : none
#
#  REQUIRES      : none
#
# ---------------------------------------------------------------------------------


#Demultiplexing rule using cutadapt with linked adapters
# -----------------------------------------------------
rule cutadapt_demux_linked:
    input:
        "ref/barcode-pairs.txt",
        f"sequences/{config['input_fastq']}"
    output:
        # all per-sample fastqs + an explicit file for unmatched
        demux=expand("demux/{s}.fastq.gz", s=adapter_names("ref/barcode-pairs.txt")),
        unknown="demux/unknown.fastq.gz"
    log:
        "logs/cutadapt/demux.log"
    shell:
        r""" 
        cutadapt \
          -g file:{input[0]} \
          --revcomp \
          -e {config[cutadapt_error_rate]} \
          -O {config[cutadapt_min_overlap]} \
          --action=none \
          -o demux/{{name}}.fastq.gz \
          --untrimmed-output {output.unknown} \
          {input[1]} \
        > {log} 2>&1
        """

# Checkpoint to count reads in each demuxed fastq.gz file and write to demux_stats.csv
# -----------------------------------------------------
checkpoint demux_stats:
    input:
        samplefiles = expand("demux/{s}.fastq.gz", s=adapter_names("ref/barcode-pairs.txt")),
        unkfiles = "demux/unknown.fastq.gz"
    output:
        report("demux_stats.csv", category = "demux", labels={"type": "stats", "status": "Unfiltered"})
    run:
        import gzip

        with open(output[0], "w") as outf:
            outf.write("sample,n_reads\n")
            print(input.samplefiles)
            for samplefile in input.samplefiles + [input.unkfiles]:
                sample = os.path.basename(samplefile).replace(".fastq.gz", "")
                with gzip.open(samplefile, "rt") as f:
                    cntr1 = 0
                    for line in f:
                        cntr1 += 1
                n_reads = cntr1 // 4
                outf.write(f"{sample},{n_reads}\n")
                print(f"Sample {sample}: {n_reads} reads")


# Checkpoint to count reads in each demuxed fastq.gz file and write to demux_stats.csv
# -----------------------------------------------------
checkpoint move_low_depth_subreads:
    input:
        csv = "demux_stats.csv",
    output:
        (directory("demux/filtered_list"))
    params:
        threshold = config["min_depth"]
    run:
        import pandas as pd
        os.makedirs(output[0], exist_ok=True)
        df = pd.read_csv(input.csv)
        # Identify files below threshold
        low_read_files = df[df['n_reads'] < params.threshold]['sample'].tolist()
        print(low_read_files)
        # Create empty "flag" files to tell Snakemake what to move
        for f in low_read_files:
            with open(os.path.join(output[0], f + ".fastq.gz"), 'w') as out:
                out.write("")


# Function that depends on low depth checkpoint to determine which files to move
# -----------------------------------------------------
def get_moved_files(wildcards):
    # This forces Snakemake to wait until the checkpoint is done
    checkpoint_output = checkpoints.move_low_depth_subreads.get(**wildcards).output[0]
    # Return the expected paths in the final destination
    return [os.path.join("demux/.low_depth", f) for f in os.listdir(checkpoint_output) if not f.startswith('.')]

# File system I/O
# -----------------------------------------------------
rule move_file:
    input:
        known = "demux/{filename}.fastq.gz"
    output:
        known = "demux/.low_depth/{filename}.fastq.gz"
    shell:
        "mv {input.known} {output.known}"

# Rule to finalize demultiplexing by removing placeholder files
# -----------------------------------------------------
rule finalize_demux:
    input:
        get_moved_files
    output:
        ".demux_finished.flag"
    shell:
        """
        mv demux/unknown.fastq.gz demux/.low_depth/unknown.fastq.gz
        touch {output}
        rm -r demux/filtered_list
        """