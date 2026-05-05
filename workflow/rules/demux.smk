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


# Demultiplexing rule using cutadapt with linked adapters
# -----------------------------------------------------
rule filter_reads_by_length:
    input:
        fastq=f"sequences/{config['input_fastq']}",
    output:
        filtered=f"sequences/filtered_{config['input_fastq']}",
    log:
        "logs/filter/filter_reads.log",
    conda:
        "../envs/demux.yaml"
    params:
        min_length=config["filter-size"],
    shell:
        r"""
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting read length filtering" >> {log}
        seqkit seq -m {params.min_length} {input.fastq} 2>> {log} | gzip -c > {output.filtered}
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Read filtering complete" >> {log}
        """


# Count statistics on filtered reads
# --------------------
rule count_filtered_reads:
    input:
        source=f"sequences/{config['input_fastq']}",
        filtered=f"sequences/filtered_{config['input_fastq']}",
    output:
        stats="logs/filter/filter_stats.csv",
    log:
        logf="logs/filter/count_filtered_reads.log",
    conda:
        "../envs/demux.yaml"
    params:
        min_length=config["filter-size"],
    script:
        "../scripts/count_filtered_reads.py"


# Demultiplexing rule using cutadapt with linked adapters
# -----------------------------------------------------
rule cutadapt_demux_linked:
    input:
        barcodes=use_debug(),
        seq=f"sequences/filtered_{config['input_fastq']}",
    output:
        # all per-sample fastqs + an explicit file for unmatched
        demux=expand("demux/{s}.fastq.gz", s=adapter_names(use_debug())),
        unknown="demux/unknown.fastq.gz",
    log:
        "logs/cutadapt/demux.log",
    conda:
        "../envs/demux.yaml"
    params:
        error_rate=config["cutadapt_error_rate"],
        min_overlap=config["cutadapt_min_overlap"],
    shell:
        r""" 
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting cutadapt demultiplexing" >> {log}
        cutadapt \
          -g file:{input.barcodes} \
          --revcomp \
          -e {params.error_rate} \
          -O {params.min_overlap} \
          --action=none \
          -o demux/{{name}}.fastq.gz \
          --untrimmed-output {output.unknown} \
          {input.seq} \
        >> {log} 2>&1
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Cutadapt demultiplexing complete" >> {log}
        """


# Checkpoint to count reads in each demuxed fastq.gz file and write to demux_stats.csv
# -----------------------------------------------------
checkpoint demux_stats:
    input:
        samplefiles=expand("demux/{s}.fastq.gz", s=adapter_names(use_debug())),
        unkfiles="demux/unknown.fastq.gz",
    output:
        report(
            "demux_stats.csv",
            category="demux",
            labels={"type": "stats", "status": "Unfiltered"},
        ),
    log:
        logf="logs/cutadapt/demux_stats.log",
    conda:
        "../envs/demux.yaml"
    script:
        "../scripts/demux_stats.py"


# Checkpoint to count reads in each demuxed fastq.gz file and write to demux_stats.csv
# Includes error handling if no barcode groups pass min_depth threshold
# Incorporates filter statistics for comprehensive error reporting
# -----------------------------------------------------
checkpoint move_low_depth_subreads:
    input:
        csv="demux_stats.csv",
        filter_stats="logs/filter/filter_stats.csv",
    output:
        out=(directory("demux/filtered_list")),
    log:
        logf="logs/cutadapt/move_low_depth_subreads.log",
    conda:
        "../envs/demux.yaml"
    params:
        threshold=config["min_depth"],
    script:
        "../scripts/move_low_depth_subreads.py"


# File system I/O
# -----------------------------------------------------
rule move_file:
    input:
        known="demux/{filename}.fastq.gz",
    output:
        known="demux/.low_depth/{filename}.fastq.gz",
    log:
        "logs/cutadapt/move_{filename}.log",
    conda:
        "../envs/demux.yaml"
    shell:
        "mv {input.known} {output.known} >> {log} 2>&1"


# Rule to finalize demultiplexing by removing placeholder files
# -----------------------------------------------------
rule finalize_demux:
    input:
        get_moved_files,
    output:
        ".demux_finished.done",
    log:
        "logs/cutadapt/finalize_demux.log",
    conda:
        "../envs/demux.yaml"
    shell:
        """
        mv demux/unknown.fastq.gz demux/.low_depth/unknown.fastq.gz >> {log} 2>&1
        touch {output} >> {log} 2>&1
        """
