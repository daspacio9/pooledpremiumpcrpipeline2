# \HEADER\-------------------------------------------------------------------------
#
#  CONTENTS      : Snakemake nanopore data pipeline
#
#  DESCRIPTION   : remove large files and move analysis files to timestamped folder
#                   so that analysis may be repeated
#
#  RESTRICTIONS  : none
#
#  REQUIRES      : none
#
# ---------------------------------------------------------------------------------

# imports
# -----------------------------------------------------
import os, sys, glob
from datetime import datetime


# local rules
localrules:
    make_pppp_output_dir,
    alignment_clean,
    demux_clean,
    consensus_clean,
    logs_clean,
    report_clean,
    clean,


config["timestamp"] = datetime.now().strftime("%Y_%m_%d-%H%M%S")


# clean up sequences folder: remove combined .fastq.gz files (but not basecalled batches), move UMI stats files
# -----------------------------------------------------
rule make_pppp_output_dir:
    output:
        touch(".make_pppp_output_dir.done"),  # changed: consistent touch filename
    log:
        "cleanlogs/make_pppp_output_dir.log",
    conda:
        "../envs/clean.yaml"
    params:
        timestampDir=lambda wildcards: config["timestamp"] + "-ppppOutputs",
    shell:
        """
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Creating output directory" >> {log}
        mkdir -p {params.timestampDir} >> {log} 2>&1
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Output directory created" >> {log}
        """


# clean up compute batches alignment
# -----------------------------------------------------
rule alignment_clean:
    input:
        rules.make_pppp_output_dir.output,
    output:
        touch(".alignment_clean.done"),
    log:
        "cleanlogs/alignment_clean.log",
    conda:
        "../envs/clean.yaml"
    shell:
        """
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting alignment directory cleanup" >> {log}
        if [ -d aln ]; then
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Removing aln directory" >> {log}
            rm -rv aln >> {log} 2>&1
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Alignment directory removed" >> {log}
        fi
        """


# clean up compute batches demux
# -----------------------------------------------------
rule demux_clean:
    input:
        rules.make_pppp_output_dir.output,
    output:
        touch(".demux_clean.done"),
    log:
        "cleanlogs/demux_clean.log",
    conda:
        "../envs/clean.yaml"
    params:
        timestampDir=lambda wildcards: config["timestamp"] + "-ppppOutputs",
    shell:
        """
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting demux directory cleanup" >> {log}
        if [ -d demux ]; then
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Removing demux directory" >> {log}
            rm -rv demux >> {log} 2>&1
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Demux directory removed" >> {log}
        fi
        """


# clean up consensus outputs
# -----------------------------------------------------
rule consensus_clean:
    input:
        rules.make_pppp_output_dir.output,
    output:
        touch(".consensus_clean.done"),
    log:
        "cleanlogs/consensus_clean.log",
    conda:
        "../envs/clean.yaml"
    params:
        timestampDir=lambda wildcards: config["timestamp"] + "-ppppOutputs",
    shell:
        """
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting consensus directory cleanup" >> {log}
        if [ -d consensus ]; then
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Creating output directories for consensus" >> {log}
            mkdir -p {params.timestampDir}/consensus >> {log} 2>&1
            mkdir -p {params.timestampDir}/consensus_split >> {log} 2>&1
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Moving consensus directories" >> {log}
            mv consensus {params.timestampDir}/consensus/ >> {log} 2>&1
            mv consensus_split {params.timestampDir}/consensus_split/ >> {log} 2>&1
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Consensus directories moved" >> {log}
        fi
        """


# clean up logs  rm -r {params.timestampDir}/consensus_split/.*.done
# -----------------------------------------------------
rule logs_clean:
    input:
        rules.make_pppp_output_dir.output,
    output:
        touch(".logs_clean.done"),
    log:
        "cleanlogs/logs_clean.log",
    conda:
        "../envs/clean.yaml"
    params:
        timestampDir=lambda wildcards: config["timestamp"] + "-ppppOutputs",
    shell:
        """
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting logs directory cleanup" >> {log}
        if [ -d logs ]; then
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Moving logs directory" >> {log}
            mkdir -p {params.timestampDir}/logs
            mv logs {params.timestampDir}/ >> {log} 2>&1
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Logs directory moved" >> {log}
        fi
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Copying config files" >> {log}
        mkdir -p {params.timestampDir}/config
        cp config/*.yaml {params.timestampDir}/config/ >> {log} 2>&1
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Copying reference files" >> {log}
        cp -r ref {params.timestampDir}/ref/ >> {log} 2>&1
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Logs cleanup complete" >> {log}
        """


# clean up reports
# -----------------------------------------------------
rule report_clean:
    input:
        rules.make_pppp_output_dir.output,
    output:
        touch(".report_clean.done"),
    log:
        "cleanlogs/report_clean.log",
    conda:
        "../envs/clean.yaml"
    params:
        timestampDir=lambda wildcards: config["timestamp"] + "-ppppOutputs",
    shell:
        """
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting report directory cleanup" >> {log}
        if [ -d report ]; then
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Moving report directory" >> {log}
            mv report {params.timestampDir}/report/ >> {log} 2>&1
        fi
        if [ -f demux_stats.csv ]; then
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Moving demux_stats.csv" >> {log}
            mv demux_stats.csv {params.timestampDir}/ >> {log} 2>&1
        fi
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Report cleanup complete" >> {log}
        """


rule ab1_clean:
    input:
        rules.make_pppp_output_dir.output,
    output:
        touch(".ab1_clean.done"),
    log:
        "cleanlogs/ab1_clean.log",
    conda:
        "../envs/clean.yaml"
    params:
        timestampDir=lambda wildcards: config["timestamp"] + "-ppppOutputs",
    shell:
        """
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting ab1 directory cleanup" >> {log}
        if [ -d ab1 ]; then
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Moving ab1 directory" >> {log}
            mv ab1/ {params.timestampDir}/ >> {log} 2>&1
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] AB1 directory moved" >> {log}
        fi
        """


# clean up intermediate files
# -----------------------------------------------------
rule clean:
    input:
        rules.make_pppp_output_dir.output,
        rules.alignment_clean.output,
        rules.demux_clean.output,
        rules.consensus_clean.output,
        rules.logs_clean.output,
        rules.report_clean.output,
        rules.ab1_clean.output,
    log:
        "cleanlogs/clean.log",
    conda:
        "../envs/clean.yaml"
    params:
        timestampDir=lambda wildcards: config["timestamp"] + "-ppppOutputs",
    shell:
        """
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting final cleanup" >> {log}
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Removing intermediate files" >> {log}
        rm {input} >> {log} 2>&1
        if [ -d cleanlogs ]; then
            mv cleanlogs {params.timestampDir}
        fi
        zip -r -m {params.timestampDir}.zip {params.timestampDir}
        """
