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
localrules: make_pppp_output_dir, alignment_clean, demux_clean, consensus_clean, logs_clean, report_clean, clean

config['timestamp'] = datetime.now().strftime("%Y_%m_%d-%H%M%S")
    
# clean up sequences folder: remove combined .fastq.gz files (but not basecalled batches), move UMI stats files
# -----------------------------------------------------
rule make_pppp_output_dir:
    output:
        touch('.make_pppp_output_dir.done')   # changed: consistent touch filename
    log:
        "cleanlogs/make_pppp_output_dir.log"
    params:
        timestampDir = lambda wildcards: config['timestamp'] + '-ppppOutputs'
    shell:
        """
        mkdir -p {params.timestampDir} >> {log} 2>&1
        """

# clean up compute batches alignment
# -----------------------------------------------------
rule alignment_clean:
    input:
        rules.make_pppp_output_dir.output
    output:
        touch('.alignment_clean.done')
    log:
        "cleanlogs/alignment_clean.log"
    shell:
        """
        if [ -d aln ]; then
            rm -rv aln >> {log} 2>&1
        fi
        """

# clean up compute batches demux
# -----------------------------------------------------
rule demux_clean:
    input:
        rules.make_pppp_output_dir.output
    output:
        touch('.demux_clean.done')
    log:
        "cleanlogs/demux_clean.log"
    params:
        timestampDir = lambda wildcards: config['timestamp'] + '-ppppOutputs'
    shell:
        """
        if [ -d demux ]; then
            rm -rv demux >> {log} 2>&1
        fi
        """

# clean up consensus outputs
# -----------------------------------------------------
rule consensus_clean:
    input:
        rules.make_pppp_output_dir.output
    output:
        touch('.consensus_clean.done')
    log:
        "cleanlogs/consensus_clean.log"
    params:
        timestampDir = lambda wildcards: config['timestamp'] + '-ppppOutputs'
    shell:
        """
        if [ -d consensus ]; then
            mkdir -p {params.timestampDir}/consensus >> {log} 2>&1
            mkdir -p {params.timestampDir}/consensus_split >> {log} 2>&1
            mv consensus {params.timestampDir}/consensus/ >> {log} 2>&1
            mv consensus_split {params.timestampDir}/consensus_split/ >> {log} 2>&1
        fi
        """

# clean up logs  rm -r {params.timestampDir}/consensus_split/.*.done
# -----------------------------------------------------
rule logs_clean:
    input:
        rules.make_pppp_output_dir.output
    output:
        touch('.logs_clean.done')
    params:
        timestampDir = lambda wildcards: config['timestamp'] + '-ppppOutputs'
    shell:
        """
        if [ -d logs ]; then
            mkdir -p {params.timestampDir}/logs
            mv logs {params.timestampDir}/
        fi
        mkdir -p {params.timestampDir}/config
        cp config/*.yaml {params.timestampDir}/config/
        cp -r ref {params.timestampDir}/ref/
        """

# clean up reports
# -----------------------------------------------------
rule report_clean:
    input:
        rules.make_pppp_output_dir.output
    output:
        touch('.report_clean.done')
    params:
        timestampDir = lambda wildcards: config['timestamp'] + '-ppppOutputs'
    log:
        "cleanlogs/report_clean.log"
    shell:
        """
        if [ -d report ]; then
            mv report {params.timestampDir}/report/ >> {log} 2>&1
        fi
        if [ -f demux_stats.csv ]; then
            mv demux_stats.csv {params.timestampDir}/ >> {log} 2>&1
        fi
        """
rule ab1_clean:
    input:
        rules.make_pppp_output_dir.output
    output:
        touch('.ab1_clean.done')
    params:
        timestampDir = lambda wildcards: config['timestamp'] + '-ppppOutputs'
    log:
        "cleanlogs/ab1_clean.log"
    shell:
        """
        if [ -d ab1 ]; then
            mv ab1/ {params.timestampDir}/ >> {log} 2>&1
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
        rules.ab1_clean.output
    params:
        timestampDir = lambda wildcards: config['timestamp'] + '-ppppOutputs'
    shell:
        """
        rm {input}
        if [ -d cleanlogs ]; then
            mv cleanlogs {params.timestampDir}
        fi
        zip -r -m {params.timestampDir}.zip {params.timestampDir}
        """
