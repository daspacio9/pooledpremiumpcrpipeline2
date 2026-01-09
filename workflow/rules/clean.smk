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
        "logs/clean/make_pppp_output_dir.log"
    params:
        timestampDir = lambda wildcards: config['timestamp'] + '-ppppOutputs'
    shell:
        """
        mkdir -p {params.timestampDir} >> {log} 2>&1
        """

# clean up compute batches alignment
# -----------------------------------------------------
rule alignment_clean:
    output:
        touch('.alignment_clean.done')
    log:
        "logs/clean/alignment_clean.log"
    shell:
        """
        if [ -d aln ]; then
            rm -rv aln >> {log} 2>&1
        fi
        """

# clean up compute batches demux
# -----------------------------------------------------
rule demux_clean:
    output:
        touch('.demux_clean.done')
    log:
        "logs/clean/demux_clean.log"
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
    output:
        touch('.consensus_clean.done')
    log:
        "logs/clean/consensus_clean.log"
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
    output:
        touch('.logs_clean.done')
    params:
        timestampDir = lambda wildcards: config['timestamp'] + '-ppppOutputs'
    log:
        "logs/clean/logs_clean.log"
    shell:
        """
        if [ -d logs ]; then
            mkdir -p {params.timestampDir}/logs >> {log} 2>&1
            mv logs {params.timestampDir}/ >> {log} 2>&1
        fi
        cp config/*.yaml {params.timestampDir}/ >> {log} 2>&1
        cp -r ref {params.timestampDir}/ >> {log} 2>&1
        """

# clean up reports
# -----------------------------------------------------
rule report_clean:
    output:
        touch('.report_clean.done')
    params:
        timestampDir = lambda wildcards: config['timestamp'] + '-ppppOutputs'
    log:
        "logs/clean/report_clean.log"
    shell:
        """
        if [ -d report ]; then
            mv report {params.timestampDir} >> {log} 2>&1
        fi
        if [ -f demux_stats.csv ]; then
            mv demux_stats.csv {params.timestampDir} >> {log} 2>&1
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
        rules.report_clean.output
    params:
        timestampDir = lambda wildcards: config['timestamp'] + '-ppppOutputs'
    log:
        "logs/clean/clean.log"
    shell:
        """
        rm {input}
        zip -r -m {params.timestampDir}.zip {params.timestampDir} >> {log} 2>&1
        """
