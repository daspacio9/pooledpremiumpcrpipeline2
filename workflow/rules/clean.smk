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


onstart:
    from pathlib import Path
    from os.path import join, dirname, realpath

    here = os.getcwd()
    print(here)
    # 2. Construct the full path pattern
    search_pattern = join(here,'.*.done')
    print(f"Searching for .done files with pattern: {search_pattern}")

    files_to_remove = Path(here).glob('.*.done')
    files_to_remove = [str(f) for f in files_to_remove if f.is_file()]
    # Remove flag files to force re-run of cleanup
    if files_to_remove:
        for f in files_to_remove:
            os.remove(f)
    else:
        print("No .done flag files found to clean up.")

    
# clean up sequences folder: remove combined .fastq.gz files (but not basecalled batches), move UMI stats files
# -----------------------------------------------------
rule make_pppp_output_dir:
    output:
        touch('.make_pppp_output_dir.done')   # changed: consistent touch filename
    params:
        timestampDir = lambda wildcards: config['timestamp'] + '-ppppOutputs'
    shell:
        """
        mkdir -p {params.timestampDir}
        """

# clean up compute batches alignment
# -----------------------------------------------------
rule alignment_clean:
    output:
        touch('.alignment_clean.done')
    shell:
        """
        if [ -d aln ]; then
            rm -r aln
        fi
        """

# clean up compute batches demux
# -----------------------------------------------------
rule demux_clean:
    output:
        touch('.demux_clean.done')
    params:
        timestampDir = lambda wildcards: config['timestamp'] + '-ppppOutputs'
    shell:
        """
        if [ -d demux ]; then
            rm -r demux
        fi
        """

# clean up consensus outputs
# -----------------------------------------------------
rule consensus_clean:
    output:
        touch('.consensus_clean.done')
    params:
        timestampDir = lambda wildcards: config['timestamp'] + '-ppppOutputs'
    shell:
        """
        if [ -d consensus ]; then
            mkdir -p {params.timestampDir}/consensus
            mv consensus/*.fastq {params.timestampDir}/consensus/
            rm -r consensus
        fi
        """

# clean up logs
# -----------------------------------------------------
rule logs_clean:
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
        cp config/*.yaml {params.timestampDir}/
        cp -r ref {params.timestampDir}/
        """

# clean up reports
# -----------------------------------------------------
rule report_clean:
    output:
        touch('.report_clean.done')
    params:
        timestampDir = lambda wildcards: config['timestamp'] + '-ppppOutputs'
    shell:
        """
        if [ -d report ]; then
            mv report {params.timestampDir}
        fi
        if [ -f demux_stats.csv ]; then
            mv demux_stats.csv {params.timestampDir}
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
    shell:
        """
        rm {input}
        zip -r -m {params.timestampDir}.zip {params.timestampDir}
        """
