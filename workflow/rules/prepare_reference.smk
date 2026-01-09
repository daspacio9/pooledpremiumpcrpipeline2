# \HEADER\-------------------------------------------------------------------------
#
#  CONTENTS      : Snakemake nanopore data pipeline
#
#  DESCRIPTION   : prepare the references files for cutadapt when using different
#                   primers to barcode samples
#
#  RESTRICTIONS  : none
#
#  REQUIRES      : none
#
# ---------------------------------------------------------------------------------


import sys
import os
from pathlib import Path

# Get the base directory of the current Snakefile
# -----------------------------------------------------
snakefile_basedir = Path(workflow.basedir)

# add PATHs relative to Snakefile:
# -----------------------------------------------------
relative_path = "scripts" 

absolute_path = f"{snakefile_basedir}/{relative_path}"

# Add the absolute path to sys.path
# -----------------------------------------------------
if str(absolute_path) not in sys.path:
    sys.path.append(absolute_path)

import common
configfile: "config/config.yaml"

### Demultiplexing rules and preparation of inputs by primer structure for cutadapt
# -----------------------------------------------------

rule prepare_barcode_reference:
    input: f"ref/{config['primer_consensus']}",
            f"ref/{config['barcode_groups']}",
            f"ref/{config['barcodes']}"
    output: out1 = "ref/barcode-pairs.txt"
    script: 
        "../scripts/prepare_barcode-pairs.py"



    

