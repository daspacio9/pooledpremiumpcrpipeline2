# \HEADER\-------------------------------------------------------------------------
#
#  CONTENTS      : Snakemake nanopore data pipeline
#
#  DESCRIPTION   : Workflow rules for writing the ab1 files from the pypileup chunks. Ab1 writer from Vincent Hu.
#  RESTRICTIONS  : none
#
#  REQUIRES      : none
#
# ---------------------------------------------------------------------------------


### Chunking the pypileup file for the ab1 writer at 5 kb. ab1 format has a max length of ~8 kb according to plasmidsaurus, but 5kb is a nice round number.
# -----------------------------------------------------
checkpoint chunk_pypileup:
    input:
        "report/{sample}_pypileup.tsv",
    output:
        chunk_flag=touch("report/.chunks_{sample}.done"),
    log:
        "logs/aln/{sample}_chunk_pypileup.log",
    conda:
        "../envs/ab1writer.yaml"
    params:
        chunk_size=5000,
    script:
        "../scripts/chunk_pypileup.py"


rule write_ab1:
    input:
        chunk_flag="report/.chunks_{sample}.done",
        chunk_file="report/pileup_chunks/{sample}_pypileup_chunk_{chunk}.tsv",
    output:
        flag=touch("ab1/.ab1_{sample}_{chunk}.done"),
        ab1_file="ab1/{sample}_{chunk}.ab1",
    log:
        logf="logs/ab1writer/{sample}_{chunk}_write_ab1.log",
    conda:
        "../envs/ab1writer.yaml"
    script:
        "../scripts/write_ab1.py"


rule ab1_done:
    input:
        ab1_files=lambda wildcards: get_all_ab1_files(wildcards),
    output:
        touch(".ab1_done.done"),
    log:
        logf="logs/aln/ab1_done.log",
    conda:
        "../envs/ab1writer.yaml"
    run:
        with open(output[0], "w") as f:
            f.write("AB1 writing completed.\n")
