
"""
Demultiplexing rules and preparation of inputs by primer structure for cutadapt.
Combines barcode groups with primer contexts to create adapter pairs for cutadapt.
"""
import os
from Bio import SeqIO
import pandas as pd
from common import reverse_complement, log_msg
import re

input1 = snakemake.input[0]  # Input FASTA file with primer contexts
input2 = snakemake.input[1]  # Input CSV file with barcode groups
input3 = snakemake.input[2]  # Input FASTA file with barcode sequences
output1 = snakemake.output.out1  # Output text file with barcode pairs

with open(snakemake.log[0], "w") as logf:
    log_msg(logf, "Starting barcode pair preparation")
    log_msg(logf, f"Input primer contexts: {input1}")
    log_msg(logf, f"Input barcode groups: {input2}")
    log_msg(logf, f"Input barcode sequences: {input3}")
    log_msg(logf, f"Output adapter pairs: {output1}")
    
    # Read the barcode groups CSV file
    log_msg(logf, "Reading barcode groups CSV file")
    barcode_groups = pd.read_csv(input2, header='infer')
    log_msg(logf, f"Loaded {len(barcode_groups)} barcode group(s)")
    
    # Create a dictionary to map sample names to barcodes
    log_msg(logf, "Building sample to barcode mapping")
    sample_barcodes = {}
    for index, row in barcode_groups.iterrows():
        sample = row[0]
        barcodes = row[1:].dropna().tolist()
        
        # Validate that each sample has exactly 2 barcodes (forward and reverse)
        if len(barcodes) != 2:
            error_msg = f"ERROR: Sample '{sample}' has {len(barcodes)} barcode(s), expected exactly 2 (forward and reverse). Found: {barcodes}"
            log_msg(logf, error_msg)
            raise ValueError(error_msg)
        
        sample_barcodes[sample] = barcodes
        log_msg(logf, f"  Sample {sample}: barcodes {barcodes}")
    
    # Read the input FASTA file with primer contexts
    log_msg(logf, "Reading primer context FASTA file")
    records = list(SeqIO.parse(input1, "fasta"))
    log_msg(logf, f"Loaded {len(records)} primer record(s)")
    
    # Validate we have at least 2 records (forward and reverse primer contexts)
    if len(records) < 2:
        error_msg = f"ERROR: Expected at least 2 primer records in {input1}, but found {len(records)}"
        log_msg(logf, error_msg)
        raise ValueError(error_msg)
    
    # Extract the primer sequences
    f_context = str(records[0].seq)
    r_context = str(records[1].seq)
    log_msg(logf, f"Forward primer context length: {len(f_context)} bp")
    log_msg(logf, f"Reverse primer context length: {len(r_context)} bp")
    
    # Split on one or more Ns to extract flanking sequences
    log_msg(logf, "Splitting primer contexts on N positions")
    f_parts = [p for p in re.split(r"(?i)N+", f_context) if p]
    r_parts = [p for p in re.split(r"(?i)N+", r_context) if p]
    log_msg(logf, f"Forward flanking parts: {len(f_parts)}, Reverse flanking parts: {len(r_parts)}")
    
    # Validate that we extracted at least one flanking part from each primer
    if not f_parts:
        error_msg = f"ERROR: Forward primer context contains no flanking sequences after splitting on N positions: {f_context}"
        log_msg(logf, error_msg)
        raise ValueError(error_msg)
    if not r_parts:
        error_msg = f"ERROR: Reverse primer context contains no flanking sequences after splitting on N positions: {r_context}"
        log_msg(logf, error_msg)
        raise ValueError(error_msg)
    
    # Read barcode sequences from FASTA file
    log_msg(logf, "Reading barcode sequences FASTA file")
    records = list(SeqIO.parse(input3, "fasta"))
    barcode_seqs = {record.id: str(record.seq) for record in records}
    log_msg(logf, f"Loaded {len(barcode_seqs)} barcode sequence(s)")
    
    # Write the adapter pairs to a text file
    log_msg(logf, "Writing adapter pairs to output file")
    pair_count = 0
    with open(output1, "w") as f:
        for s in sample_barcodes:
            fkey = sample_barcodes[s][0]
            rkey = sample_barcodes[s][1]
            
            if fkey not in barcode_seqs or rkey not in barcode_seqs:
                log_msg(logf, f"WARNING: Missing barcode sequence for {s} ({fkey} or {rkey})")
                continue
            
            f_barcode = barcode_seqs[fkey]
            r_barcode_rev = reverse_complement(barcode_seqs[rkey])
            f_flank = f_parts[0]
            r_flank_rev = reverse_complement(r_parts[0])
            
            adapter_pair = f"^{f_flank}{f_barcode}...{r_barcode_rev}{r_flank_rev}$"
            f.write(f">{s}\n")
            f.write(f"{adapter_pair}\n")
            
            log_msg(logf, f"  Sample {s}: {fkey}+{rkey} adapter pair written")
            pair_count += 1
    
    log_msg(logf, f"Barcode pair preparation complete. Created {pair_count} adapter pair(s)")
