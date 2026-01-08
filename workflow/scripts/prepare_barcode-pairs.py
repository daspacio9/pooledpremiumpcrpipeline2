
### Demultiplexing rules and preparation of inputs by primer structure for cutadapt
input1 = snakemake.input[0]  # Input FASTA file with primer contexts
input2 = snakemake.input[1]  # Input CSV file with barcode groups
input3 = snakemake.input[2]  # Input FASTA file with barcode sequences
output1 = snakemake.output.out1  # Output text file with barcode pairs
import os
from Bio import SeqIO
import pandas as pd
from common import reverse_complement
import re
# Read the barcode groups CSV file
barcode_groups = pd.read_csv(input2, header='infer')
# Create a dictionary to map sample names to barcodes
sample_barcodes = {}
for index, row in barcode_groups.iterrows():
    sample = row[0]
    barcodes = row[1:].dropna().tolist()
    sample_barcodes[sample] = barcodes

# Read the input FASTA file
records = list(SeqIO.parse(input1, "fasta"))

# Extract the primer sequences
f_context = str(records[0].seq)
r_context = str(records[1].seq)

# Split on one or more Ns
f_parts = [p for p in re.split(r"(?i)N+", f_context) if p]   # remove empties
r_parts = [p for p in re.split(r"(?i)N+", r_context) if p]
print(f_parts, r_parts) 
print(f_context, r_context)

#replace the barcode name with the barcode sequence in barcode pairs file
records = list(SeqIO.parse(input3, "fasta"))
barcode_seqs = {record.id: str(record.seq) for record in records}
#print(barcode_seqs["M527"])

# Write the adapter pairs to a text file
with open(output1, "w") as f:
    for s in sample_barcodes:
        fkey = sample_barcodes[s][0]
        rkey = sample_barcodes[s][1]
        print(barcode_seqs[fkey], barcode_seqs[rkey])
        f.write(f">{s}\n")
        f.write(f"^{f_parts[0]}{barcode_seqs[fkey]}...{reverse_complement(barcode_seqs[rkey])}{reverse_complement(r_parts[0])}$\n")
