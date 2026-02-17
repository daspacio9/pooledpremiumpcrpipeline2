import pandas as pd
import re
from datetime import datetime


def now():
    """Return the current date and time as a string formatted as YY-MM-DD-HH-MM-SS."""
    return datetime.now().strftime("%y-%m-%d-%H%M%S")


def reverse_complement(seq):
            comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
            return seq.translate(comp)[::-1]



def parse_reads_row(read_str, ref_base):
    """Parse the complex pileup read string for a single row."""
    counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'indel': 0, 'ref_match': 0}
    i = 0
    while i < len(str(read_str)):
        char = read_str[i]
        if char in ['.', ',']:
            # Base matches the reference base
            counts['ref_match'] += 1
            counts[ref_base.upper()] += 1
        elif char in ['A', 'C', 'G', 'T', 'a', 'c', 'g', 't']:
            # Mismatch
            counts[char.upper()] += 1
        elif char in ['+', '-']:
            # Indel, skip over the indel length and sequence
            counts['indel'] += 1
            i += 1
            indel_len_str = ""
            while i < len(read_str) and read_str[i].isdigit():
                indel_len_str += read_str[i]
                i += 1
            indel_len = int(indel_len_str)
            i += indel_len - 1
        elif char in ['^', '$', '*']:
            # Read start (^), read end ($), or deletion placeholder (*)
            if char == '^':
                i += 1  # Skip the mapping quality character
            if char == '*':
                counts['indel'] += 1
        i += 1
    return counts

def parse_mpileup(file_path):
    """
    Parses a mpileup file and returns a pandas DataFrame with counts of A, C, G, T,
    ref matches, and indel events at each position.
    """
    column_names = ['chrom', 'pos', 'ref', 'coverage', 'reads', 'qualities']
    try:
        mpileup_df = pd.read_csv(file_path, sep='\t', header=None, names=column_names)
    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
        return None

    # Add columns for each base, indels, and reference matches
    for base in ['A', 'C', 'G', 'T']:
        mpileup_df[base] = 0
    mpileup_df['indel'] = 0
    mpileup_df['ref_match'] = 0

    # Iterate through the DataFrame and parse the read strings for each row
    for row_index, row in mpileup_df.iterrows():
        read_counts = parse_reads_row(row['reads'], row['ref'])
        for base in ['A', 'C', 'G', 'T']:
            mpileup_df.at[row_index, base] = read_counts.get(base, 0)
        mpileup_df.at[row_index, 'indel'] = read_counts.get('indel', 0)
        mpileup_df.at[row_index, 'ref_match'] = read_counts.get('ref_match', 0)

    return mpileup_df


# Example usage with your pileup file

# df_pileup = parse_mpileup('cole1-group-C1_mpileup.txt')
# df_pileup.to_csv("cole1-group-C1_pypileup.tsv", sep='\t', index=False)

import struct
import math

class AbifWriter:
    def __init__(self, filename):
        self.filename = filename
        self.entries = []
        self.data_buffer = bytearray()
        # Data starts immediately after the 128-byte header
        self.current_offset = 128
        
    def add_tag(self, name, number, element_type, element_size, num_elements, data):
        """
        Adds a directory entry.
        """
        # Name must be 4 bytes
        name_bytes = name.encode('ascii') if isinstance(name, str) else name
        data_size = len(data)
        
        # ABIF Rule: Data <= 4 bytes is stored DIRECTLY in the offset field
        if data_size <= 4:
            # Pad to 4 bytes (left-aligned)
            padded_data = data.ljust(4, b'\x00')
            offset_val = struct.unpack('>I', padded_data)[0]
            entry = struct.pack('>4s I h h I I I I',
                name_bytes, number, element_type, element_size, num_elements, data_size, offset_val, 0)
            self.entries.append(entry)
        else:
            # Data > 4 bytes is stored in the data block
            entry = struct.pack('>4s I h h I I I I',
                name_bytes, number, element_type, element_size, num_elements, data_size, self.current_offset, 0)
            self.entries.append(entry)
            self.data_buffer.extend(data)
            self.current_offset += data_size

    def write(self):
        with open(self.filename, 'wb') as f:
            # 1. Header: Signature + Version
            f.write(b'ABIF')
            f.write(struct.pack('>H', 101))
            
            # 2. Root Directory Entry (Located at OFFSET 6)
            # This points to the directory we will write at the end of the file
            dir_location = self.current_offset
            tag_entry_size = 28
            root_entry = struct.pack('>4s I h h I I I I',
                b'tdir', 1, 1023, tag_entry_size, len(self.entries),
                len(self.entries) * tag_entry_size, dir_location, 0)
            f.write(root_entry)
            
            # 3. Header Padding (Fill up to 128 bytes)
            f.write(b'\x00' * (128 - f.tell()))
            
            # 4. Data Block
            f.write(self.data_buffer)
            
            # 5. Directory (The map of tags)
            for entry in self.entries:
                f.write(entry)

### AB1 Writing Function
def write_ab1(filename, sequence, trace_g, trace_a, trace_t, trace_c):
    """
    Writes a valid ABIF file from raw sequence and trace lists.
    Adheres to the 4-point-per-base density required by the sample format.
    """
    writer = AbifWriter(filename)
    
    seq_len = len(sequence)
    trace_len = len(trace_g)
    
    # --- 1. VALIDATION ---
    # The sample file format STRICTLY uses 4 points per base.
    expected_len = seq_len * 4
    if trace_len != expected_len:
        print(f"Warning: Trace length ({trace_len}) is not 4x sequence length ({expected_len}). "
              "This may cause 'Low Quality' errors in Benchling.")

    # --- 2. PACKING DATA (Big-Endian) ---
    
    # Traces: Type 4 (Short / 2 bytes)
    packed_g = struct.pack(f'>{trace_len}h', *trace_g)
    packed_a = struct.pack(f'>{trace_len}h', *trace_a)
    packed_t = struct.pack(f'>{trace_len}h', *trace_t)
    packed_c = struct.pack(f'>{trace_len}h', *trace_c)
    
    # Sequence: Type 2 (Char / 1 byte)
    packed_seq = sequence.encode('ascii')
    
    # Quality Scores: Type 2 (Char / 1 byte) <-- CRITICAL FIX
    # Previous attempts failed because we used Type 1 (Byte). The sample uses Type 2.
    # We default to 50 (High Quality).
    packed_pcon = struct.pack(f'>{seq_len}B', *[50]*seq_len)
    
    # Peak Locations: Type 4 (Short / 2 bytes)
    # For 4-point density, the peak is usually at index 2 of the 0-3 block.
    ploc = [int(i * 4 + 2) for i in range(seq_len)]
    packed_ploc = struct.pack(f'>{seq_len}h', *ploc)

    # --- 3. ADDING TAGS ---
    # The order and Types must match the working sample exactly.
    
    # PBAS 1: Sequence (Type 2, Size 1)
    writer.add_tag('PBAS', 1, 2, 1, seq_len, packed_seq)
    
    # PLOC 1: Peak Locations (Type 4, Size 2)
    writer.add_tag('PLOC', 1, 4, 2, seq_len, packed_ploc)
    
    # PBAS 2: Sequence Copy (Type 2, Size 1)
    writer.add_tag('PBAS', 2, 2, 1, seq_len, packed_seq)
    
    # DATA 9-12: Analyzed Traces (Type 4, Size 2)
    writer.add_tag('DATA', 9,  4, 2, trace_len, packed_g)
    writer.add_tag('DATA', 10, 4, 2, trace_len, packed_a)
    writer.add_tag('DATA', 11, 4, 2, trace_len, packed_t)
    writer.add_tag('DATA', 12, 4, 2, trace_len, packed_c)
    
    # FWO_ 1: Base Order (Type 2, Size 1)
    writer.add_tag('FWO_', 1, 2, 1, 4, b'GATC')
    
    # PCON 1: Quality Scores (Type 2, Size 1) <-- CRITICAL FIX
    writer.add_tag('PCON', 1, 2, 1, seq_len, packed_pcon)
    
    # PCON 2: Quality Scores Copy (Type 2, Size 1) <-- CRITICAL FIX
    writer.add_tag('PCON', 2, 2, 1, seq_len, packed_pcon)
    
    # --- 4. WRITE FILE ---
    writer.write()
    #print(f"Successfully generated {filename}")

def parse_pypileup(file_path, cols = ['A', 'C', 'G', 'T']):
    """Reads the pypileup TSV file and returns a DataFrame with just the A, C, G, T columns."""
    df = pd.read_csv(file_path, sep='\t')
    df_filtered = df[cols]
    
    # Find the max within just those 4 columns
    subset_max = df_filtered.max().max()

    # Update the columns in place
    df_norm = (df_filtered / subset_max)
    #print(df_norm[0:6])
    dfmax = df_filtered.idxmax(axis=1)
    cons_seq = "".join(dfmax.astype(str)).strip()
    return df_norm, cons_seq

def generate_trace(df, seq):
    """Generates a trace coordinate list for each base from a pypileup basecounts parse, with a spike at the peak location (index 2 of the block) and 0s on either side."""
    total_pts = len(seq) * 4
    g, a, t, c = [0]*total_pts, [0]*total_pts, [0]*total_pts, [0]*total_pts

    for i, _ in enumerate(seq):
        # Create a spike at the peak location (index 2 of the block)
        peak_idx = i * 4 + 2

        # Simple triangle shape 0-1000-0
        val_peak = 1000
        val_side = 0
        
        # Scale the peak value by the normalized base count for that position and base
        a[peak_idx] = int(val_peak * df["A"][i])
        g[peak_idx] = int(val_peak * df["G"][i])
        c[peak_idx] = int(val_peak * df["C"][i])
        t[peak_idx] = int(val_peak * df["T"][i])

        a[peak_idx+1] = int(val_peak * df["A"][i])
        g[peak_idx+1] = int(val_peak * df["G"][i])
        c[peak_idx+1] = int(val_peak * df["C"][i])
        t[peak_idx+1] = int(val_peak * df["T"][i])
        
        a[peak_idx-1] = a[peak_idx+1] = val_side
        g[peak_idx-1] = g[peak_idx+1] = val_side
        c[peak_idx-1] = c[peak_idx+1] = val_side
        t[peak_idx-1] = t[peak_idx+1] = val_side
    return {'G': g, 'A': a, 'T': t, 'C': c}
# 3. Write
#write_ab1("generated_from_scratch.ab1", test_seq, g, a, t, c)