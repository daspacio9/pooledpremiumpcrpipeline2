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

