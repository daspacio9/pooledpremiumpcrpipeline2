# Helper to get all ab1 files for all passed samples and all chunks using the chunk_pypileup checkpoint
def get_all_ab1_files(wildcards):
    ab1_files = []
    # Get all samples that passed
    for sample in get_passed_samples(wildcards):
        # Use the chunk_pypileup checkpoint to get the output directory
        chunk_dir = os.path.join("report", "chunks")
        import glob
        chunk_files = glob.glob(os.path.join(chunk_dir, f"{sample}_pypileup_chunk_*.tsv"))
        for chunk_file in chunk_files:
            chunk_id = os.path.basename(chunk_file).split('_chunk_')[-1].replace('.tsv', '')
            ab1_files.append(f"ab1/{sample}_{chunk_id}.ab1")
    return ab1_files
# \HEADER\-------------------------------------------------------------------------
#
#  CONTENTS      : Snakemake nanopore data pipeline
#
#  DESCRIPTION   : Common rules and functions for the pipeline
#
#  RESTRICTIONS  : none
#
#  REQUIRES      : none
#
# ---------------------------------------------------------------------------------

# Function that depends on low depth checkpoint to determine which files to move
# -----------------------------------------------------
def get_moved_files(wildcards):
    # This forces Snakemake to wait until the checkpoint is done
    checkpoint_output = checkpoints.move_low_depth_subreads.get(**wildcards).output[0]
    # Return the expected paths in the final destination
    return [os.path.join("demux/.low_depth", f) for f in os.listdir(checkpoint_output) if not f.startswith('.')]

# Function that overrides the default reference file for debugging
# -----------------------------------------------------
def use_debug():
    if config["debug"]:
        return "ref/barcode-pairs_debug.txt"
    else:
        return "ref/barcode-pairs.txt"

# Function to get passed samples after demux checkpoint is complete
# -----------------------------------------------------
def get_passed_samples(wildcards): #wildcards is required argument for all snakemake dynamic input/output objects to be accessed in the function
# This function reads the demux_stats.csv file to get samples that passed the depth filter
    passed_samples = []
    demux_stats_file = checkpoints.demux_stats.get().output[0] # get the output file of the checkpoint
    with open(demux_stats_file) as f:
        next(f)  # skip header
        for line in f:
            sample, n_reads = line.strip().split(",")
            if sample != "unknown" and int(n_reads) >= config["min_depth"]:
                passed_samples.append(sample)
    return passed_samples

# Function to read adapter (primer) names from the barcode-pairs.txt file. Functions called in input need to be in the snakefile or smk file that uses them.
# -----------------------------------------------------
def adapter_names(path):
    names = []
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith(("#",";")):
                continue
            # take the token after '>'
            m = s.startswith(">")
            if m:
                names.append(s[1:])
    return names

# Helper function that determines which files are produced by the split_fastq checkpoint and passes them to snakemake only after files are made
# The specific files to pass to the downstream rules can be specified using the prefix and suffix arguments
# -----------------------------------------------------
def get_splits(wildcards, prefix= '', suffix = ''):
    # 1. Ensure the checkpoint has finished
    checkpoint_output = checkpoints.split_fastq.get(**wildcards).output[0]
    # 2. List the files that were actually created
    # This glob_wildcards looks inside the directory produced by the checkpoint
    filenames = glob_wildcards(os.path.join(checkpoint_output, "{sample}_consensus.fastq")).sample
    #print(filenames)
    # 3. Return the full paths to the next rule
    return expand(os.path.join(checkpoint_output, f"{prefix}{{sample}}{suffix}"), sample=filenames)

# Helper to get all ab1 files for all passed samples and all chunks
def get_all_ab1_files(wildcards):
    import glob
    ab1_files = []
    for sample in get_passed_samples(wildcards):
        # Find all chunk files for this sample
        chunk_pattern = f"report/chunks/{sample}_pypileup_chunk_*.tsv"
        chunk_files = glob.glob(chunk_pattern)
        for chunk_file in chunk_files:
            chunk_id = chunk_file.split('_chunk_')[-1].replace('.tsv', '')
            ab1_files.append(f"ab1/{sample}_{chunk_id}.ab1")
    return ab1_files


