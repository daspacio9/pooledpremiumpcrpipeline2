# Plan: Fix Snakemake Linting & Clean Up README

## Overview
Update README with real project info (remove template boilerplate), add rule-specific conda environments for all 15+ rules missing them, create log directives for 3 rules, migrate 7 complex run directives to separate Python scripts, and fix the path concatenation in consensus.smk. This improves code maintainability, reproducibility, and Snakemake best practices compliance.

## Implementation Steps

### 1. Clean up README & documentation
- Replace template placeholders in README.md (`<owner>`, `<repo>`, `<name>`, `<description>`)
- Update workflow description with actual purpose (demultiplexing + consensus generation)
- Fix working directory documentation section
- Remove references to unused features (simulate_reads, get_genome examples)
- Update config/README.md if it contains template boilerplate

**Status:** ✅ COMPLETED

### 2. Create rule-specific conda environments
Create 6 new YAML files in `workflow/envs/`:

- **clean.yaml** – no external tools needed (shell commands only, Python built-in)
- **demux.yaml** – cutadapt, seqkit
- **consensus.yaml** – medaka, biopython, pysam
- **align.yaml** – minimap2, samtools, bwa, pandas, matplotlib
- **ab1writer.yaml** – biopython, pysam (with common.py module)
- **prepare_ref.yaml** – biopython, python (for prepare_barcode-pairs.py)

**Status:** ✅ COMPLETED

### 3. Add conda directives to clean.smk rules
File: `workflow/rules/clean.smk`

Add `conda: "../envs/clean.yaml"` to these 6 rules:
- Line 25: `make_pppp_output_dir`
- Line 39: `alignment_clean`
- Line 55: `demux_clean`
- Line 73: `consensus_clean`
- Line 94: `logs_clean`
- Line 114: `report_clean`
- Line 132: `ab1_clean`
- Line 149: `clean`

**Status:** ⏳ NOT STARTED

### 4. Add conda directives to demux.smk rules
File: `workflow/rules/demux.smk`

Add appropriate conda directives to:
- Line 15: `filter_reads_by_length` → `conda: "../envs/demux.yaml"`
- Line 31: `cutadapt_demux_linked` → `conda: "../envs/demux.yaml"`
- Line 112: `move_file` → `conda: "../envs/demux.yaml"`
- Line 124: `finalize_demux` → `conda: "../envs/demux.yaml"`

**Status:** ⏳ NOT STARTED

### 5. Add conda directives to other rules
**consensus.smk:**
- Line 42: `medaka_consensus_from_subreads` → `conda: "../envs/consensus.yaml"`

**align.smk:**
- Line 16: `aln_to_consensus` → `conda: "../envs/align.yaml"`

**ab1writer.smk:**
- Rules without conda spec → `conda: "../envs/ab1writer.yaml"`

**prepare_reference.smk:**
- Line 40: `prepare_barcode_reference` → `conda: "../envs/prepare_ref.yaml"`

**Status:** ⏳ NOT STARTED

### 6. Add log directives to 3 rules
Add `log: "logs/{description}/{rule_name}.log"` pattern:

**clean.smk:**
- Line 94: `logs_clean` → `log: "logs/clean/logs_clean.log"`
- Line 149: `clean` → `log: "logs/clean/clean.log"`

**prepare_reference.smk:**
- Line 40: `prepare_barcode_reference` → `log: "logs/prepare_reference/barcode_reference.log"`

**Status:** ⏳ NOT STARTED

### 7. Migrate 7 long run directives to scripts
Create new Python scripts in `workflow/scripts/` and update rules:

| Rule | Source | New Script | Location |
|------|--------|-----------|----------|
| `fastq_to_fastq_subreads` | consensus.smk:16 | `fastq_to_fastq_subreads.py` | Line 16 |
| `split_fastq` | consensus.smk:74 | `split_fastq.py` | Line 74 |
| `consensus_summary_csv` | consensus.smk:104 | `consensus_summary_csv.py` | Line 104 |
| `demux_stats` | demux.smk:60 | `demux_stats.py` | Line 60 |
| `move_low_depth_subreads` | demux.smk:87 | `move_low_depth_subreads.py` | Line 87 |
| `plot_coverage` | align.smk:54 | `plot_coverage.py` | Line 54 |
| `chunk_pypileup` | ab1writer.smk:14 | `chunk_pypileup.py` | Line 14 |

Each migrated script must:
- Use `snakemake.input`, `snakemake.output`, `snakemake.wildcards`, `snakemake.params`, `snakemake.log` as needed
- Preserve all original logic
- Include docstring explaining purpose

**Status:** ⏳ NOT STARTED

### 8. Fix path composition in consensus.smk
**File:** `workflow/rules/consensus.smk`
**Line:** 92

Replace string concatenation with f-string or `pathlib.Path`:
- Current: Uses path composition with `+` operator
- Fix: Use f-string formatting or pathlib for readability

**Status:** ⏳ NOT STARTED

### 9. Optional cleanup
Consider for future (not blocking):
- Remove or archive unused environment files:
  - `workflow/envs/get_genome.yaml`
  - `workflow/envs/simulate_reads.yaml`
  - `workflow/envs/validate_genome.yaml`

**Status:** ⏳ DEFERRED

## Verification Checklist

- [ ] Run `snakemake -s workflow/Snakefile -j 4 --lint` → expect 0 warnings
- [ ] Check all 6 new conda env files exist and have valid YAML syntax
- [ ] Verify all 7 migrated scripts can access snakemake object properties
- [ ] Test workflow dry-run: `snakemake -s workflow/Snakefile -n`
- [ ] Check git status for untracked files (new scripts, envs)

## Decision Summary

- **Conda strategy:** Rule-specific environments for better isolation and documented tool dependencies per-step
- **README focus:** Core fixes only (remove placeholders, update actual content); defer deep restructuring
- **Script migration:** All 7 run directives → separate scripts for readability and testing
