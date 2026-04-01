"""
Generate coverage and mismatch frequency plots from pypileup data.
Creates two PDF files: coverage plot and substitution frequency plot.
"""
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Read reference match data
df_ref = pd.read_csv(snakemake.input[0], sep="\t")
df_ref['freq_mismatch'] = (df_ref['coverage'] - df_ref['ref_match']) / df_ref['coverage']
df_ref['freq_mismatch'] = df_ref['freq_mismatch'].fillna(1)  # Handle NaN values if coverage is 0

# Plot coverage
plt.figure(figsize=(10, 4))
plt.plot(df_ref["pos"], df_ref["coverage"], color='blue', linewidth=0.2)
plt.plot(df_ref["pos"], df_ref["ref_match"], color='green', linewidth=0.4, alpha=0.5)
plt.title(f"Coverage Plot for {snakemake.wildcards.sample}")
plt.xlabel("Position")
plt.ylabel("Depth")
plt.yscale("linear")
plt.ylim(0, round(df_ref["coverage"].max()*1.1))  # Add some headroom
plt.legend(["Coverage", "Ref Match"])
plt.grid(False)
plt.savefig(snakemake.output[0])
plt.close()

# Plot frequency of mismatches
plt.figure(figsize=(10, 4))
plt.plot(df_ref["pos"], df_ref["freq_mismatch"], color='red', linewidth=0.5)
plt.fill_between(df_ref["pos"], df_ref["freq_mismatch"], color='red', alpha=0.3)
plt.title(f"Substitution Plot for {snakemake.wildcards.sample}")
plt.xlabel("Position")
plt.ylabel("freq")
plt.yscale("linear")
plt.ylim(0, (df_ref["freq_mismatch"]).max()*1.1)
plt.legend(["freq_mismatch"])
plt.grid(False)
plt.savefig(snakemake.output[1])
plt.close()

with open(snakemake.output.plot_flag, "w") as _f:
    _f.write("Plotting completed.\n")