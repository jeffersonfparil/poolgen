import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
from pathlib import Path
import sys
warnings.filterwarnings("ignore")
plt.style.use('ggplot')

plt.rcParams.update({
    "font.size": 10,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 7,
    "figure.titlesize": 10,
})

filename = sys.argv[1]

gwas = pd.read_csv(filename, index_col=0)
gwas_phenotypes = {name: group for name, group in gwas.groupby('phenotype')}

output_paths = ""

for gwas_key in gwas_phenotypes:
    gwas = gwas_phenotypes[gwas_key]
    gwas["chromosome"] = gwas.index
    gwas = gwas[gwas["chromosome"] != "intercept"]
    gwas["position"] = gwas["pos"]
    gwas["log_pvalue"] = -np.log10(gwas["pvalue"])

    chromosomes = sorted(gwas["chromosome"].unique())
    color_map = {chrom: plt.cm.tab10.colors[i] if i < 10 else "gray"
        for i, chrom in enumerate(chromosomes)}

    gwas = gwas.reset_index(drop=False)

    x_ticks = []
    x_labels = []
    current_position = 0

    plt.figure(figsize=(10, 3))

    for i, chrom in enumerate(chromosomes):
        chrom_data = gwas[gwas["chromosome"] == chrom].copy()
        chrom_data["x"] = chrom_data["position"] + current_position
        
        plt.scatter(chrom_data["x"], chrom_data["log_pvalue"], 
            color=color_map[chrom], alpha=0.7, label=chrom, s=12)
        
        mid_position = chrom_data["x"].median()
        x_ticks.append(mid_position)
        x_labels.append(chrom)
        current_position = chrom_data["x"].max() + 1e6

    sig_threshold = 0.05/len(gwas) # bonferonni correction
    plt.axhline(y=-np.log10(sig_threshold), color="red", linestyle="--")

    plt.xlabel("Chromosome")
    plt.ylabel("$-\\log_{10}$(p-value)")
    plt.title("Manhattan plot (" + Path(filename).stem + "," + gwas_key + ")")
    plt.xticks(x_ticks, x_labels, rotation=90)
    max_legend_items = 10

    handles, labels = plt.gca().get_legend_handles_labels()

    if len(labels) > max_legend_items:
        from matplotlib.lines import Line2D
        handles = handles[:max_legend_items - 1]
        labels = labels[:max_legend_items - 1]
        handles.append(Line2D([0], [0], color='gray', linestyle='dotted'))
        labels.append(f"... {len(chromosomes) - (max_legend_items - 1)} more")

    plt.legend(handles, labels, title="Chromosome", bbox_to_anchor=(1.05, 1), loc="upper left")

    plt.tight_layout()
    output_path = Path(filename).stem + "_" + gwas_key + "_manhattan.png"
    output_paths += output_path + "\n"
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches='tight', dpi=300)

output_paths = output_paths[:-1]
print("File created: " + output_paths, end="")
