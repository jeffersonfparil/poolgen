import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import warnings
from pathlib import Path 
warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser(prog="Manhattan plotter")
parser.add_argument("filename")
args = parser.parse_args()

gwas = pd.read_csv(args.filename, index_col=0)
gwas_phenotypes = {name: group for name, group in gwas.groupby('phenotype')}

output_paths = ""

for gwas_key in gwas_phenotypes:
    gwas = gwas_phenotypes[gwas_key]
    gwas["chromosome"] = gwas.index
    gwas = gwas[gwas["chromosome"] != "intercept"]
    gwas["position"] = gwas["pos"]
    gwas["log_pvalue"] = -np.log10(gwas["pvalue"])

    chromosomes = sorted(gwas["chromosome"].unique())
    colors = plt.cm.tab10(np.linspace(0, 1, len(chromosomes)))

    gwas = gwas.reset_index(drop=False)

    x_ticks = []
    x_labels = []
    current_position = 0

    plt.figure(figsize=(10, 3))

    for i, chrom in enumerate(chromosomes):
        chrom_data = gwas[gwas["chromosome"] == chrom].copy()
        chrom_data["x"] = chrom_data["position"] + current_position
        sns.scatterplot(x=chrom_data["x"], y=chrom_data["log_pvalue"], 
                    color=colors[i], alpha=0.5, label=chrom, linewidth = 0)
        mid_position = chrom_data["x"].median()
        x_ticks.append(mid_position)
        x_labels.append(chrom)
        current_position = chrom_data["x"].max() + 1e6

    sig_threshold = 0.05/len(gwas) # bonferonni correction
    plt.axhline(y=-np.log10(sig_threshold), color="red", linestyle="--")

    plt.xlabel("Chromosome")
    plt.ylabel("$-\\log_{10}$(p-value)")
    plt.title("Manhattan plot (" + Path(args.filename).stem + "," + gwas_key + ")")
    plt.xticks(x_ticks, x_labels, rotation=90)
    plt.legend(title="Chromosome", bbox_to_anchor=(1.05, 1), loc="upper left")

    plt.tight_layout()
    output_path = Path(args.filename).stem + "_" + gwas_key + "_manhattan.png"
    output_paths += output_path + "\n"
    plt.savefig(output_path, dpi=300)

output_paths = output_paths.removesuffix("\n")
print(output_paths, end="")
