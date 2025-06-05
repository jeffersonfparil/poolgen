import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
from pathlib import Path 
from scipy.stats import kstest
import sys
warnings.filterwarnings("ignore")
plt.style.use('ggplot')

filename = sys.argv[1]

gwas = pd.read_csv(filename, index_col=0)
gwas = gwas.reset_index(drop=True)
gwas_phenotypes = {name: group for name, group in gwas.groupby('phenotype')}

output_paths = ""

for gwas_key in gwas_phenotypes:
    gwas = gwas_phenotypes[gwas_key]

    pvalues = gwas["pvalue"].sort_values().values
    n = len(pvalues)
    expected = -np.log10(np.linspace(1/n, 1, n))
    observed = -np.log10(pvalues)

    ks_stat, ks_pvalue = kstest(pvalues, 'uniform')
    ks_text = f"Kolmogorov-Smirnov statistic: {ks_stat:.4f}"

    plt.figure(figsize=(6, 6))
    plt.scatter(expected, observed, alpha=0.7)
    plt.plot([0, max(expected)], [0, max(expected)], color="red", linestyle="--")
    plt.xlabel("Expected $-\\log_{10}p$")
    plt.ylabel("Observed $-\\log_{10}p$")
    plt.title("Q-Q Plot of GWAS p-values uniformity (" + Path(filename).stem + "," + gwas_key + ")")
    plt.text(0.05, 0.95, ks_text, transform=plt.gca().transAxes,
             fontsize=9, verticalalignment='top',
             bbox=dict(boxstyle="round", facecolor="white", alpha=0.6))
    plt.grid(True)
    plt.tight_layout()
    output_path = Path(filename).stem + "_" + gwas_key + "_qq.png"
    output_paths += output_path + "\n"
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)

output_paths = output_paths[:-1]
print(output_paths, end="")
