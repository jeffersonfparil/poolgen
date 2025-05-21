import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import warnings
from pathlib import Path 
warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser(prog="QQ plotter")
parser.add_argument("filename")
args = parser.parse_args()

gwas = pd.read_csv(args.filename, index_col=0)
gwas = gwas.reset_index(drop=True)

pvalues = gwas["pvalue"].sort_values().values
n = len(pvalues)
expected = -np.log10(np.linspace(1/n, 1, n))
observed = -np.log10(pvalues)

df_qq = pd.DataFrame({"expected": expected, "observed": observed})

plt.figure(figsize=(6, 6))
sns.scatterplot(data=df_qq, x="expected", y="observed", alpha=0.5, edgecolor=None)
plt.plot([0, max(expected)], [0, max(expected)], color="red", linestyle="--")
plt.xlabel("Expected $-\\log_{10}p$")
plt.ylabel("Observed $-\\log_{10}p$")
plt.title("Q-Q Plot of GWAS p-values uniformity (" + Path(args.filename).stem + ")")
plt.grid(True)
plt.tight_layout()
output_path = Path(args.filename).stem + "_qq.png"
plt.savefig(output_path, dpi=300)
print(output_path, end = "")
