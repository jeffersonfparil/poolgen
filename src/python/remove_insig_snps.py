import pandas as pd
import warnings
from pathlib import Path 
import sys
warnings.filterwarnings("ignore")


filename = sys.argv[1]
gwas = pd.read_csv(filename, index_col=0)
num_phenotypes = gwas['phenotype'].nunique()

gwas["chromosome"] = gwas.index
gwas = gwas[gwas["chromosome"] != "intercept"]
gwas["position"] = gwas["pos"]

sig_threshold = 0.05/(len(gwas)/num_phenotypes) # bonferonni correction

gwas = gwas[gwas["pvalue"] < sig_threshold]

output_path = Path(filename)
gwas.to_csv(output_path)
print("DONE: Insignificant SNPs removed from " + output_path.name + "", end="")
