import pandas as pd
import warnings
from pathlib import Path 
import sys
warnings.filterwarnings("ignore")

filename = sys.argv[1]
gwas = pd.read_csv(filename)
num_phenotypes = gwas['phenotype'].nunique()

gwas = gwas[gwas["#chr"] != "intercept"]

sig_threshold = 0.05/(len(gwas)/num_phenotypes) # bonferonni correction

gwas = gwas[gwas["pvalue"] < sig_threshold]

output_path = Path(filename)
gwas.to_csv(output_path, index=False)
print("File modified: SNPs with insignificant associations removed from " + output_path.name + "", end="")
