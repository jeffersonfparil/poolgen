import pandas as pd
import argparse
import warnings
from pathlib import Path 
warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser(prog="SNP siever")
parser.add_argument("filename")
args = parser.parse_args()

gwas = pd.read_csv(args.filename, index_col=0)
num_phenotypes = gwas['phenotype'].nunique()

gwas["chromosome"] = gwas.index
gwas = gwas[gwas["chromosome"] != "intercept"]
gwas["position"] = gwas["pos"]

sig_threshold = 0.05/(len(gwas)/num_phenotypes) # bonferonni correction

gwas = gwas[gwas["pvalue"] < sig_threshold]

output_path = Path(args.filename)
gwas.to_csv(output_path)
