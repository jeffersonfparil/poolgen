import pandas as pd
import warnings
from pathlib import Path 
import sys
warnings.filterwarnings("ignore")

gwas_filename = sys.argv[1]
gwas = pd.read_csv(gwas_filename, index_col=0)
num_phenotypes = gwas['phenotype'].nunique()

gwas["chromosome"] = gwas.index
gwas = gwas[gwas["chromosome"] != "intercept"]

sig_threshold = 0.05/(len(gwas)/num_phenotypes) # bonferonni correction

gwas = gwas[gwas["pvalue"] < sig_threshold]

gff_filename = sys.argv[2]
gff_window_size = int(sys.argv[3])
gff = pd.read_csv(gff_filename, sep="\t", header=None, comment='#')

gff.columns = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]

gwas['key'] = 1
gff['key'] = 1
gwas['pos'] = gwas['pos'].astype(int)
merged = pd.merge(gwas, gff, on='key').drop(columns='key')

filtered = merged[
    (merged['start'] <= merged['pos'] + gff_window_size) & 
    (merged['end'] >= merged['pos'] - gff_window_size)
]

gff = filtered.drop(columns='pos').drop_duplicates().reset_index(drop=True)

output_path = Path(gff_filename)
output_path = output_path.with_name(output_path.stem + "_filtered" + output_path.suffix)
gff.to_csv(output_path, sep="\t", header=False, index=False, quoting=3)

print("File created: " + output_path.name, end="")
