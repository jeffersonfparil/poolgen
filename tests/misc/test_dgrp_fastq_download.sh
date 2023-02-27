#!/bin/bash

### Download DRGP raw sequences (List of SRA experiment IDs from: https://www.hgsc.bcm.edu/arthropods/dgrp-lines)
cd poolgen/

echo '#!/bin/bash
i=$1
f=$2
dirbase=$3
echo "###############################################"
line=$(head -n${i} ${f} | tail -n1)
echo $line
name=$(echo "$line" | cut -d"," -f1)
urls=( $(echo "$line" | rev | cut -d"," -f1 | rev | sed "s/ and / /g") )
# printf "%s\n" ${urls[@]}
# echo "${name}---${urls}"
for u in "${urls[@]}"
do
    # u=${urls[0]}
    outdir=${dirbase}/${u}/
    prefetch $u --output-directory ${outdir}
    for srr in $(find ${outdir} -name "SRR*" | grep ".sra$")
    do
        # srr=$(find ${outdir} -name "SRR*" | grep ".sra$" | head -n1)
        srr_name=$(basename $srr)
        bam_name=${srr_name%.sra*}.bam
        sam-dump ${srr} | samtools view -bS - > ${dirbase}/${name}-${bam_name}
    done
    rm -R ${outdir}
done
' > download_dgrp_for_parallelisation.sh
chmod +x download_dgrp_for_parallelisation.sh

### Execute in parallel
f="tests/misc/test_dgrp_fastq_urls.csv"
d="tests/misc"
time \
parallel -j 32 \
./download_dgrp_for_parallelisation.sh \
    {} \
    ${f} \
    ${d} \
    ::: $(seq 2 $(cat $f | wc -l))

### Download phenotype data
wget -O DGRP_phenotype_lifetime_fecundity.tmp.xlsx  \
    http://dgrp2.gnets.ncsu.edu/data/website/Durham_DGRP_Line_Means_NComm.xlsx
ssconvert DGRP_phenotype_lifetime_fecundity.tmp.xlsx \
    DGRP_phenotype_lifetime_fecundity.tmp.csv
tail -n+3 DGRP_phenotype_lifetime_fecundity.tmp.csv | \
    cut -d',' -f1,13 > DGRP_phenotype_lifetime_fecundity.csv
rm *.tmp*

### Determin the pooling based on the phenotypes
R
args = commandArgs(trailingOnly=TRUE)
args = c("DGRP_phenotype_lifetime_fecundity.csv", "1", "2", "Lifetime_fecundity", "5", "0.2", "0.2", "0.2", "0.2", "0.2")
fname = args[1]
line_col = as.numeric(args[2])
phen_col = as.numeric(args[3])
phen_name = args[4]
n_pools = as.numeric(args[5])
pool_sizes = as.numeric(args[6:(n_pools+5)])

dat = read.csv(fname, header=FALSE)
dat = data.frame(Line=dat[,line_col], Pheno=dat[,phen_col])
colnames(dat) = c("Line", phen_name)
dat = dat[order(dat[,2], decreasing=FALSE), ]
n = nrow(dat)
pool_sizes = cumsum(floor(pool_sizes * n))
if (pool_sizes[n_pools] < n) {
    idx = ceiling(n_pools/2)
    add = n - pool_sizes[n_pools]
    pool_sizes[idx:n_pools] = pool_sizes[idx:n_pools] + add
}




### Merge bam files so that each pool correspond to trait-based groupings 
for name in $(ls ${d} | cut -d'-' -f1 | sort | uniq)
do
    samtools merge \
        ${name}.bam \
        ${name}-*.bam
done

