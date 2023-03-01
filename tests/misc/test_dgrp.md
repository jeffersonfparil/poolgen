# Using DGRP reference panel 2 data

1. Download DRGP raw sequences (List of SRA experiment IDs from: https://www.hgsc.bcm.edu/arthropods/dgrp-lines)
```shell
cd poolgen/tests/misc

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
f="test_dgrp_fastq_urls.csv"
d="."
time \
parallel -j 32 \
./download_dgrp_for_parallelisation.sh \
    {} \
    ${f} \
    ${d} \
    ::: $(seq 2 $(cat $f | wc -l))
```

2. Filter bam files

- Remove empty bam files
```shell
for f in $(ls *.bam)
do
    n=$(head $f | wc -l)
    if [ $n -lt 10 ]
    then
        rm $f
    fi
done
```

- Index bam files
```shell
parallel samtools index {} ::: $(ls *bam)
```

- Remove poorly indexed bam files
```shell
for f in $(ls *.bam.bai)
do
    # f=$(ls *.bam.bai | head -n1 | tail -n1)
    n=$(head $f | wc -l)
    if [ $n -lt 10 ]
    then
        rm ${f%.bam*}*
    fi
done
```

3. Download phenotype data
```shell
wget -O DGRP_phenotype_lifetime_fecundity.tmp.xlsx  \
    http://dgrp2.gnets.ncsu.edu/data/website/Durham_DGRP_Line_Means_NComm.xlsx
ssconvert DGRP_phenotype_lifetime_fecundity.tmp.xlsx \
    DGRP_phenotype_lifetime_fecundity.tmp.csv
tail -n+3 DGRP_phenotype_lifetime_fecundity.tmp.csv | \
    cut -d',' -f1,13 > DGRP_phenotype_lifetime_fecundity.csv
rm *.tmp*
```

4. Download the GWAS results which used individual genotype data
```shell
wget https://static-content.springer.com/esm/art%3A10.1038%2Fncomms5338/MediaObjects/41467_2014_BFncomms5338_MOESM1247_ESM.xlsx
```

5. Determine the pooling based on the phenotypes
```shell
Rscript test_dgrp_pooling.r DGRP_phenotype_lifetime_fecundity.csv \
                            1 \
                            2 \
                            Lifetime_fecundity \
                            5 \
                            0.2 0.2 0.2 0.2 0.2
```

6. Fix name prefix from RAL to DGRP
```shell
for f in $(ls *_Pool_*_bam_list.txt)
do
    sed -i 's/RAL/DGRP/g' $f
done
```

7. Subsample and merge bam files so that each pool correspond to trait-based groupings 
```shell
echo '#!/bin/bash
f=$1
i=$2
n=$3
# f=$(ls *_Pool_*_bam_list.txt | head -n3 | tail -n1)
# i=2
# n=100000
wild_bam=$(head -n${i} $f | tail -n1)
bam=$(ls ${wild_bam} | head -n1) ### use only 1 bam per line
echo "##########################"
echo $bam
frac=$( samtools idxstats ${bam} | cut -f3 | awk -v n=$n @BEGIN {total=0} {total += $1} END {frac=n/total; if (frac > 1) {print 1} else {print frac}}@ )
samtools view -bs ${frac} ${bam} > ${bam%.bam*}-subsample.bam
' | sed "s/@/'/g" > subsample_parallel.sh
chmod +x subsample_parallel.sh
n=10000000
for f in $(ls *_Pool_*_bam_list.txt)
do
    # f=$(ls *_Pool_*_bam_list.txt | head -n3 | tail -n1)
    name=${f%_bam_*}
    ### Downsample
    parallel ./subsample_parallel.sh ${f} {} ${n} ::: $(seq 1 $(cat $f | wc -l))
    ### Merge
    samtools merge \
        ${name}.bam \
        $(cat $f | sed 's/*.bam/*-subsample.bam/g')\
    ### Cleanup
    rm $(cat $f | sed 's/*.bam/*-subsample.bam/g')
done
```

8. Download the Drosophila reference genome and annotation
```shell
wget -O -  ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r5.13_FB2008_10/fasta/dmel-all-chromosome-r5.13.fasta.gz | gunzip - > Drosophila_melanogaster.fna
wget -O - http://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r5.13_FB2008_10/gff/dmel-all-r5.13.gff.gz | gunzip - > Drosophila_melanogaster.gff
```

9. Pileup
```shell
ls DGRP_phenotype_lifetime_fecundity_Pool_*.bam > DGRP_phenotype_lifetime_fecundity_BAM_LIST.txt
samtools mpileup \
    --min-MQ 20 \
    --min-BQ 20 \
    --fasta-ref Drosophila_melanogaster.fna \
    --bam-list DGRP_phenotype_lifetime_fecundity_BAM_LIST.txt \
    --output DGRP_phenotype_lifetime_fecundity.pileup
```

10. Synchronise
```shell
cargo run -- pileup2sync \
    -f DGRP_phenotype_lifetime_fecundity.pileup \
    -o DGRP_phenotype_lifetime_fecundity.sync \
    --pool-names DGRP_phenotype_lifetime_fecundity_BAM_LIST.txt \
    --file-format sync \
    --n-threads 32 \
    --min-cov 10

```

11. Test Fisher's exact test
```shell
cargo run -- fisher_exact_test \
    -f DGRP_phenotype_lifetime_fecundity.sync \
    --n-threads 32
```