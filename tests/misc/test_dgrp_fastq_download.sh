#!/bin/bash

### Download DRGP raw sequences (List of SRA experiment IDs from: https://www.hgsc.bcm.edu/arthropods/dgrp-lines)
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

### Remove empty bam files
for f in $(ls *.bam)
do
    n=$(head $f | wc -l)
    if [ $n -lt 10 ]
    then
        rm $f
    fi
done

### Index bam files
parallel samtools index {} ::: $(ls *bam)

### Remove poorly indexed files
for f in $(ls *.bam.bai)
do
    # f=$(ls *.bam.bai | head -n1 | tail -n1)
    n=$(head $f | wc -l)
    if [ $n -lt 10 ]
    then
        rm ${f%.bam*}*
    fi
done

### Download phenotype data
wget -O DGRP_phenotype_lifetime_fecundity.tmp.xlsx  \
    http://dgrp2.gnets.ncsu.edu/data/website/Durham_DGRP_Line_Means_NComm.xlsx
ssconvert DGRP_phenotype_lifetime_fecundity.tmp.xlsx \
    DGRP_phenotype_lifetime_fecundity.tmp.csv
tail -n+3 DGRP_phenotype_lifetime_fecundity.tmp.csv | \
    cut -d',' -f1,13 > DGRP_phenotype_lifetime_fecundity.csv
rm *.tmp*

### Download the GWAS results using individual genotype data
wget https://static-content.springer.com/esm/art%3A10.1038%2Fncomms5338/MediaObjects/41467_2014_BFncomms5338_MOESM1247_ESM.xlsx

### Determine the pooling based on the phenotypes
Rscript test_dgrp_pooling.r DGRP_phenotype_lifetime_fecundity.csv \
                            1 \
                            2 \
                            Lifetime_fecundity \
                            5 \
                            0.2 0.2 0.2 0.2 0.2

### Fix name prefix from RAL to DGRP
for f in $(ls *_Pool_*_bam_list.txt)
do
    sed -i 's/RAL/DGRP/g' $f
done

### Merge bam files so that each pool correspond to trait-based groupings 
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
n=100000
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


