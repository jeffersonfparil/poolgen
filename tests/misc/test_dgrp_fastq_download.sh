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


### Merge bam files so that each pool correspond to trait-based groupings 
for name in $(ls ${d} | cut -d'-' -f1 | sort | uniq)
do
    samtools merge \
        ${name}.bam \
        ${name}-*.bam
done

