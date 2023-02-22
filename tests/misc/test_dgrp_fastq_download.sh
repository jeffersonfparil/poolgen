#!/bin/bash

### Download DRGP raw sequences

cd poolgen/
while read -r line
do
    # line=$(head -n4 "tests/misc/test_dgrp_fastq_urls.csv" | tail -n1)
    name=$(echo "$line" | cut -d',' -f1)
    urls=( $(echo "$line" | rev | cut -d',' -f1 | rev | sed 's/ and / /g') )
    # printf "%s\n" ${urls[@]}
    # echo "${name}---${urls}"
    for u in "${urls[@]}"
    do
        # u=${urls[0]}
        outdir=./tests/misc/${u}/
        prefetch $u --output-directory ${outdir}
        for srr in $(find ${outdir} -name 'SRR*' | grep ".sra$")
        do
            # srr=$(find ./${u}/ -name 'SRR*' | grep ".sra$" | head -n1)
            fastq-dump --outdir ${outdir} ${srr}
        done
    done
done < <(tail -n+2 "./tests/misc/test_dgrp_fastq_urls.csv")
