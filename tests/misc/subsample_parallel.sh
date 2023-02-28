#!/bin/bash
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
frac=$( samtools idxstats ${bam} | cut -f3 | awk -v n=$n 'BEGIN {total=0} {total += $1} END {frac=n/total; if (frac > 1) {print 1} else {print frac}}' )
samtools view -bs ${frac} ${bam} > ${bam%.bam*}-subsample.bam

