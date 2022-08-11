#!/bin/bash
BAM_LIST=$1
REFERENCE_GENOME=$2
MAPQ=$3
BASQ=$4
OUTPUT=$5
####################################################
### TEST:
BAM_LIST=/data/Lolium/Quantitative_Genetics/04_MPILEUP/ALL_MERGED_ALL_bam.list
ls /data/Lolium/Quantitative_Genetics/03_BAM/*.bam > ${BAM_LIST}
REFERENCE_GENOME=/data/Lolium/Genomics/SEQUENCES/DNA/REFERENCE_GENOMES/Reference.fasta ### lope_V1.0.fasta fixed to a have single line per scaffold
MAPQ=42
BASQ=42
OUTPUT=/data/Lolium/Quantitative_Genetics/04_MPILEUP/ALL_MERGED_ALL.bam
####################################################
# time \
# samtools mpileup \
#     -aa \
#     --min-MQ ${MAPQ} \
#     --min-BQ ${BASQ} \
#     --fasta-ref ${REFERENCE_GENOME} \
#     --bam-list ${BAM_LIST} > ${OUTPUT} ### absolutely all positions
time \
samtools mpileup \
    -b ${BAM_LIST} \
    -d 100000 \
    -q ${MAPQ} \
    -Q ${BASQ} \
    -f ${REFERENCE_GENOME} \
    -o ${OUTPUT}

### MISC:
# samtools view ACC01.bam | awk '{sum+=$5} END { print "Mean MAPQ =",sum/NR}'
# samtools view ACC01.bam | head -n1 | cut -f11 | awk -l ordchr -v RS='.{1}' '{print ord(RT)+33}'
