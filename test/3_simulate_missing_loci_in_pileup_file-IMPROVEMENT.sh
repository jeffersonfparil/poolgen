#!/bin/bash 

echo "##############################################"
echo "### SIMULATE MISSING LOCI IN A PILEUP FILE ###"
echo "##############################################"

### Input and output specifications:
while [ ! -z "$1" ]; do
  case "$1" in
    --pileup-file|-f) shift
      PILEUP=$1;;
    --fraction-missing-pools|-p) shift
      FRAC_MISSING_POOLS=$1;;
    --fraction-missing-loci|-l) shift
      FRAC_MISSING_LOCI=$1;;
    --read-length|-r) shift
      READ_LENGTH=$1;;
    --output|-o) shift
        OUTPUT=$1;;
    --help|-h) shift
      echo " --pileup-file|-f               [STR] pileup file."  
      echo " --fraction-missing-pools|-p    [FLOAT] fraction of the pools with missing loci."
      echo " --fraction-missing-loci|-l     [FLOAT] fraction of the loci with missing data."
      echo " --output|-o                    [STR; DEFAULT='out_simissing.pileup'] filename of the output pileup file ."
      echo " --help|-h                      help documentation."
      echo ""
      exit 0
      ;;
    *)
      echo "What is this?! $1"
      exit 1
      ;;
  esac
shift
done

#############################
### TEST
# PILEUP=test.pileup
# FRAC_MISSING_POOLS=0.10
# FRAC_MISSING_LOCI=0.10
# READ_LENGTH=25
#############################

### Just a fun pure bash progress bar
function fun_progress_bar {
  n_int_current=$1
  n_int_max=$2
  n_int_length=$3
  ########################
  ### TEST
  # n_int_current=23
  # n_int_max=100
  # n_int_length=30
  ########################
  let "n_factor = $n_int_max / $n_int_length"
  if [ $n_factor -gt 0 ]
  then
    let "n_int_current = $n_int_current / $n_factor"
    let "n_int_pending = $n_int_length - $n_int_current"
  else
    let "n_factor = $n_int_length / $n_int_max"
    let "n_int_current = $n_int_current * $n_factor"
    let "n_int_pending = $n_int_length - $n_int_current"
  fi
  let "n_int_percent_done = n_int_current * 100 / $n_int_length"
  if [ $n_int_percent_done -gt 100 ]
  then
    n_int_current=$n_int_length
    n_int_pending=0
    n_int_percent_done=100
  fi
  vec_str_progress=$(printf '=%.0s' $(seq 1 $n_int_current))
  vec_str_pending=$(printf ' %.0s' $(seq 1 $n_int_pending))
  printf "\r[$vec_str_progress$vec_str_pending] $n_int_percent_done%%"
}

### If no output filename was set use the default output filename.
if [ -z $OUTPUT]
then 
    OUTPUT="out_simissing.pileup"
fi

### Initialise output file
cp $PILEUP "$OUTPUT"

### Extract input file statistics
N_LOCI=$(cat ${PILEUP} | wc -l)
N_COLUMNS=$(head -n1 ${PILEUP} | awk '{ print NF }')
N_POOLS=$(echo "($N_COLUMNS / 3) - 1" | bc)

### Define the number of missing pools and loci to simulate
N_MISSING_LOCI=$(printf "%.0f\n" $(echo "$N_LOCI * $FRAC_MISSING_LOCI" | bc))
N_MISSING_POOLS=$(printf "%.0f\n" $(echo "$N_POOLS * $FRAC_MISSING_POOLS" | bc))

### Define the clusters of reads and the number of clusters we will be sampling
N_CHUNKS=$(echo "($N_LOCI + $READ_LENGTH - 1) / $READ_LENGTH" | bc)
N_MISSING_CHUNKS=$(echo "($N_MISSING_LOCI + $READ_LENGTH - 1) / $READ_LENGTH" | bc)

### Random sampling of missing loci
FNAME_LINES_OF_MISSING_LOCI="missing_loci.temp"
# seq 1 $N_LOCI | shuf -n $N_MISSING_LOCI | sort -V > "$FNAME_LINES_OF_MISSING_LOCI"
VEC_CHUNK_END_POSITIONS=$(seq 1 $N_CHUNKS | shuf -n $N_MISSING_CHUNKS | sort -V)
for chunk in $VEC_CHUNK_END_POSITIONS
do
  N_START=$(echo "(($chunk - 1) * $READ_LENGTH) + 1" | bc)
  if [ $chunk -eq $N_CHUNKS ]
  then
    N_END=$N_LOCI ### if we reach the last chunk then the last position should be limited to the total number of loci to account for mod(N_LOCI, N_CHUNKS) != 0
  else
    N_END=$(echo "$chunk * $READ_LENGTH" | bc)
  fi
  seq $N_START $N_END >> "$FNAME_LINES_OF_MISSING_LOCI"
done

### Random sampling of missing pools
FNAME_COLUMNS_OF_MISSING_POOLS="missing_pools.temp"
seq 1 $N_POOLS | shuf -n $N_MISSING_POOLS | sort -V > "$FNAME_COLUMNS_OF_MISSING_POOLS"
for line in $(cat "$FNAME_COLUMNS_OF_MISSING_POOLS")
do
    TRUE_COLUMN_NUMBER=$(echo "(($line - 1) * 3) + 1 + 3" | bc)
    sed -i "s/^${line}$/${TRUE_COLUMN_NUMBER}/g" "$FNAME_COLUMNS_OF_MISSING_POOLS"
done

### Modify the output file to include the missing loci
i=0
for locus in $(cat "$FNAME_LINES_OF_MISSING_LOCI")
do
    for pool_depth in $(cat "$FNAME_COLUMNS_OF_MISSING_POOLS")
    do
        pool_SNP=$(echo "$pool_depth + 1" | bc)
        pool_quality=$(echo "$pool_depth + 2" | bc)
        ### make sure we output with tab delimiters with -v OFS='\t'
        gawk -v FS='\t' -v OFS='\t' -i inplace "FNR==${locus}{\$${pool_depth}=0;\$${pool_SNP}=\$${pool_quality}=\"*\"};1" $OUTPUT
    done
    let "i = $i + 1"
    fun_progress_bar $i $N_MISSING_LOCI 40    
done

### Clean-up
rm *.temp

### Closing message
echo "##############################################"
echo "### Please find the output: '${OUTPUT}'"
echo "##############################################"
