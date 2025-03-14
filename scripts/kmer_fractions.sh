#!/bin/bash

KMC_PATH="${HOME}/kmc/bin/kmc"
if [ ! -x $KMC_PATH ]; then
   echo "Error: ${KMC_PATH} does not exist"
  exit 1

fi


# Initialize variables for arguments
input_fasta=""
output_name=""
min_k=""
max_k=""
sequence_number=""
step=""

# Parse command-line arguments
while getopts "i:o:k:K:s:" opt; do
  case ${opt} in
    i) input_fasta=${OPTARG} ;;
    o) output_name=${OPTARG} ;;
    k) min_k=${OPTARG} ;;
    K) max_k=${OPTARG} ;;
    s) step=${OPTARG} ;;
    \?)
      echo "Usage: $0 -i <input_fasta_file> -o <output_dir_name> -k <min_k_value> -K <max_k_value> -s <step>"
      exit 1
      ;;
  esac
done

# Check that all required arguments are provided
if [[ -z "$input_fasta" || -z "$output_name" || -z "$min_k" || -z "$max_k" ]]; then
  echo "Error: Missing required argument(s)"
  echo "Usage: $0 -i <input_fasta_file> -o <output_dir_name> -k <min_k_value> -K <max_k_value> -s <step>"
  exit 1
fi


if [ -d ${output_name} ]; then
  echo "Error: ${output_name} already exists. Exiting."
  exit 1
fi


mkdir $output_name
mkdir ${output_name}/kmc_tmp

# Print out the values (optional)
echo "Input FASTA file: $input_fasta"
echo "Output directory: $output_name"
echo "Minimum k value: $min_k"
echo "Maximum k value: $max_k"
echo "Step:" $step
# Your script logic here


KMC_DB=${output_name}/$(basename ${input_fasta%%.fa})_kmc
KMC_TMP=${output_name}/kmc_tmp 
KMC_OUT=${output_name}/kmc_out.txt
KMC_LOG=${output_name}/kmc.log

echo -e "kmer\trare\tunique\ttotal" >> $output_name/results.tsv

for i in $(seq $min_k $step $max_k); do
	echo $i;
	$KMC_PATH -ci2 -cs400000 -k$i -fm $input_fasta $KMC_DB $KMC_TMP > $KMC_OUT 2>> $KMC_LOG;
	RARE=$(grep 'below' $KMC_OUT | awk '{print $NF}');
	UNIQUE=$(grep 'unique k-mers' $KMC_OUT | awk '{print $NF}');
	ALL=$(grep 'Total no. of k-mers' $KMC_OUT | awk '{print $NF}');

	rare_frac=$(echo "scale=4; ${RARE} / ${ALL}" | bc);
	echo -e "$i\t$rare_frac\t$UNIQUE\t$ALL" >> ${output_name}/results.tsv;

done

rm $KMC_OUT
rm ${KMC_DB}.kmc_pre
rm ${KMC_DB}.kmc_suf
rmdir $KMC_TMP

echo "Done"
