#!/bin/bash

# Initialize variables for arguments
input_fasta=""
output_tsv=""
min_k=""
max_k=""
sequence_number=""
step=""

# Parse command-line arguments
while getopts "i:o:k:K:n:s:" opt; do
  case ${opt} in
    i) input_fasta=${OPTARG} ;;
    o) output_tsv=${OPTARG} ;;
    k) min_k=${OPTARG} ;;
    K) max_k=${OPTARG} ;;
    n) sequence_number=${OPTARG} ;;
    s) step=${OPTARG};;
    \?)
      echo "Usage: $0 -i <input_fasta_file> -o <output_tsv_file> -k <min_k_value> -K <max_k_value> -n <sequence_number> -s <step>"
      exit 1
      ;;
  esac
done

# Check that all required arguments are provided
if [[ -z "$input_fasta" || -z "$output_tsv" || -z "$min_k" || -z "$max_k" || -z "$sequence_number" ]]; then
  echo "Error: Missing required argument(s)"
  echo "Usage: $0 -i <input_fasta_file> -o <output_tsv_file> -k <min_k_value> -K <max_k_value> -n <sequence_number> -s <step>"
  exit 1
fi

if [ -f ${output_tsv} ]; then
  echo "Error: ${output_tsv} already exists. Exiting."
  exit 1
fi

# Print out the values (optional)
echo "Input FASTA file: $input_fasta"
echo "Output TSV file: $output_tsv"
echo "Minimum k value: $min_k"
echo "Maximum k value: $max_k"
echo "Sequence number: $sequence_number"
echo "Step:" $step
# Your script logic here

KMC_PATH="${HOME}/kmc/bin/kmc"
TOOL_PATH="${HOME}/kmc/bin/kmc_tools"
DUMP_PATH="${HOME}/kmc/bin/kmc_dump"


filter=$(($sequence_number + 1))
KMC_DB=$(basename ${input_fasta%%.fa})_kmc



echo -e "kmer\trare\tover-represented\ttotal" >> $output_tsv

for i in $(seq $min_k $step $max_k); do
	echo $i;
	$KMC_PATH -ci2 -cs400000 -k$i -fm $input_fasta $KMC_DB kmc_tmp > kmc_out.txt 2>> kmc.log;
	RARE=$(grep 'below' kmc_out.txt | awk '{print $NF}');
	UNIQUE=$(grep 'unique k-mers' kmc_out.txt | awk '{print $NF}');
	ALL=$(grep 'Total no. of k-mers' kmc_out.txt | awk '{print $NF}');

	$DUMP_PATH -ci${filter} $KMC_DB over.txt;
	OVER=$(awk '{sum += $NF} END {print sum}' over.txt);

	rare_frac=$(echo "scale=4; ${RARE} / ${ALL}" | bc);
	over_frac=$(echo "scale=4; ${OVER} / ${ALL}" | bc);

	echo -e "$i\t$rare_frac\t$over_frac\t\t$ALL" >> $output_tsv;

done

echo "Done"
