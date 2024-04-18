#!/bin/bash

# check if the dependencies are present and available on path
command -v tabix >/dev/null || { echo "This shell does not have access to tabix. This software requires tabix in path"; exit; }
command -v vcftools >/dev/null || { echo "This shell does not have access to vcftools. This software requires vcftools in path"; exit; }

# defining functions

gene_region_finder () {
  # re-assign the variables here to use inside the function
  chromosome="$1"
  vcf_file_holder="$2"
  loci="$3"
  output="$4"

  # Check if the output directory exists
  if [[ ! -d "$output" ]]; then
    read -p "Output directory '$output' does not exist. Would you like to create it? (y/n) " yn
    case $yn in
        [Yy]* ) mkdir -p "$output";;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no."; exit;;
    esac
  fi

  base_name="$(basename -- $vcf_file)"
  base_name="${base_name%.*}" # this should be the base of the file without the extension

  # Check if the .tbi file exists
  if [[ ! -f "${vcf_file}.tbi" ]]; then
    read -p "The .tbi file for the given vcf file does not exist. Would you like to generate it? (y/n) " yn
    case $yn in
        [Yy]* ) tabix -p vcf $vcf_file;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no."; exit;;
    esac
  fi

  # Crunching the numbers for the linkage loci range from the input
  snp_start_ld=$((loci-500000))
  snp_stop_ld=$((loci+500000))

  ## make sure that the start LD is not below zero

  if (( ${snp_start_ld} < 1)); then
    snp_start_ld=0
    echo "the given loci was too close to the start and needed to be reset to 0"
  fi

  # generate a base name for the output file using the base name of the input file and the ld range
  tabix_output_file_name="${base_name}_${snp_start_ld}_${snp_stop_ld}"

  tabix_output_file="${output}/${tabix_output_file_name}.txt"

  tabix -h ${vcf_file} ${chromosome}:${snp_start_ld}-${snp_stop_ld} > ${tabix_output_file}

  cat <<- HERE > "${chromosome}_${loci}_temp_region_file.txt"
  Chr Snp
  $chromosome $loci
HERE

  vcftools --vcf $tabix_output_file \
    --geno-r2-positions "${chromosome}_${loci}_temp_region_file.txt" \
    --ld-window 500 \
    --out $output/${tabix_output_file_name}
  
  # remove unwanted files here
  rm "${chromosome}_${loci}_temp_region_file.txt" # the temp file that holds the `--geno-r2-positions` input data for vcftools, was created using the HERE doc.
  
  rm $output/${tabix_output_file_name}.log # the log file that is completely unnecessary for us.

  rm $tabix_output_file

  vcf_file="$output/${tabix_output_file_name}.list.geno.ld"
  
  awk_file_output="$output/${base_name}_${snp_start_ld}_${snp_stop_ld}.ld_data"

  awk '$6+0 >= 0.1' $vcf_file | grep -v "nan" | grep -v "N_INDV" | sort -k4,4n > $awk_file_output

  # Print the first item of the 4th column
  start=$(awk '{if(NR==1) print $4}' $awk_file_output)
  
  # Print the last item of the 4th column
  stop=$(awk '{last=$4} END{print last}' $awk_file_output)

  printf "%s %s %s\n" "$chromosome" "$start" "$stop" > $awk_file_output
  
  rm $vcf_file # remove unnecessary intermediate file
  
 }


# CLI flags
while test $# -gt 0; do
  case "$1" in
    --file)
      shift
      input_file="$1"
      shift
      ;;
    --output)
      shift
      output="${1:-"panvar_run"}"
      shift
      ;;
    *)
      # Handle other flags or break out of the loop
      break
      ;;
  esac
done

if [[ -n "$input_file" ]]; then
  if [[ ! -f "$input_file" ]]; then
    echo "Input file '$input_file' not found."
    exit 1
  fi

  # Process the TSV file
  while IFS=$' ' read -r chromosome vcf_file loci; do
    # Use the parsed values here
    #DEBUG echo "first debug Processing: $chromosome $vcf_file $loci"
    #DEBUG echo "$output"
    # pass arguments to the function that wrangles the data
    gene_region_finder "$chromosome" "$vcf_file" "$loci" "$output" 
  done < "$input_file"
else
  # Reset the argument pointer if necessary
  # set -- 

  while test $# -gt 0; do
    case "$1" in
      --chromosome)
        shift
        chromosome=$1
        shift
        ;;
      --vcf_file)
        shift
        vcf_file=$1
        shift
        ;;
      --output)
        shift
        output=$1
        shift
        ;;
      --loci)
        shift
        loci=$1
        shift
        ;;
      *)
        # Handle unknown option or break
        break
        ;;
    esac
  done
  # Validate and use the arguments here
  gene_region_finder "$chromosome" "$vcf_file" "$loci" "$output"
fi
