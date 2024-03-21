#!/bin/bash

# check if the dependencies are present and available on path
command -v tabix >/dev/null || { echo "This shell does not have access to tabix. This software requires tabix in path"; exit; }
command -v vcftools >/dev/null || { echo "This shell does not have access to vcftools. This software requires vcftools in path"; exit; }

# defining functions

gene_region_finder () {
  # re-assign the variables here to use inside the function
  chromosome="$1"
  vcf_file_holder="$2"
  distance="$3"
  output="$4"

  base_name="$(basename -- $vcf_file)"
  base_name="${base_name%.*}" # this should be the base of the file without the extension

  # Crunching the numbers for the linkage distance range from the input
  snp_start_ld=$((distance-500000))
  snp_stop_ld=$((distance+500000))

  ## make sure that the start LD is not below zero

  if (( ${snp_start_ld} < 1)); then
    snp_start_ld=0
    echo "the given distance was too close to the start and needed to be reset to 0"
  fi

  # generate a base name for the output file using the base name of the input file and the ld range
  tabix_output_file_name="${base_name}_${snp_start_ld}_${snp_stop_ld}"

  tabix_output_file="${output}/${tabix_output_file_name}.txt"

  tabix -h ${vcf_file} ${chromosome}:${snp_start_ld}-${snp_stop_ld} > ${tabix_output_file}

  cat <<- HERE > temp_region_file.txt
  Chr Snp
  $chromosome $distance
HERE

  vcftools --vcf $tabix_output_file \
    --geno-r2-positions temp_region_file.txt  \
    --ld-window 500 \
    --out $output/${tabix_output_file_name}
  
  # remove unwanted files here
  rm temp_region_file.txt # the temp file that holds the `--geno-r2-positions` input data for vcftools, was created using the HERE doc.
  rm $output/${tabix_output_file_name}.log # the log file that is completely unnecessary for us.

  rm $tabix_output_file

  vcf_file="$output/${tabix_output_file_name}.list.geno.ld"
  
  awk_file_output="$output/${base_name}_${snp_start_ld}_${snp_stop_ld}.ld_data"

  awk '$6 >= "0.1"' $vcf_file | grep -v "nan" | grep -v "N_IDV" | sort -k4,4n > $awk_file_output

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
  while IFS=$' ' read -r chromosome vcf_file distance; do
    # Use the parsed values here
    #DEBUG echo "first debug Processing: $chromosome $vcf_file $distance"
    #DEBUG echo "$output"
    # pass arguments to the function that wrangles the data
    gene_region_finder "$chromosome" "$vcf_file" "$distance" "$output" 
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
      --distance)
        shift
        distance=$1
        shift
        ;;
      *)
        # Handle unknown option or break
        break
        ;;
    esac
  done
  # Validate and use the arguments here
  gene_region_finder "$chromosome" "$vcf_file" "$distance" "$output"
fi

