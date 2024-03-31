#!/bin/bash

while [[ "$#" -gt 0 ]]; do
    case "$1" in
        -g|--gene_name) # This is the name of the gene that you want searched
          GENE_NAME="$2"
          shift
          shift
          ;;
        -loc|--gene_location) # This is the file with the gene name locations extended
          GENE_LOCATION="$2"
          shift
          shift
          ;;
        -v|--vcf) # This is the vcf file that needs to be searched
          VCF_FILE="$2"
          shift
          shift
          ;;
        -o|--output) # This is the directory where the output will be piped to
          OUTPUT_DIR="$2"
          shift
          shift
          ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

# Check if the gene_location file exists
if [[ ! -f "$GENE_LOCATION" ]]; then
    echo "The gene_location file does not exist. Please enter the correct path:"
    read GENE_LOCATION
fi

# Check if the vcf_file exists
if [[ ! -f "$VCF_FILE" ]]; then
    echo "The vcf file does not exist. Please enter the correct path:"
    read VCF_FILE
fi

# Check if the output_dir exists
if [[ ! -d "$OUTPUT_DIR" ]]; then
    echo "The output directory does not exist. Do you want to create it? (yes/no):"
    read RESPONSE
    if [[ "$RESPONSE" == "yes" ]]; then
        mkdir -p "$OUTPUT_DIR"
        echo "Output directory created at $OUTPUT_DIR"
    else
        echo "Exiting the program."
        exit 1
    fi
fi
