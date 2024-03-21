#!/bin/bash

# input

# 1. gene list file
# This will be a sample file that usually looks like this
# Chr02 6115430 6176639

# Function to process a single gene
process_gene() {
    local gene_id="$1"
    local chromosome="$2"
    local lower="$3"
    local upper="$4"
    local output_path="$5"
    local input_file="$6"

    # Get relevant info for gene of interest
    gene_info=$(grep -w "$gene_id" "$input_file")

    if [[ -z "$gene_info" ]]; then
        echo "Gene $gene_id not found in the input file."
        return
    fi

    gene_cord=$(echo "$gene_info" | awk '{print $1":"$2"-"$3}')

    # Extract VCF info using tabix
    tabix -h "$input_file" "$gene_cord" | grep -E "#|$gene_id" > "${output_path}/${gene_id}_temp.vcf"

    # Convert to single line file
    cat "${output_path}/vcf.header" > "${output_path}/temp_${gene_id}_diversity_output.txt"

    cat "${output_path}/${gene_id}_temp.vcf" | "${BASE_DIR}/resources/scripts/vcfEffOnePerLine.pl" | grep -E "#|$gene_id" | \
    java -jar "${BASE_DIR}/resources/scripts/SnpSift.jar" extractFields - CHROM POS "ANN[*].FEATUREID" REF ALT "ANN[*].EFFECT" "ANN[*].AA" "ANN[*].IMPACT" "GEN[*].GT" | \
    grep Chr >> "${output_path}/temp_${gene_id}_diversity_output.txt"

    # Remove multiple/non-relevant gene SNPs
    grep -v "-" "${output_path}/temp_${gene_id}_diversity_output.txt" > "${output_path}/individual_files_full/${gene_id}_diversity_output.txt"

    # Clean up
    rm "${output_path}/${gene_id}_temp.vcf"
    rm "${output_path}/temp_${gene_id}_diversity_output.txt"

    # Exchange allele scores for numeric values
    sed -i 's|0/0|0|g; s|0/1|1|g; s|1/1|2|g; s|./.|NA|g' "${output_path}/individual_files_full/${gene_id}_diversity_output.txt"

    # For filtering only on meaningful variants
    grep -E 'CHROM|LOW|MODERATE|HIGH' "${output_path}/individual_files_full/${gene_id}_diversity_output.txt" > "${output_path}/individual_files_impactful/${gene_id}_diversity_output_impactful.txt"
}

while [[ "$#" -gt 0 ]]; do
    case "$1" in
        -f|--gene_file) # This is the path to gene files such as Sbicolor_313_v3_1_gene_locations_extended.txt
          GENE_FILE="$2"
          shift
          shift
          ;;
        -g|--gene)
          GENE_ID="$2"
          shift
          shift
          ;;
        -c|--chromosome)
            if [[ -n "$2" ]]; then
                chromosome="$2"
                shift
            else
                echo "Error: --chromosome requires a value"
                exit 1
            fi
            ;;
        -l|--lower)
            if [[ -n "$2" ]]; then
                lower="$2"
                shift
            else
                echo "Error: --lower requires a value"
                exit 1
            fi
            ;;
        -u|--upper)
            if [[ -n "$2" ]]; then
                upper="$2"
                shift
            else
                echo "Error: --upper requires a value"
                exit 1
            fi
            ;;
        -o|--output)
            if [[ -n "$2" ]]; then
                output_path="$2"
                shift
            else
                echo "Error: --output requires a value"
                exit 1
            fi
            ;;
        -i|--input) # this is for the chromosome file
            if [[ -n "$2" ]]; then
                input_file="$2"
                shift
            else
                echo "Error: --input requires a value"
                exit 1
            fi
            ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

# Check if lower is smaller than upper
while [[ "$lower" -ge "$upper" ]]; do
    echo "Error: Lower value must be smaller than upper value"
    read -p "Enter the lower value: " lower
    read -p "Enter the upper value: " upper
done

# Check if the output path exists
while [[ ! -d "$output_path" ]]; do
    read -p "Output path does not exist. Do you want to create it? [y/n]: " create_output
    if [[ "$create_output" == "y" || "$create_output" == "Y" ]]; then
        if mkdir "$output_path" 2>/dev/null; then
            echo "Output path created: $output_path"
        else
            echo "Error: You do not have permission to create the output path."
            read -p "Please enter a valid output path: " output_path
        fi
    elif [[ "$create_output" == "n" || "$create_output" == "N" ]]; then
        read -p "Please enter a valid output path: " output_path
    else
        echo "Invalid input. Please enter 'y' or 'n'."
    fi
done

# Check if the input file exists
while [[ ! -f "$input_file" ]]; do
    read -p "Input file does not exist. Please enter a valid file path: " input_file
done

# Check if the gene file exists
if [[ ! -f "$GENE_FILE" ]]; then
  echo "Gene file not found: $GENE_FILE"
  echo "Please provide a valid gene file using the --gene_file or -f flag."
  exit 1
fi

# Process a single gene or multiple genes from a file
if [[ -n "$GENE_ID" ]]; then
    if grep -Fq "$GENE_ID" "$GENE_FILE"; then
        process_gene "$GENE_ID" "$chromosome" "$lower" "$upper" "$output_path" "$input_file"
    else
        echo "Gene $GENE_ID not found in the gene file."
    fi
else
    echo "Please provide either a single gene using the --gene or -g flag or a gene file using the --gene_file or -f flag."
    exit 1
fi
