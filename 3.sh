#!/bin/bash

#!/bin/bash

process_gene(){

    gene_name="$1"
    gene_location="$2"
    vcf_file="$3"
    output_dir="$4"
    
    # Create a sub-directory for the gene
    gene_dir="${output_dir}/${gene_name}"
    mkdir -p "$gene_dir"

    # Check if gene_name exists in gene_location file
    gene_line=$(grep -F "$gene_name" "$gene_location")
    if [[ -z "$gene_line" ]]; then
        echo "Warning: $gene_name does not exist in $gene_location"
        exit 1
    fi

    # If gene_name exists in gene_location file, use awk to create the required string
    gene_string=$(echo "$gene_line" | awk '{print $1 ":" $2 "-" $3}') # this creates the query for tabix to use agains the vcf file

    tabix -h $vcf_file $gene_string | \
        grep -E "#|${gene_name}" > "${gene_dir}/${gene_name}_region_filtered.vcf"

    awk '{if(!/^##/) {print $0; exit}}' $vcf_file > "${gene_dir}/${gene_name}.vcf_header"

    # New line added here
    cat  ${gene_dir}/${gene_name}_region_filtered.vcf | vcfEffOnePerLine.pl | grep -E "#|${gene_name}" | \
    java -jar SnpSift.jar extractFields - CHROM POS "ANN[*].FEATUREID" REF ALT "ANN[*].EFFECT" "ANN[*].AA" "ANN[*].IMPACT" "GEN[*].GT" | \
    grep Chr >> "${gene_dir}/${gene_name}_processed.vcf"

    # Rijan: remove multiple/non-relevant gene snps
    grep -v "-" "${gene_dir}/${gene_name}_processed.vcf" > "${gene_dir}/${gene_name}_single_processed.vcf"

    sed 's|0\/0|0|g; s|0\/1|1|g; s|1\/1|2|g; s|.\/.|NA|g' "${gene_dir}/${gene_name}_single_processed.vcf" > "${gene_dir}/${gene_name}_numerical_allele_scores.vcf"

    grep -E 'CHROM|LOW|MODERATE|HIGH' "${gene_dir}/${gene_name}_numerical_allele_scores.vcf" > "${gene_dir}/${gene_name}_impactful.vcf"
}



while [[ "$#" -gt 0 ]]; do
    case "$1" in
        -g|--gene_name) # This is the name of the gene that you want searched
          shift
          GENE_NAME="$1"
          shift
          ;;
        -loc|--gene_location) # This is the file with the gene name locations extended
          shift
          GENE_LOCATION="$1"
          shift
          ;;
        -v|--vcf) # This is the vcf file that needs to be searched
          shift
          VCF_FILE="$1"
          shift
          ;;
        -o|--output) # This is the directory where the output will be piped to
          shift
          OUTPUT_DIR="$1"
          shift
          ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
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
