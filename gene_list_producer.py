import os
import pandas as pd
import argparse

# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("-l", "--gene_location_file", help="Gene location file")
parser.add_argument("-start", "--start_point", type=int, help="Start point")
parser.add_argument("-stop", "--stop_point", type=int, help="Stop point")
parser.add_argument("-c", "--chromosome", help="Chromosome")
parser.add_argument("-o", "--output_file", help="Output file (including path)")
parser.add_argument("-b", "--bulk_input_file", help="Bulk input file as a tsv.")

args = parser.parse_args()

def region_file_parser(file_name: str, chromosome: str, start: int, stop: int) -> pd.DataFrame:
    gene_location = pd.read_csv(file_name, sep="\t")
    region_gene_subset = gene_location[
        (gene_location['Chrom'] == chromosome) &
        (gene_location['Ext_Start'] >= start) &
        (gene_location['Ext_Start'] <= stop)
    ].drop_duplicates(subset='Gene')
    return region_gene_subset

def write_output_files(gene_subsets: pd.DataFrame, output_file: str):
    output_dir = os.path.dirname(output_file)
    if not output_dir:
        output_dir = os.getcwd()  # Set to current working directory if output_dir is empty
    file_base_name = os.path.splitext(os.path.basename(output_file))[0]

    # Concatenate all subsets and drop duplicates
    final_gene_subset = pd.concat(gene_subsets).drop_duplicates(subset='Gene')

    # Write the final DataFrame to a CSV file
    gene_list_file_path = f"{output_dir}/{file_base_name}_gene_table.csv"
    final_gene_subset.to_csv(gene_list_file_path, sep="\t", index=False)

    # Write the gene list to a text file
    gene_table_file_path = f"{output_dir}/{file_base_name}_gene_list.txt"
    with open(gene_table_file_path, 'w') as f:
        for gene in final_gene_subset['Gene'].unique():
            f.write(f"{gene}\n")

if args.bulk_input_file is None:
    subset = region_file_parser(args.gene_location_file, args.chromosome, args.start_point, args.stop_point)
    write_output_files([subset], args.output_file)
else:
    bulk_file = pd.read_csv(args.bulk_input_file, sep="\t")
    subsets = []
    for i in bulk_file.itertuples():
        gene_file = getattr(i, 'genloc', args.gene_location_file)
        output_file = getattr(i, 'output_file', args.output_file)
        subset = region_file_parser(gene_file, i.Chrom, i.start, i.stop)
        subsets.append(subset)
    write_output_files(subsets, output_file)
