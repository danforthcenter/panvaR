import os
import pandas as pd
import argparse
import glob

# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("--gene_location_file", help="Gene location file")
parser.add_argument("--start_point", type=int, help="Start point")
parser.add_argument("--stop_point", type=int, help="Stop point")
parser.add_argument("--chromosome", help="Chromosome")
parser.add_argument("--output_dir", default=os.getcwd(), help="Output directory (default is current directory)")
parser.add_argument("--bulk_input_file", help="If you want to supply a bulk set of inputs as a tsv file. Supply the same thing as single input but each arg should be a named column and each row is then a set.")

args = parser.parse_args() # collect the arguments

def region_file_parser(file_name: str, chromosome: str, start: int, stop: int):

    # Load gene location data
    gene_location = pd.read_csv(file_name, sep="\t")
    
    file_base_name = os.path.splitext(os.path.basename(file_name))[0]

    # This is only filtering the region file for rows and SNPs of interest.
    region_gene_subset = gene_location[(gene_location['Chrom'] == chromosome) & (gene_location['Ext_Start'] >= start) & (gene_location['Ext_Start'] <= stop)]

    if (len(region_gene_subset) <= 2):

        new_start = start + 1000
        new_stop = stop  + 1000

        if (new_start > min(gene_location['Ext_Start']) & new_stop < max(gene_location['Ext_Stop'])):

            region_file_parser(file_name, chromosome, new_start, new_stop)
        else:
            print(f"For file {file_name} No genes found in the region {start} and {stop} -ran out of regions to explore.")

    output_file_name = f"{args.output_dir}/{file_base_name}_{chromosome}_{start}_{stop}.gene_list"

    region_gene_subset.to_csv(output_file_name,sep="\t",header=True,index=False)



if args.bulk_input_file is None:

    region_file_parser(file_name = args.gene_location_file, chromosome = args.chromosome, start=args.start_point,stop=args.stop_point)

elif args.bulk_input_file is not None:

    bulk_file = pd.read_csv(args.bulk_input_file,sep="\t")

    for i in bulk_file.itertuples():
        region_file_parser(file_name=i.File, chromosome=i.Chromosome,start=i.Ext_Start,stop=i.Ext_Stop)
