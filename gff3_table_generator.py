import gffutils
import argparse
import os

def main():
    parser = argparse.ArgumentParser(description='Extract gene features from a GFF3 file')
    parser.add_argument('gff3_file', help='Path to the GFF3 file')
    parser.add_argument('-o', '--output', help='Output file name (default: print to console)')
    parser.add_argument('-d', '--db_file', default='default.db', help='Path to the output database file (default: default.db)')
    args = parser.parse_args()

    # Create a GFF database from the GFF3 file
    db = gffutils.create_db(args.gff3_file, args.db_file)

    # Extract gene features (locations)
    genes = db.features_of_type('gene')

    # Write output to file or print to console
    if args.output:
        with open(args.output, 'w+') as f:
            f.write("Chrom\tExt_Start\tExt_Stop\tGene\n")  # header
            for gene in genes:
                f.write(f"{gene.chrom}\t{gene.start}\t{gene.end}\t{gene.attributes['Name'][0]}\n")
    else:
        print("Chrom\tExt_Start\tExt_Stop\tGene")  # header
        for gene in genes:
            print(f"{gene.chrom}\t{gene.start}\t{gene.end}\t{gene.attributes['Name'][0]}")

    # Delete the database file if it's the default one
    if args.db_file == 'default.db':
        os.remove(args.db_file)

if __name__ == '__main__':
    main()
