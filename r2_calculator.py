import polars as pl
import pysam
import pandas
import numpy as np
from itertools import compress

def add_if_not_none(item):
    '''
    Input:
        A sub-array from a numpy array
    Does:
        If all the items in the received array are Integers - it adds them and returns them. Else if something in the array is None, it returns None
    Returns:
        Integer or None
    '''
    # If all the items are integers add them, else return None
    if all(item != None):
        return sum(item)
    else:
        return 


def fetch_alleles(vcf, chrom, loci):
    """
    Fetches alleles from a VCF file for a given chromosome and locus.

    Args:
        vcf: A VCF file object
        chrom: The chromosome to fetch alleles for
        loci: The locus to fetch alleles for

    Returns:
        A numpy array of alleles

    Raises:
        ValueError: If the locus does not have any variation
    """
    
    # Check if the locus of interest exists
    if ((loci - 1) >= 0):
        # Fetch the records for the given chromosome and locus
        current_loci = vcf.fetch(chrom, loci - 1 , loci)
        
        # Convert the records to a list, as it's not a list by default
        current_loci = list(current_loci)
        
        # Check if the locus has any variation
        if current_loci[0].pos != loci:
            raise ValueError("The loci does not have any variation!")
    
    # Initialize an empty list to store the alleles
    alleles = []
    
    # Iterate over the records
    for rec in current_loci:
        # Extract the alleles from the record
        alleles.extend([s.alleles for s in rec.samples.values()])
    
    # Convert the list of alleles to a numpy array
    alleles = np.array(alleles)
    
    # Return the numpy array of alleles
    return alleles

def allele_to_ints(allele_array, ref):

    """
    Converts Biallellic Alleles in Char form to Int form

    Args:
        the array of alleles (bi-allelic)
        reference Allele

    Returns:
        Array in interger format
    """

    conditions = [allele_array == ref, (allele_array != ref) & (allele_array != None)]
    choices = [0, 1]

    int_array = np.select(conditions, choices, default=None)

    return int_array

def reference_allele(vcf_file, chrom, loci):
    """
    Args:
        A vcf_file opened using pysam
        Chrom:str :- The chromosome at the loci
        loci:str :- The loci of interest

    Does:
        Creates a list from the iterator that queries the vcf_file and gets the reference allele
    Returns:
        allele:str
    """

    iterator = list(vcf_file.fetch(chrom, loci - 1, loci))

    allele = iterator[0].ref

    return allele

def genotype_correlation(X:list, Y:list):
    """
    Args:
        Array of 0,1,2s at one loci of interest
        Array of 0,1,2s at the second loci of interest
    
    Does:
        Linkage disequilibrium calculation
    
    Returns:
        A float that has the r^2 value for the loci
        
    """
    cov_xy = np.cov(X, Y)[0, 1]
    
    # Calculate the variances of X and Y
    var_x = np.var(X)
    var_y = np.var(Y)
    
    # Calculate the genotype correlation (r^2)
    r_squared = (cov_xy ** 2) / (var_x * var_y)
    
    return r_squared


# use vcf_path to pass path to vcf_file, else just pass vcf_data
def rtwo_calculator(chromosome:str, loci:list, vcf_path:str = None, vcf_data = None):
    """
    This function calculates r^2 values between the first locus and multiple other loci.

    Args:
        chromosome (str): The chromosome in the file of interest.
        loci (list): A list of loci of interest. The first locus is used as the reference.
        vcf_path (str): The path to the VCF file of interest.
        vcf_data: An open VCF file object (optional).

    Returns:
        dict: A dictionary with locus pairs as keys and r^2 values as values.
    """
    # if vcf_path is supplied open the vcf_file
    if vcf_path != None:
        vcf_data = pysam.VariantFile(vcf_path)
    # else just use vcf_dat as supplied

    # get the phased genotypes in [<char>,<char>] format for the first locus
    genotype_array_ref = fetch_alleles(vcf_data, chromosome, loci[0])

    # convert the phased genotypes into numerical format for the first locus
    ref_allele = reference_allele(vcf_data, chromosome, loci[0])
    genotype_array_ref_numerical = allele_to_ints(genotype_array_ref, ref_allele)

    # compact the numerical genotype data for the first locus
    genotype_array_ref_compacted = list(map(add_if_not_none, genotype_array_ref_numerical))

    rtwo_values = {}

    for locus in loci[1:]:
        # get the phased genotypes in [<char>,<char>] format for the current locus
        genotype_array = fetch_alleles(vcf_data, chromosome, locus)

        # convert the phased genotypes into numerical format for the current locus
        locus_ref = reference_allele(vcf_data, chromosome, locus)
        genotype_array_numerical = allele_to_ints(genotype_array, locus_ref)

        # compact the numerical genotype data for the current locus
        genotype_array_compacted = list(map(add_if_not_none, genotype_array_numerical))

        n_indv_list = []
        for i in zip(genotype_array_ref_compacted, genotype_array_compacted):
            if all([j is not None for j in i]):
                n_indv_list.append(True)
            else:
                n_indv_list.append(False)

        # filter the genotype data for the first locus and the current locus
        genotype_ref_filtered = list(compress(genotype_array_ref_compacted, n_indv_list))
        genotype_filtered = list(compress(genotype_array_compacted, n_indv_list))

        # calculate the r^2 value between the first locus and the current locus
        rtwo_value = genotype_correlation(genotype_ref_filtered, genotype_filtered)

        # store the r^2 value with the locus pair as the key
        rtwo_values[(loci[0], locus)] = rtwo_value

    return rtwo_values

if __name__ == "__main__":
    import sys, argparse

    parser = argparse.ArgumentParser(description = "Calculate the r^2 values between the first locus and multiple other loci in a VCF file.")

    parser.add_argument('-v','--vcf_file', nargs= 1, default=None, help = "The path to the VCF file of interest.")

    parser.add_argument('-l','--loci', nargs= '+', default=None, type=int , help = "The loci of interest. The first locus is used as the reference.")

    parser.add_argument('-c','--chrom',nargs = 1, help = "The chromosome in the file of interest.")

    #get a persistent iterator
    args = parser.parse_args()

    rtwo_values = rtwo_calculator(
        vcf_path = args.vcf_file[0],
        chromosome = args.chrom[0],
        loci = args.loci
    )

    for locus_pair, rtwo in rtwo_values.items():
        print(f"The value for r^2 between loci {locus_pair[0]} and {locus_pair[1]} is: {rtwo}")
