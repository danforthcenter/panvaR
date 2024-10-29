#' pca_martix_generator
#'
#' @description 
#' This function takes either a vcf file or a bed file and generates a PCA matrix.
#'
#' @param genotype_data Path to the input genotype data.
#' Either a vcf file (vcf/gz) or bed file
#' @param pc_count (optional) How many principle components do you want?
#' Defaults to 5
#' @param maf (optional) The minor allele frequency that will be applied to the genotype data.
#' Defaults to 0.05
#' @param missing_rate (optional) The missing_rate filter that will be applied to the genotype data.
#' Defaults to 0.1
#' @return A matrix with your principle components.
#'
#' @examples
#' pca_matrix_generator("path/to/genotype_data",pc_count = 5, maf = 0.05, missing_rate = 0.1)
#' 
#' @import tidyverse
#' @import data.table
#' @import sys
#' @import parallel
#' @import bigsnpr
#' @import modelr
#' 
#' @export

# A function to make PCAs for your genotype data.

pca_matrix_generator <- function(genotype_data,pc_count = 5, maf = 0.05, missing_rate = 0.1) {

	# Get the core count from the ergonomics set of code
	core_count = good_core_count()

	# get the extention of the genotype_data file
	genotype_data_format <- extention_func(genotype_data)

	# if the format is gz or vcf send to plink2
	# else send the genotype data to be filtered for maf and missing rate filter
	if(genotype_data_format == "gz"){
		genotype_data_right_format_base <- vcf_to_plink2(genotype_data) # This returns a list of bed and bim file
		genotype_data_right_format <- genotype_data_right_format_base$bed
		genotype_data_right_format <- genotype_data_right_format_base$bed

		# Now clean up the data
		cleaned_up_bed_file <- bed_file_clean_up(genotype_data_right_format,maf = maf, missing_rate = missing_rate)
	} else if (genotype_data_format == "vcf"){
		genotype_data_right_format_base <- vcf_to_plink2(genotype_data) # This returns a list of bed and bim file
		genotype_data_right_format <- genotype_data_right_format_base$bed

		# Now clean up the data
		cleaned_up_bed_file <- bed_file_clean_up(genotype_data_right_format,maf = maf, missing_rate = missing_rate)
	} else {

		# if the data is already in plink2 format
		# Clean up for the missing rate and the maf
		cleaned_up_bed_file <- bed_file_clean_up(genotype_data,maf = maf, missing_rate = missing_rate)
	}

	# Read the bam file into using `snp_readBed`
	bed_file_readout <- snp_readBed(
	    cleaned_up_bed_file
	)

	rds_file_base <- base_name_func(cleaned_up_bed_file,include_dir = TRUE)

	# the path to the rds file made by `rds_file_base`
	rds_file_path <- paste0(c(rds_file_base,".rds"),collapse = "")

	# Read the genotype using `snp_attach`
	genotype_rds_data <- snp_attach(rds_file_path)

	# The genotypes
	the_genotypes <- genotype_rds_data$genotypes

	# The vector of chromosomes
	the_chromosomes <- genotype_rds_data$map$chromosome

	# The vector of input basepairs
	the_bp <- genotype_rds_data$map$physical.pos

	# How take the list of existing chromosomes 
	# Make a table that mapes the existing chromsomes to a unqiue numerical ID
	# This saves us the trouble chars causing issues
	chromosomes_as_ints <- the_chromosomes %>%
	  enframe(name = NULL, value = "chromosome") %>%
	  mutate(chr_int = dense_rank(chromosome)) %>%
	  pull(chr_int)

	# Make the exclusion index  
	# Rijan: Reading material clumping and pruning here https://www.biostars.org/p/343818/

	exclusion_index <- snp_indLRLDR(
	    infos.chr = chromosomes_as_ints, 
	    infos.pos = the_bp
	)

	# Rijan: Reading material clumping and pruning here https://www.biostars.org/p/343818/
	inclusion_index <- snp_clumping(
	    the_genotypes, 
	    infos.chr = chromosomes_as_ints, 
	    exclude = exclusion_index,
	    ncores = core_count
	)

	# Apply a PCA using an algorithm optimized for large file backed matrices
	big_random_pca <- big_randomSVD(
	    the_genotypes, 
	    fun.scaling = snp_scaleBinom(),
	    ind.col = inclusion_index,
	    ncores = core_count
	)

	the_PCs <- predict(big_random_pca)

	the_covariates <- the_PCs[,1:pc_count]

	return(the_covariates)
}