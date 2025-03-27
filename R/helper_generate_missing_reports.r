options(scipen=999) # This makes sure that scientific notation does not supplant natural numbers
# This R script is meant to set the default values for the supplied list of inputs

# function will generate and return two data.tables that hold the missing information for the vcf files that you supply

#' generate_missing_reports
#' @description
#' This function takes a vcf file and generate the missing reports for both lines and SNPs
#'
#' @param path_to_vcf_file Path to the vcf file.
#' @return A list containing the following elements:
#'   \item{table_for_snps}{A table that shows the missing values for the SNPS in the given VCF file.}
#'   \item{table_for_lines}{A table that shows the missing values for the lines in the given VCF file.}
#' 
#' @import tidyverse
#' @import data.table
#' @import sys
#' @import parallel
#' @import bigsnpr
#' @import modelr
#'
#' @export
#'
#' @examples
#' generate_missing_report("/path/to/your/vcf/file.vcf")

generate_missing_report <- function(path_to_vcf_file){

	# Check if the vcf file supplied has a tbi file
	proper_tbi(path_to_vcf_file) # This function was defined in `general_functions`

	# The name of the binary that will be passed to the sys package, relevant down the file
	plink2_call <- "plink2"

	# varaint file directory
	missing_reports_directory <- temporary_directory()

	# general output basename
	# This is where we will receive the report files
	# you can pass just the basename becasue plink2 will automatically append the relevant extension
	# .vmiss for variant
	# .smiss for sample

	plink2_basename <- paste0(missing_reports_directory, "/", "plink2")

	# Generate the information for the variants at hand

	plink_args_varaints <- c(
    	"--vcf", path_to_vcf_file,
    	"--missing", "variant-only",
    	"--out", plink2_basename
	)

	tryCatch(
    	{
    	    try <- exec_wait(
    	        plink2_call,
    	        args = plink_args_varaints,
    	        std_out = FALSE,
    	        std_err = TRUE
    	    )
    	},
    	error = function(e){
    	    # Custom error message
    	    print(paste("The attempt to generate the missing report for variant failed, please check this error message:-", e$message))
    	    1 # Return 1 on error
    	}
    )

	# args for samples
	plink_args_samples <- c(
    	"--vcf", path_to_vcf_file,
    	"--missing", "sample-only",
    	"--out", plink2_basename
	)

	tryCatch(
		{
		    try <- exec_wait(
		        plink2_call,
		        args = plink_args_samples,
		        std_out = FALSE,
		        std_err = TRUE
		    )
		},
		error = function(e){
		    # Custom error message
		    print(paste("The attempt to generate the missing report for sample failed, please check this error message:-", e$message))
		    1 # Return 1 on error
		}
	)

	# generate paths for the vmiss and smiss files that now sit in the temporary directory
	vmiss_path <- paste0(plink2_basename, ".vmiss")
	smiss_path <- paste0(plink2_basename, ".smiss")

	variant_table <- fread(vmiss_path)
	sample_table <- fread(smiss_path)

	return(list(table_for_snps = variant_table, table_for_lines = sample_table))
}