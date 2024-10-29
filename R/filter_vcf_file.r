#' filter_vcf_file
#'
#' @description 
#' This function takes a (tab-delimited) table of formatted with 
#' the chromosome in the first tables and the bp in the second.
#' It takes that table and filters a vcf file for those bps.
#' 
#' Title
#' A function to return a list of BPs with LD beyond assigned r2 value
#'
#' @param vcf_file_path Path to the vcf file that should have the tag SNP.
#' @param keep_table The table of BPs to keep, either a table or a tabular object
#'
#' @return Path to a vcf file that with only the filtered set of BPs
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
#' ld_filtered_snp_list("<vcf_file_path>", chrom = "<chorm>", bp = <bp_value>, r2_threshold = 0.6)

# A general function to convert vcf files to Plink

filter_vcf_file <- function(vcf_file_path, keep_table, output_prefix = NA){

	proper_tbi(vcf_file_path)

	 # make a prefix for output
    if(is.na(output_prefix)){

		# tempoarary directory to store the files

		vcf_file_dir <- temporary_directory()
		
        base_name_prefix <- 
            base_name_func(vcf_file_path)

		subset_table_name <- paste(sample(LETTERS, 5, replace = TRUE), collapse = "") # A random letter generator	

		output_path = paste0(vcf_file_dir,"/",base_name_prefix,subset_table_name,".vcf")
		
    } else {
        output_path <- output_prefix
    }

	# because VCF files can be either gz or unzipped
	if(endsWith(vcf_file_path,".gz")){
		file_type = "--gzvcf"
	} else if(endsWith(vcf_file_path,".vcf")){
		file_type = "--vcf"
	}
    

	keep_table_path <- keep_table_sanitizer(keep_table)

	# get the name of the snp using `return_snplist_for_bp`
	# this is necessary because plink2 does not have the option to supply 
	# just bp and chromosomes
	
    binary_call <- "vcftools"
	
	# vcftools 

	binary_args <- c(
		file_type, vcf_file_path,
		"--positions", keep_table_path,
		"--recode", "--recode-INFO-all",
		"--out", output_path
	)
	
    tryCatch(
        {
			error_message <- tempfile()
            try <- exec_wait(
                binary_call,
                args = binary_args,
                std_out = TRUE,
                std_err = error_message
            )
        },
        error = function(e){
            # Custom error message
            print(paste("Execution attempt produced error:-", e$message))
            1 # Return 1 on error
        }
    )

    if(try == 0){ 
		final_path = paste0(output_path,".recode.vcf") # this is something VCF tools adds unfortunately
		return(final_path)
    } else{

		print("There were errors when calculating ld for this set of inputs.")
		print("Please read the error message and re-try")
		return(readLines(error_message))
	}
}