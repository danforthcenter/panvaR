#' ld_filtered_snp_list
#'
#' @description 
#' This function will take a Chrom, a BP, a Bed file, a r2 threshold and return
#' for that tag bp/tag SNP, a list of BPs that have are in LD beyond the r2 threshold.
#' Title
#' A function to return a list of BPs with LD beyond assigned r2 value
#'
#' @param path_to_bed_file Path to the bed file that should have the tag SNP.
#' @param chrom The chromosome of the tag SNP.
#' @param bp The bp of the snp
#' @param (optional) The r2 threshold that you want to work with.
#' Defaults to 
#' @return A table with the list of BPs that cross the assigned threshold.
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
#' ld_filtered_snp_list("<bed_file_path>", chrom = "<chorm>", bp = <bp_value>, r2_threshold = 0.6)

# A general function to convert bed files to Plink

ld_filtered_snp_list <- function(path_to_bed_file, chrom, bp, r2_threshold = 0.5){

	bed_file_input = base_name_func(path_to_bed_file, super_name = TRUE, include_dir = TRUE) # to more accurately get the bed file path
    
	output_file <- temporary_directory()

	r2_table_name <- paste(sample(LETTERS, 5, replace = TRUE), collapse = "") # A random letter generator

	output_path <- paste0(output_file,"/",r2_table_name)

	# get the name of the snp using `return_snplist_for_bp`
	# this is necessary because plink2 does not have the option to supply 
	# just bp and chromosomes
	
	snp_name <- return_snplist_for_bp(path_to_bed_file=path_to_bed_file, chrom = chrom, bp = bp)

	print(paste0("snp_name is",snp_name))
	print(paste0("snp length is", length(snp_name)))
	
    binary_call <- "plink2"
	
	binary_args <- c(
		"--allow-extra-chr",
		"--bfile",bed_file_input,
		"--ld-snp", snp_name,
		"--r2-phased", "--ld-window-kb",
		9999999, "--ld-window-r2", 
		r2_threshold, "--out",
		output_path
	)

	# Rijan: Given how legible Plink2's error messages are
	# Rijan: I think a simple STDOUT catch is good enough

	# Rijan: By the looks of it Plink2 can handle .gz and .bed paths just fine
	
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

		final_ld_table = fread(paste0(output_path,".vcor"))

		colnames(final_ld_table) = c(
			"Tag_snp_chrom",
			"Tag_snp_bp",
			"Tag_snp_id",
			"Subject_snp_chrom",
			"Subject_snp_bp",
			"Subset_snp_id",
			"Phased_r2"
		)

		# Just to safe-guard against being bombarded with all the snps
		confirmed_final_ld_table <- final_ld_table %>% 
			filter(Tag_snp_bp == bp)
		
		return(confirmed_final_ld_table)
    } else{

		print("There were errors when calculating ld for this set of inputs.")
		print("Please read the error message and re-try")
		return(readLines(error_message))
	}
}