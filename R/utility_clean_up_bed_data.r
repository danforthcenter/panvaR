#' bed_file_clean_up
#' A function to take a Bed file and filter for the missing rate and minor allele frequency.
#'
#' @param path_to_bed_file Path to the bed file.
#' @param missing_rate (optional) What is the missing rate that you want to filter for.
#' Defaults to 0.1
#' @param maf (optional) What is the minor allele frequncy to filter for?
#' 
#' @param output_prefix (optional) What should be the new basename or prefix for the new file be?
#' Defaults to {base_name}_filtered{missing_rate}.bed
#' @return The path to the new bed file with the missing rate filter applied.
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
#' bed_file_clean_up('<path_to_your_bed_file>')

# A general function to convert bed files to Plink

bed_file_clean_up <- function(path_to_bed_file, missing_rate = 0.1, maf = 0.05, output_prefix = NA){

	bed_file_input = base_name_func(path_to_bed_file, super_name = TRUE, include_dir = TRUE)
	
    # make a prefix for output
    if(is.na(output_prefix)){

		# tempoarary directory to store the files

		cleanup <- temporary_directory()
		
        base_name_prefix <- 
            base_name_func(path_to_bed_file)

		output_path = paste0(cleanup,"/",base_name_prefix,"_cleaned") # saves the trouble of 
		
    } else {
        output_path <- output_prefix
    }
	# TODO :- If the user supplies a name copy the bed file from scratch to  
	# create a temporary directory 
    binary_call <- "plink2"

	binary_args <- c(
        "--allow-extra-chr",
		"--bfile",bed_file_input,
		"--geno", missing_rate,
		"--maf", maf,
		"--make-bed",
		"--out", output_path
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

		final_output_path = paste0(output_path,".bed")
		
        print(
            paste("The Plink files are available in",final_output_path)
        )
		return(final_output_path)
    } else{

		print("There were errors when applying MAF and missing rate filter to your BED file.")
		print("Please read the error message and re-try")
		return(readLines(error_message))
	}
}