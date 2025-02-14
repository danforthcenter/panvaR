#' return_snplist_for_bp
#'
#' @description 
#' This function will take a Chrom, a BP, a Bef file and return the name of the SNP.
#' This functionality makes it so that we can get the name of the snp accurately and not blindly rely on the BP and Chrom,
#' making our request more reliable.
#' Title
#' A function to return the names of SNPs at specific locations.
#'
#' @param path_to_bed_file Path to the bed file that should have the tag SNP.
#' @param chrom The chromosome of the tag SNP.
#' @param bp The bp of the snp
#' @return The names of the SNP at the bp.
#' 
#' @import tidyverse
#' @import data.table
#' @import parallel
#' @importFrom bigstatsr big_randomSVD
#' @import modelr
#' @export
#'
#' @examples
#' subset_around_tag('<path_to_your_bed_file>')

# A general function to convert bed files to Plink

return_snplist_for_bp <- function(path_to_bed_file, chrom, bp){

	bed_file_input = base_name_func(path_to_bed_file, super_name = TRUE, include_dir = TRUE) # to more accurately get the bed file path
    
	output_file <- temporary_directory()

	snp_list_name <- paste(sample(LETTERS, 5, replace = TRUE), collapse = "")

	output_path <- paste0(output_file,"/",snp_list_name)
	
    binary_call <- "plink2"

	binary_args <- c(
		"--allow-extra-chr",
		"--bfile",bed_file_input,
		"--chr", chrom,
		"--from-bp", bp,
		"--to-bp", bp,
		"--write-snplist",
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

		final_output_list = readLines(paste0(output_path,".snplist"))
		
		return(final_output_list)
    } else{

		print("There were errors when subsetting the bed file around the tag SNP.")
		print("Please read the error message and re-try")
		return(readLines(error_message))
	}
}