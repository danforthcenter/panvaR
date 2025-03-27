#' subset_aroung_tag
#' A function to take a Bed file and make a window around the tag snp.
#'
#' @param path_to_bed_file Path to the bed file that should have the tag SNP.
#' @param chrom The chromosome of the tag SNP.
#' @param bp The bp of the tag snp
#' @param Window (optional) The radius of the window around the tag SNP
#' Defaults to 500kb
#' @param output_prefix (optional) What should be the new basename or prefix for the new file be?
#' Defaults to {base_name}_windowed
#' @return The path to the new bed file with the missing rate filter applied.
#' 
#' @import tidyverse
#' @import data.table
#' @import parallel
#' @importFrom bigstatsr big_randomSVD
#' @import modelr
#' 
#' @export
#'
#' @examples
#' subset_around_tag('<path_to_your_bed_file>')

# A general function to convert bed files to Plink

subset_around_tag <- function(path_to_bed_file, chrom, bp, window = 500000,output_prefix = NA){

	bed_file_input = base_name_func(path_to_bed_file, super_name = TRUE, include_dir = TRUE)
	
    # make a prefix for output
    if(is.na(output_prefix)){

		# tempoarary directory to store the files

		windowed <- temporary_directory()
		
        base_name_prefix <- 
            base_name_func(path_to_bed_file)

		output_path = paste0(windowed,"/",base_name_prefix,"_windowd_",bp) # Saves the trouble of duplicate names
		
    } else {
        output_path <- output_prefix
    }

	# calculate the start and the stop value of the window 
    snp_start_ld = bp - window
    snp_stop_ld = bp + window
    
    # make sure the start is not less than 0
    if (snp_start_ld < 0) {
        snp_start_ld = 0
    }
	
    binary_call <- "plink2"

	binary_args <- c(
        "--allow-extra-chr",
		"--bfile",bed_file_input,
		"--chr", chrom,
		"--from-bp", snp_start_ld,
		"--to-bp", snp_stop_ld,
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

		print("There were errors when subsetting the bed file around the tag SNP.")
		print("Please read the error message and re-try")
		return(readLines(error_message))
	}
}