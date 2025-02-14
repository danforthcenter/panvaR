#' vcf_to_plink
#'
#' @description 
#' This is a function that takes a vcf file and converts it to Plink format.
#'
#' @param vcf_file Path to the input VCF file.
#' @param output The path where the Plink file should be written.
#' Defaults to writing the output file in a `tempdir`. If you want the file to be kept around longer supply an accessible path.
#' 
#' @return The path to the new Bed file.
#'
#' @examples
#' vcf_window_subset("path/to/your_file.vcf",chrom = "Chr_001", base_snp = 6857045, output_name = "windowed_vcf_file")
#' 
#' @import tidyverse
#' @import data.table
#' @import parallel
#' @importFrom bigstatsr big_randomSVD
#' @import modelr
#' 
#' @export

# A general function to convert vcf files to Plink

vcf_to_plink2 <- function(vcf_file_path, output_prefix = NA){

    # make a prefix for output
    if(is.na(output_prefix)){

		# tempoarary directory to store the files

		plink2_tempdir <- temporary_directory()
		
        base_name_prefix <- 
            base_name_func(vcf_file_path)

		output_path = paste0(plink2_tempdir,"/",base_name_prefix)
		
    } else {
        output_path <- output_prefix
    }
	# TODO :- If the user supplies a name copy the bed file from scratch to  
	# create a temporary directory 


    # check if bed and bim files already exist
    bim_file_path = paste0(output_path,".bim")
    bed_file_path = paste0(output_path,".bed")

    if(file.exists(bim_file_path) && file.exists(bed_file_path)){
        # return the path to the bed file
        print("It looks like your bed and bim files already exist. This was done with heuristics and regex.")

        print("The heuristics could be wrong, if that's the case please delete them and try again.")

        return(list(bed = bed_file_path,bim = bim_file_path))
    }
    binary_call <- "plink2"

	binary_args <- c(
        "--allow-extra-chr",
		"--vcf",
		vcf_file_path, "--make-bed",
		"--set-all-var-ids", "Chr_@_BP_#",
		"--out", output_path
	)

	# Rijan: Given how legible Plink2's error messages are
	# Rijan: I think a simple STDOUT catch is good enough

	# Rijan: By the looks of it Plink2 can handle .gz and .vcf paths just fine
	
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

		bed_file_path = paste0(output_path,".bed")

        bim_file_path = paste0(output_path,".bim")
		
        print(
            paste("The Plink files are available in",bed_file_path)
        )
		return(list(bed = bed_file_path,bim = bim_file_path))
    } else{

		print("Your VCF file could not be converted to Plink2's native format.")
		print("Please read the error message and re-try")
		return(readLines(error_message))
	}
}
