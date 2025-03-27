#' apply_maf
#' A function to take a vcf file and filter applies a minor allele frequency filter.
#'
#' @param path_to_vcf_file Path to the vcf file.
#' @param missing_rate (optional) What is the Minor Allele frequency that you want apply?
#' Defaults to 0.05
#' @param new_name (optional) What should be the new basename or prefix for the new file be?
#' Defaults to {base_name}_filtered_{maf_rate}.vcf
#' @return The path to the new vcf file with the missing rate filter applied.
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
#' apply_maf("<path_to_vcf_file>", missing_rate = 0.05)
#' apply_maf <- function(path_to_vcf_file, missing_rate = 0.05, output_name = "new_file")

apply_maf <- function(path_to_vcf_file, maf = 0.05, output_name = NA){

    proper_tbi(path_to_vcf_file) # Make sure you have a tbi file for your VCF file

    # For now I don't think I have a way around the verbose tryCatch wraps for system calls

    # If the output name is not NA then generate a name for the output file

    if(is.na(output_name)){

		base_name = sub("\\.[^.]*$", "", basename(path_to_vcf_file)) # use this code to get the basename of the vcf file without the extension

		dir_name <- dirname(path_to_vcf_file)

		output_name = paste0(dir_name,"/",base_snp, "_window_",window,".vcf",sep = "") # so the new name will be {base_name}_filtered{missing_rate}.vcf

	} else {

		output_name = output_name
	}

	std_output_name = tempfile()

    # calculate the start and the stop value of the window 
    snp_start_ld = base_snp - window
    snp_stop_ld = base_snp + window
    
    # make sure the start is not less than 0
    if (snp_start_ld < 0) {
        snp_start_ld = 0
    }

    # The bash syntax that we want to automate is
    #`tabix -h ${vcf_file} ${chromosome}:${snp_start_ld}-${snp_stop_ld} > ${tabix_output_file}`
    binary_name <- "tabix"

    the_args = c(
       "-h", path_to_vcf_file,
       paste0(chrom,":", snp_start_ld, "-", snp_stop_ld)
    )
    
    tryCatch(
      {
        std_err_file <- tempfile()  # Create a temporary file to store standard error
        try <- exec_wait(
          binary_name,
          args = the_args,
          std_out = std_output_name, # This can be massive and can crash the console 
          std_err = std_err_file # Redirect standard error to a file
        )
        
        # After exec_wait, read the standard error from the file
        std_err_content <- readLines(std_err_file)
        if (length(std_err_content) > 0) {
          print("Standard Error:")
          print(std_err_content)
        }
        
        # Clean up: Remove the temporary file
        unlink(std_err_file)
      },
      error = function(e){
        # Custom error message
        print("There was an error with running vcftools' LD calculation step on your input.")
        print(paste("Execution attempt produced this error:", conditionMessage(e)))
        stop()
        1
      }
    )
	
	if(write_file){
		file.copy(from = std_output_name, to = output_name) # Rijan:- File copy is safer than file move
	} else{
		return(list(std_output_name,output_name))
	}
}