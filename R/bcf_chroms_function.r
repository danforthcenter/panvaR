bcf_chroms_func <- function(vcf_file_path){

	# make a temporary file to catch the standard output

	std_output_name <- temp_file(create_file=TRUE) # Defined in the file R/temporary.r
	# The shell call that I want to run is:-
	# `bcftools index -s file.vcf.gz``

	binary_call <- "bcftools"

	binary_args <- c(
		"index","-s",
		vcf_file_path
	)

	tryCatch(
		{
			try <- exec_wait(
				binary_call,
				args = binary_args,
				std_out = std_output_name,
				std_err = FALSE
			)
		},
		error = function(e){
			# Custom error message
			print(paste("Execution attempt produced error:-", e$message))
			1 # Return 1 on error
		}
	)

	if (try == 0){

		final_index_table = fread(std_output_name)

		# get the column names
        bcf_chroms <- final_index_table %>%
            pull(1)

		return(bcf_chroms)
	} else {

		print("There were errors when calculating index for this set of inputs.")
		print("Please read the error message and re-try")
		return(readLines(std_output_name))
	}
    
}