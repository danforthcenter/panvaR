# This R script is meant to set the default values for the supplied list of inputs

# function will generate and return two data.tables that hold the missing information for the vcf files that you supply

generate_missing_report <- function(path_to_vcf_file){

	# Check if the vcf file supplied has a tbi file
	proper_tbi(path_to_vcf_file) # This function was defined in `general_functions`

	# The name of the binary that will be passed to the sys package, relevant down the file
	plink2_call <- "plink2"

	# varaint file directory
	missing_reports_directory <- tempdir()

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

	return(list(v_table = variant_table, s_table = sample_table))
}

# A function to apply default missing values rates to vcf files using plink2

healthy_maf_values <- function(path_to_vcf_file, missing_rate = 0.05, new_name = NA){

	# Check if the vcf file supplied has a tbi file
	proper_tbi(path_to_vcf_file) # This function was defined in `general_popgen_functions`

	# The code to remove missing values from vcf files is
	# if new_name is NA then make a new name for the output that will be produced
	if(new_name == NA){

		base_name = sub("\\.[^.]*$", "", basename(current_file)) # use this code to get the basename of the vcf file without the extension

		new_name = paste0(base_name, "_filtered",missing_rate,".vcf",sep = "") # so the new name will be {base_name}_filtered{missing_rate}.vcf
	} else {

		new_name = new_name
	}
	
	# The basg code to remove missing values from vcf files is
	# `plink2 --vcf <input_file.vcf> --geno 0.1 --make-pgen --out <filtered_output_name>`
	# we need to translate this such that the sys package can use it

	# The name of the binary that will be passed to the sys package, relevant down the file
	plink2_call <- "plink2"

	# Generate the information for the variants at hand

	plink_maf_args <- c(
    	"--vcf", path_to_vcf_file,
    	"--geno", missing_rate,
    	"--make-pgen", "--out",
		new_name
	)

	tryCatch(
    	{
    	    try <- exec_wait(
    	        plink2_call,
    	        args = plink_maf_args,
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

	# Return the new name of the file
	print(paste0(
		"The new name of the file is: ",
		new_name
	))

	# Return the new name
	return(new_name)

}

vcf_window_subset <- function(path_to_vcf_file, chrom, base_snp, window = 500000, output_name = NA){

    proper_tbi(path_to_vcf_file) # Make sure you have a tbi file for your VCF file

    # For now I don't think I have a way around the verbose tryCatch wraps for system calls

    # If the output name is not NA then generate a name for the output file

    if(is.na(output_name)){

		base_name = sub("\\.[^.]*$", "", basename(path_to_vcf_file)) # use this code to get the basename of the vcf file without the extension

		dir_name <- dirname(path_to_vcf_file)

		output_name = paste0(dir_name,"/",base_snp, "_window_",window,sep = "") # so the new name will be {base_name}_filtered{missing_rate}.vcf

	} else {

		output_name = output_name
	}


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
       "--h", path_to_vcf_file,
       paste0(chrom,":", snp_start_ld, "-", snp_stop_ld), 
       ">", output_name
    )
    
    tryCatch(
      {
        std_err_file <- tempfile()  # Create a temporary file to store standard error
        try <- exec_wait(
          binary_name,
          args = the_args,
          std_out = FALSE, 
          std_err = std_err_file # Redirect standard error to a file
        )
        
        # After exec_wait, read the standard error from the file
        std_err_content <- readLines(std_err_file)
        print("Execution completed.")
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


    return(output_name)
}