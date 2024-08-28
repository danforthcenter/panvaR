# A function to apply default missing values rates to vcf files using plink2

#' Title
#' A function to take a vcf file and filter it for missing values.
#'
#' @param path_to_vcf_file Path to the vcf file.
#' @param missing_rate (optional) What is the missing rate that you want to filter for.
#' Defaults to 0.1
#' @param new_name (optional) What should be the new basename or prefix for the new file be?
#' Defaults to {base_name}_filtered{missing_rate}.vcf
#' @return The path to the new vcf file with the missing rate filter applied.
#' @export
#'
#' @examples
#' healthy_maf_values("your_vcf_file_here.vcf", missing_rate = 0.10, new_name = "tidy_organized_name.vcf")

healthy_maf_values <- function(path_to_vcf_file, missing_rate = 0.1, new_name = NA){

	# Check if the vcf file supplied has a tbi file
	proper_tbi(path_to_vcf_file) # This function was defined in `general_popgen_functions`

	# The code to remove missing values from vcf files is
	# if new_name is NA then make a new name for the output that will be produced
	if(is.na(new_name)){

		base_name = base_name_func(path_to_vcf_file) # use this code to get the basename of the vcf file without the extension

		new_name = paste0(base_name, "_filtered_for_",missing_rate,sep = "") # so the new name will be {base_name}_filtered{missing_rate}.vcf
	} else {

		new_name = new_name
	}
	print(base_name)
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
    	        std_out = TRUE,
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