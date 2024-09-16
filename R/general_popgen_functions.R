# Make a tbi file for a given VCF file if it does not already exist
proper_tbi <- function(vcf_file) {

  # get the extension of the supplied vcf file
  vcf_file_extenion <- extention_func(path_to_file = vcf_file)

  if (vcf_file_extenion != "gz") {
    stop(
      "The file you have supplied is not compressed using BGZIP and BCFtools cannot work with it. \n 
      Please compress your file using BGZIP and try again.
      "
    )
  }

  if (!file.exists(paste0(vcf_file, ".tbi"))) {
    # ask the user if they want to generate the .tbi file
    message("The .tbi file for the given vcf file does not exist - a .tbi file is needed for panvaR to function. Would you like to generate it? (y/n) ")
    answer <- tolower(readline())
    if (answer == "y") {
      # generate the .tbi file
      message("Asking BCFtools to generate the .tbi file...")
      
      binary_name <- "bcftools"
      the_args <- c("index", "-t", vcf_file)
      
      tryCatch(
        {
          std_err_file <- tempfile()  # Create a temporary file to store standard error
          try <- exec_wait(
            binary_name,
            args = the_args,
            std_out = NULL,  # We don't need to capture standard output
            std_err = std_err_file  # Redirect standard error to a file
          )
          
          # After exec_wait, read the standard error from the file
          std_err_content <- readLines(std_err_file)
          if (length(std_err_content) > 0) {
            message("Please read this error produced by tabix:")
            message(paste(std_err_content, collapse = "\n"))
          }
          
          # Clean up: Remove the temporary file
          unlink(std_err_file)
        },
        error = function(e) {
          # Custom error message
          message("There was an error while generating the .tbi file.")
          message(paste("Execution attempt produced this error:", conditionMessage(e)))
          stop("Failed to generate .tbi file. Please check the error message and try again.")
        }
      )
      return(TRUE)
    } else {
      stop("The .tbi file for the given vcf file does not have a tabix index. Please generate it and try again.")
    }
  } else {
    return(TRUE)
  }
}


# A general function to convert vcf files to Plink
vcf_to_bed <- function(vcf_file_path, output_prefix = NA){

    # make a prefix for output
    if(is.na(output_prefix)){
        base_name_prefix <- 
            strsplit(vcf_file_path,".", fixed = TRUE)[[1]][1]
    } else {
        base_name_prefix <- output_prefix
    }
    
    plink_call <- "plink2"

    plink_args <- c("--vcf",vcf_file_path, "--make-bed","--out", base_name_prefix)

    tryCatch(
        {
            try <- exec_wait(
                plink_call,
                args = plink_args,
                std_out = TRUE,
                std_err = FALSE
            )
        },
        error = function(e){
            # Custom error message
            print(paste("Execution attempt produced error:-", e$message))
            1 # Return 1 on error
        }
    )

    if(try == 0){
        print(
            paste("The Plink files are available in",base_name_prefix)
        )
    }
}

# to call grep
tabix_grep_call <- function(subset_file, gene_name){

	grep_call <- "grep"

    grep_output <- tempfile()

	the_args <- c(
	  "-E",  # Extended regex patterns
	  paste0("#|", gene_name),  # Pattern to match
        subset_file
	)

	tryCatch(
	  {
	    try <- exec_wait(
	      grep_call,
	      args = the_args,
	      std_out = grep_output,
	      std_err = TRUE
	    )
	  },
	  error = function(e) {
	    print(paste("There was an error with running grep on your input."))
	    print(paste("Execution attempt produced error:-", e$message))
	    stop()
	    1  # Return 1 on error
	  }
	)
    
    return(grep_output)
}


# To call tabix
tabix_subset_polymorphisms <- function(start, stop, chrom, gene_name, vcf_file){

    proper_tbi(vcf_file) # It is possible for the vcf file to not have a tbi file at this point. 

    first_tabix_output <- tempfile()
    
    tabix_call <- "tabix"
	
    the_args = c(
    "-h", vcf_file,
    paste0(chrom,":",start, "-", stop)
    )

	tryCatch(
		{
			try <- exec_wait(
				tabix_call,
				args = the_args,
				std_out = first_tabix_output,
				std_err = TRUE
			)
		},
		error = function(e){
            # Custom error message
			print("There was an error with running grep on your input.")
            print(paste("Execution attempt produced error:-", e$message))
			stop()
            1 # Return 1 on error
        }
	)

    grep_filtered <- tabix_grep_call(first_tabix_output, gene_name)

    return(grep_filtered)
    
}