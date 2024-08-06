# Is the binary you want to use available in the shell?
locate_bin_on_shell <- function(command){

    the_command <- paste0(command)

    result <- tryCatch(
        {
            test <- sys::exec_wait(command,
            args = "-v",
            std_out = TRUE,
            std_err = FALSE
            )
            0 # Return 0 on success
        },
        error = function(e) {
            # Custom error message
            print(paste0(command, " is either not installed or not available in this shell."))
            print(paste("Execution attempt produced error:-", e$message))
            1 # Return 1 on error
        }
    )

    if (result == 0){
        return(TRUE)
    }
	else{
		return(FALSE)
	}
}


# Make a tbi file for a given VCF file if it does not already exist
proper_tbi <- function(vcf_file) {
  if (!file.exists(paste0(vcf_file, ".tbi"))) {
    # ask the user if they want to generate the .tbi file
    message("The .tbi file for the given vcf file does not exist. Would you like to generate it? (y/n) ")
    answer <- tolower(readline())
    if (answer == "y") {
      # generate the .tbi file
      message("Asking tabix to generate the .tbi file...")
      system(paste0("tabix -p vcf ", vcf_file))
      return(TRUE)
    } else {
      stop("The .tbi file for the given vcf file does not have a tabix index. Please generate it and try again.")
    }
  } else {
    return(TRUE)
  }
}


# Count and provide the number of cells that should be used
good_core_count <- function(){
	
	core_count = nb_cores() - 1 
	if(core_count < 1){
	    core_count = 1
	}

	return(core_count)
}


# A general function to convert vcf files to Plink
vcf_to_bed <- function(vcf_file_path, output_prefix = "auto"){

    # make a prefix for output
    if(output_prefix == "auto"){
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