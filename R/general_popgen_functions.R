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
      
      binary_name <- "tabix"
      the_args <- c("-p", "vcf", vcf_file)
      
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

#######
# A function to generate a table with LD against a tag SNP
generate_ld_table <- function(path_to_vcf_file, chrom, snp_bp, window = 500000) {

	# The goal of this function is to take a vcf file
	# a chrom and a bp position
	# and return a table of LD values

	proper_tbi(path_to_vcf_file) # Make sure that the vcf file has a tbi file

	# Here doc file for tabix snp calculation
  	dt <- data.table(
  	  Chr = c(chrom),
  	  Snp = c(snp_bp)
  	)

	current_geno_r2_positions_table <- tempfile()

  	dt %>%
  	  fwrite(current_geno_r2_positions_table, sep = "\t", col.names = TRUE) # We will use this as the "Here doc" that we pass to tabix for the snp file

	# set the system call to vcftools and all the variables needed for it to work.
	tabix_file_subset <- tempfile() # The file that will hold the output from vcftools

	print(tabix_file_subset)

	binary_name <- "vcftools"

    the_args = c(
       "--vcf", path_to_vcf_file,
       "--geno-r2-positions", current_geno_r2_positions_table, 
       "--ld-window-bp", window,
       "--out", tabix_file_subset
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
        # Clean up: Remove the temporary file
        unlink(std_err_file)
      },
      error = function(e){
        # Custom error message
        print("There was an error with running vcftools' LD calculation step on your input.")
        print(paste("Execution attempt produced this error:", conditionMessage(e)))
      }
    )

	# read the table in tabix_file_subset but we need to add the right extenion to it
	# the automated extension addition is a quirk of vcftools
	current_ld_file_path = paste0(tabix_file_subset, ".list.geno.ld")

	# read the current_ld_file_path into a data.table
	current_ld_table <- fread(current_ld_file_path)

	return(current_ld_table)
}

#### 
# A function to break a VCF file into a window supplied by the user

# The goal of this function is to take a vcf_file , a chrom, a base_snp and subset it 
# This makes PLINK operations easier and faster

vcf_window_subset <- function(path_to_vcf_file, chrom, base_snp, window = 500000, output_name = NA){

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
          std_out = output_name, # This can be massive and can crash the console 
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


    return(output_name)
}