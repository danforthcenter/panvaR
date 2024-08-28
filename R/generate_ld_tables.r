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

	# This will hold the Here Doc
	current_geno_r2_positions_table <- tempfile()
	
	# get the base of the input file
	base_prefix <- base_name_func(path_to_vcf_file, super_name = TRUE)

	# create a directory to hold and isolate current data
	output_dir <- tempdir()

	# what will the name of the ouptut file be?
	output_path_prefix <- paste0(output_dir,"/",base_prefix, sep = "")

  	dt %>%
  	  fwrite(current_geno_r2_positions_table, sep = "\t", col.names = TRUE) # We will use this as the "Here doc" that we pass to tabix for the snp file

	binary_name <- "vcftools"

    the_args = c(
       "--vcf", path_to_vcf_file,
       "--geno-r2-positions", current_geno_r2_positions_table, 
       "--ld-window-bp", window,
       "--out", output_path_prefix
    )
    
    tryCatch(
    {
      std_err_file <- tempfile()  # Create a temporary file to store standard error
      exit_status <- exec_wait(
        binary_name,
        args = the_args,
        std_out = FALSE, 
        std_err = std_err_file # Redirect standard error to a file
      )
      
      # Check exit status
      if (exit_status != 0) {
        # Read stderr content
        stderr_content <- readLines(std_err_file)
        stop(sprintf("VCFtools exited with non-zero status: %d\nStandard Error Output:\n%s", 
                     exit_status, paste(stderr_content, collapse = "\n")))
      }
      
      # Verify output file exists
      if (!file.exists(current_ld_file_path)) {
        # Read stderr content
        stderr_content <- readLines(std_err_file)
        stop(sprintf("Expected output file was not created: %s\nStandard Error Output:\n%s", 
                     current_ld_file_path, paste(stderr_content, collapse = "\n")))
      }
      
      # Clean up: Remove the temporary file
      unlink(std_err_file)
    },
    error = function(e) {
      cat("There was an error with running vcftools' LD calculation step on your input.\n")
      cat("Error message:", conditionMessage(e), "\n")
      
      # If std_err_file exists and hasn't been deleted, print its contents
      if (file.exists(std_err_file)) {
        stderr_content <- readLines(std_err_file)
        cat("Standard Error Output:\n")
        cat(paste(stderr_content, collapse = "\n"), "\n")
        unlink(std_err_file)  # Clean up the file after printing
      }
    }
  )

	# read the table in tabix_file_subset but we need to add the right extenion to it
	# the automated extension addition is a quirk of vcftools
	current_ld_file_path = paste0(output_path_prefix, ".list.geno.ld")

	# read the current_ld_file_path into a data.table
	current_ld_table <- fread(current_ld_file_path)

	return(current_ld_table)
}