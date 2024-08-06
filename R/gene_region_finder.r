gene_region_finder <- function(chromosome, vcf_file, loci, output_file="/scratch/", distance=500000, window=500000, r2_threshold=0.1) {

  options(scipen=999) # This makes sure that shell commands do not get converted to scientific notation.

  if (!locate_bin_on_shell("tabix")) {
    stop("Please install tabix or make it properly available and try again.")
  }

  if (!locate_bin_on_shell("vcftools")) {
    stop("Please install bcftools or make it properly available and try again.")
  }

	if (!file.exists(vcf_file)){
		stop("The VCF file is not accessible from this location. Did you supply the whole path?")
	}
  # Check if the output file's directory exists, create if not
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # get file base_names for later use
  base_name <- basename(vcf_file)
  base_name <- sub("\\..*$", "", base_name) # this should be the base of the file without the extension

  # check if the .tbi file exists for the current vcf file and give the yes the option to generate it
  proper_tbi(vcf_file)

  # Crunching the numbers for the linkage loci range from input
  snp_start_ld <- loci - distance
  snp_end_ld <- loci + distance

  ## make sure that the start loci is not less than 1 
  if (snp_start_ld < 1) {
    snp_start_ld <- 0
    print("The given loci generated a start point less than 1. Setting it to 0.")
  }

  # generate the base name for the output file using the base name of the input file and ld range
  tabix_file_subset <- paste0(base_name, "_", snp_start_ld, "_", snp_end_ld,".txt")

	tabix_file_subset <- tempfile()

  # asking tabix to subset the file
  system(paste0("tabix -h ", vcf_file, " ", chromosome, ":", snp_start_ld, "-", snp_end_ld, " > ", tabix_file_subset))

  # check to see if the tabix output file has at least 1 line
  # that way we can know that the subsetting did not go awray
  if (length(readLines(tabix_file_subset)) == 0) {
    stop("The program was unable to subset your VCF file for the given input. Please check that the right input file, loci and chromsome name were supplied. Checking for typos might be a good idea, Chr_05 is not the same as Chr05.")
  }

  # Here doc file for tabix snp calculation
  dt <- data.table(
    Chr = c(chromosome),
    Snp = c(loci)
  )

	current_geno_r2_positions_table <- tempfile()

  dt %>%
    fwrite(current_geno_r2_positions_table, sep = "\t", col.names = TRUE)

  # Make a call to the right file

	vcftools_call <- "vcftools"
	
  vcftools_args = c(
    "--vcf", tabix_file_subset,
    "--geno-r2-positions", current_geno_r2_positions_table, # Because the file name is constant we do not need a variable here
    "--ld-window-bp", window,
    "--out", tabix_file_subset
    )

	tryCatch(
		{
			try <- exec_wait(
				vcftools_call,
				args = vcftools_args,
				std_out = TRUE,
				std_err = TRUE
			)
		},
		error = function(e){
            # Custom error message
			print("There was an error with running vcftools on your input.")
            print(paste("Execution attempt produced error:-", e$message))
			stop()
            1 # Return 1 on error
        }
	)
	

  current_ld_file_path = paste0(tabix_file_subset, ".list.geno.ld")

  # check the output of vcftools for length to see that it isn't empty
  if (length(readLines(tabix_file_subset)) == 0) {
    stop("The program was unable to calculate LD for the given loci. Please check why this happened. The program needs LD to generate the final reports.")
  }

  # TODO: delete log file and maybe move subset file
  # read the output of LD calulation into a data table

  current_ld_calculations <- fread(current_ld_file_path, header = TRUE)

  # filter and sort
  sorted_by_ld <- current_ld_calculations %>%
  filter(`R^2` >= r2_threshold, 
         !is.na(`R^2`), 
         `R^2` != "N_INDV") %>%
  arrange(`R^2`)
  
  start <- first(sorted_by_ld$POS2)

  stop <- last(sorted_by_ld$POS2)

  # add current data in a table format

  # assuming new_data is a data.table with the same structure as current_table
  current_output <- list(
    VCF_file = vcf_file,
    Chrom = chromosome,
    locus = loci,
    distance = distance,
    r2_threshold = r2_threshold,
    start = start,
    stop = stop
  )

	return(current_output)

}