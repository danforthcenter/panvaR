# # A function to return just the basename of a file
# base_name_func <- function(path_to_file,super_name = FALSE, include_dir = FALSE){
#   
#   # use this code to get the basename of the file without the extension
#   name_string_items <- strsplit(basename(path_to_file),"\\.")[[1]]
#   
#   #this gets everything other than the last item
#   # this is also a good fail-safe against supplying a file that has no extension
#   name_items <- name_string_items[-length(name_string_items)]
#   
#   # do you need the directory name of the file
#   dir_name <- dirname(path_to_file)
#   
#   if (super_name){
#     the_base = paste(name_items, collapse = ".")
#   } else {
#     the_base = name_items[1]
#   }
#   
#   if (include_dir){
#     output_name = paste0(dir_name,"/",the_base,sep = "")
#   } else {
#     output_name = the_base
#   }
#   
#   return(output_name)
# }
# 
# # This function get the extention of the file supplied
# extention_func <- function(path_to_file,super_name = FALSE, include_dir = FALSE){
#   
#   # use this code to get the basename of the file without the extension
#   name_string_items <- strsplit(basename(path_to_file),"\\.")[[1]]
#   
#   # Just get the last item in the list
#   extension <- name_string_items[length(name_string_items)]
#   
#   return(extension)
# }
# 
# 
# # Is the binary you want to use available in the shell?
# locate_bin_on_shell <- function(command){
#   
#   the_command <- paste0(command)
#   
#   result <- tryCatch(
#     {
#       test <- sys::exec_wait(command,
#                              args = "-v",
#                              std_out = TRUE,
#                              std_err = FALSE
#       )
#       0 # Return 0 on success
#     },
#     error = function(e) {
#       # Custom error message
#       print(paste0(command, " is either not installed or not available in this shell."))
#       print(paste("Execution attempt produced error:-", e$message))
#       1 # Return 1 on error
#     }
#   )
#   
#   if (result == 0){
#     return(TRUE)
#   }
#   else{
#     return(FALSE)
#   }
# }
# 
# # Count and provide the number of cells that should be used
# good_core_count <- function(){
#   
#   core_count = detectCores() - 1 
#   if(core_count < 1){
#     core_count = 1
#   }
#   
#   return(core_count)
# }
# 
# # A function to get all the BPs to keep
# snps_to_keep <- function(ld_table){
#   
#   # This function only takes the very striclty defined
#   # and shaped table as made by the function `ld_filtered_snp_list`
#   # in the file `ld_filtered_snp_list.R`
#   
#   # This function is needed because the way
#   # Plink2 calculates LD, it leaves out the tag snp
#   # and keeps everything else
#   # that is all that we are accomodating for.
#   
#   # Get other snps
#   other_snps <- ld_table %>%
#     select(.data$Subject_snp_chrom, .data$Subject_snp_bp)  %>%
#     rename(CHROM = .data$Subject_snp_chrom, BP = .data$Subject_snp_bp)
#   
#   # Get tag snps
#   tag_snp <- ld_table %>%
#     select(.data$Tag_snp_chrom, .data$Tag_snp_bp) %>%
#     distinct() %>%
#     rename(CHROM = .data$Tag_snp_chrom, BP = .data$Tag_snp_bp)
#   
#   keep_list <- rbind(other_snps,tag_snp)
#   
#   return(keep_list)
# }
# 
# # Sanitize the tables to keep
# keep_table_sanitizer <- function(table) {
#   
#   keep_table_path <- temp_file(prefix = "keep_table")
#   
#   if(is.character(table)){
#     
#     keep_table <- fread(table)
#     
#     # TODO:- How to make sure the Chrs -
#     # Come before the loci in the the -
#     # keep tables?
#     
#     keep_table %>%
#       fwrite(keep_table_path,sep = "\t",col.names = FALSE)
#   } else{
#     keep_table <- as.data.table(table,col.names=TRUE)
#     
#     keep_table <- keep_table %>%
#       select(.data$CHROM, .data$BP)
#     
#     keep_table %>%
#       fwrite(keep_table_path,sep = "\t",col.names = FALSE)
#   }
#   
#   return(keep_table_path)
#   
# }
# 
# # A function to get the tag snp out of the gwas results table
# tag_snp_func <- function(gwas_results){
#   
#   # This function assumes that the supplied tables has the following fields
#   # 1. 'the_bp`
#   # 2. 'the_chromosomes'
#   # 3. 'the_pvalues' -> -log10() of the pvalues 
#   
#   # Just double checking to make sure that the table is arranged properly
#   
#   current_table <- gwas_results %>% 
#     arrange(desc(as.numeric(.data$Pvalues)))
#   
#   tag_snp_row <- current_table %>%
#     slice(1)
#   
#   tag_snp <- list(
#     tag_snp_bp = as.integer(tag_snp_row$BP), 
#     tag_snp_chromosome = tag_snp_row$CHROM
#     
#   )
#   
#   return(tag_snp)
# }
# 
# # The user might supply the tag snp formatted as - 
# # "<chr>:<bp>". This function provides the convienience -
# # to make the split easily.
# 
# tag_snp_splitter <- function(tag_snp){
#   
#   tag_snp_items <- strsplit(tag_snp,":")[[1]]
#   
#   tag_snp_chrom <- tag_snp_items[1]
#   
#   tag_snp_bp <- as.double(tag_snp_items[2])
#   
#   return_list <- list(
#     chrom = tag_snp_chrom,
#     bp = tag_snp_bp
#   )
#   
#   return(return_list)
# }
# 
# # Make the table of LD that retains the tag SNP
# ld_table_maker <- function(ld_table){
#   
#   ld_table_subject <- ld_table %>% 
#     select(
#       .data$Subject_snp_chrom,
#       .data$Subject_snp_bp,
#       .data$PHASED_R2
#     ) %>%
#     rename(CHROM = .data$Subject_snp_chrom, BP = .data$Subject_snp_bp, LD = .data$PHASED_R2)
#   
#   ld_table_tag <- ld_table %>%
#     select(.data$Tag_snp_chrom, .data$Tag_snp_bp) %>%
#     distinct() %>%
#     rename(CHROM = .data$Tag_snp_chrom, BP = .data$Tag_snp_bp)
#   
#   ld_table <- rbind(ld_table_subject,ld_table_tag, fill = TRUE)
#   
#   return(ld_table)
#   
# }
# 
# 
# # Function:- `window_unit_func` 
# # Context:- The user can either supply a numeric and granular window in numbers or bp -
# # or supply a kilobase window as characters "500kb".
# # Goal:- Take an input from the user and return the window value.
# 
# window_unit_func <- function(window_value) {
#   
#   # Check if the input is of type character
#   if (typeof(window_value) == 'character') {
#     
#     # Convert the user-supplied input to numeric and multiply it by 1000
#     kb_numeric = as.numeric(gsub("KB", "", window_value, ignore.case = TRUE))
#     
#     # Check if conversion was successful
#     if (is.na(kb_numeric)) {
#       stop("Error: Input could not be converted to numeric.")
#     }
#     
#     kb_window = kb_numeric * 1000
#     
#     return(kb_window)
#   }
#   
#   # If the input is already numeric, return it
#   if (typeof(window_value) == 'double' || typeof(window_value) == 'integer') {
#     return(window_value)
#   }
#   
#   stop("Error: Unsupported input type.")
# }
# 
# 
# has.annotations <- function(file_path) {
#   # Set up the path to the SnpSift jar file
#   jar_path <- system.file("java", "snpSift.jar", package = "panvaR")
#   
#   # Create temporary files for output and errors
#   output_file <- temp_file(prefix = "has.annotations_output_file")
#   error_file <- temp_file(prefix = "has.annotation_error_file")
#   
#   # Set up the arguments for the Java call - just try to get ANN field
#   binary_args <- c(
#     "-jar", jar_path, "extractFields", file_path, 
#     "CHROM", "POS", "ANN[*].EFFECT"
#   )
#   
#   # Try to execute the command and determine if annotations exist
#   result <- tryCatch(
#     {
#       # Execute the system call to run SnpSift
#       exit_code <- exec_wait(
#         "java",  # Main binary call is Java
#         args = binary_args,  # Pass the arguments
#         std_out = output_file,  # Capture standard output
#         std_err = error_file   # Capture error messages
#       )
#       
#       # Check if command executed successfully and produced output
#       if (exit_code == 0) {
#         # Read the first few lines of output to check for data
#         output_content <- readLines(output_file, n = 5)
#         
#         # If we have more than just a header line and the data isn't "."
#         # (which indicates missing annotation), then annotations exist
#         has_ann <- length(output_content) > 1 && 
#           !all(grepl("^\\.$", output_content[-1], perl = TRUE))
#         
#         return(has_ann)
#       } else {
#         # Read error message
#         error_msg <- readLines(error_file)
#         
#         # If error mentions missing ANN field, return FALSE
#         if (any(grepl("ANN.*not found", error_msg, ignore.case = TRUE))) {
#           return(FALSE)
#         } else {
#           # For other errors, print them but assume no annotations
#           print(paste("SnpSift error:", paste(error_msg, collapse = "\n")))
#           return(FALSE)
#         }
#       }
#     },
#     error = function(e) {
#       print(paste("Execution attempt produced error:", e$message))
#       return(FALSE)  # Return FALSE on error (assuming no annotations)
#     },
#     finally = {
#       # Clean up temporary files
#       if (file.exists(output_file)) file.remove(output_file)
#       if (file.exists(error_file)) file.remove(error_file)
#     }
#   )
#   
#   return(result)
# }
# 
# # Check if the gwas table has the required data fields -
# # if not then stop, if yes then return the right subset
# check_gwas_table <- function(current_gwas_table){
#   
#   needed_call_names <- c("CHROM","BP","Pvalues")
#   
#   # check if the supplied table is a path or a table
#   if(is.character(current_gwas_table)){
#     if (!file.exists(current_gwas_table)) {
#       stop("The GWAS table file path provided does not exist: ", current_gwas_table)
#     }
#     current_gwas_table <- data.table::fread(current_gwas_table)
#   } else if (!is.data.frame(current_gwas_table)) {
#     stop("Input to check_gwas_table must be a file path or a data.frame/data.table.")
#   }
#   
#   # are all the fields present in the current table
#   test_fields <- all(needed_call_names %in% colnames(current_gwas_table))
#   
#   if(!test_fields){
#     missing_cols <- setdiff(needed_call_names, colnames(current_gwas_table))
#     print(paste("You need to have these fields :-", paste(needed_call_names, collapse=", ")))
#     stop(paste("The GWAS table does not have the right data fields. Missing:", paste(missing_cols, collapse=", ")))
#   } else{
#     # If true select for the three fields that we need
#     current_gwas_table <- as.data.table(current_gwas_table) %>%
#       select(all_of(needed_call_names))
#     return(current_gwas_table)
#   }
# }
# 
# # A general function to convert vcf files to Plink
# vcf_to_bed <- function(vcf_file_path, output_prefix = NA) {
#   # make a prefix for output
#   if (is.na(output_prefix)) {
#     base_name_prefix <-
#       strsplit(vcf_file_path, ".", fixed = TRUE)[[1]][1]
#   } else {
#     base_name_prefix <- output_prefix
#   }
#   
#   if (!is.null(options()$plink_path)) {
#     binary_call <- options()$plink_path
#   } else {
#     binary_call <- "plink2"
#   }
#   
#   plink_args <- c("--vcf",
#                   vcf_file_path,
#                   "--make-bed",
#                   "--out",
#                   base_name_prefix)
#   
#   tryCatch({
#     try <- exec_wait(plink_call,
#                      args = plink_args,
#                      std_out = TRUE,
#                      std_err = FALSE)
#   }, error = function(e) {
#     # Custom error message
#     print(paste("Execution attempt produced error:-", e$message))
#     1 # Return 1 on error
#   })
#   
#   if (try == 0) {
#     print(paste("The Plink files are available in", base_name_prefix))
#   }
# }
# 
# # to call grep
# tabix_grep_call <- function(subset_file, gene_name) {
#   grep_call <- "grep"
#   
#   grep_output <- tempfile()
#   
#   the_args <- c("-E", # Extended regex patterns
#                 paste0("#|", gene_name), # Pattern to match
#                 subset_file)
#   
#   tryCatch({
#     try <- exec_wait(grep_call,
#                      args = the_args,
#                      std_out = grep_output,
#                      std_err = TRUE)
#   }, error = function(e) {
#     print(paste("There was an error with running grep on your input."))
#     print(paste("Execution attempt produced error:-", e$message))
#     stop()
#     1  # Return 1 on error
#   })
#   
#   return(grep_output)
# }
# 
# 
# # To call tabix
# tabix_subset_polymorphisms <- function(start, stop, chrom, gene_name, vcf_file) {
#   proper_tbi(vcf_file) # It is possible for the vcf file to not have a tbi file at this point.
#   
#   first_tabix_output <- tempfile()
#   
#   tabix_call <- "tabix"
#   
#   the_args = c("-h", vcf_file, paste0(chrom, ":", start, "-", stop))
#   
#   tryCatch({
#     try <- exec_wait(tabix_call,
#                      args = the_args,
#                      std_out = first_tabix_output,
#                      std_err = TRUE)
#   }, error = function(e) {
#     # Custom error message
#     print("There was an error with running grep on your input.")
#     print(paste("Execution attempt produced error:-", e$message))
#     stop()
#     1 # Return 1 on error
#   })
#   
#   grep_filtered <- tabix_grep_call(first_tabix_output, gene_name)
#   
#   return(grep_filtered)
#   
# }
# 
# # The tag snp should be supplied in the as
# # "chrom:bp" format
# tag_snp_check <- function(tag_snp) {
#   if (grepl(":", tag_snp)) {
#     return(TRUE)
#   } else {
#     message("The tag snp should be supplied in the as
#              \"chrom:bp\" format")
#     return(FALSE)
#   }
# }
# 
# plink2_bcftools_chroms_dictionary <- function(vcf_file_path,bim_file_path) {
#   
#   # Function to extract integers
#   extract_integers <- function(x) {
#     as.numeric(gsub("\\D", "", x))
#   }
#   
#   bcf_chroms <- bcf_chroms_func(vcf_file_path)
#   
#   bim_file_chroms <- fread(bim_file_path) %>% 
#     pull(1) %>% 
#     unique()
#   
#   bcf_chroms_table <- data.table(
#     chrom_names = bcf_chroms,
#     chrom_ints = sapply(bcf_chroms, extract_integers)
#   )
#   
#   chrom_test <- all(bim_file_chroms %in% bcf_chroms)
#   
#   if(!chrom_test){
#     plink2_chroms <- as.data.table(bim_file_chroms)
#     
#     vcf_plink2_dictionary_table <- left_join(
#       plink2_chroms,
#       bcf_chroms_table, 
#       by = c( "bim_file_chroms" = "chrom_ints" )
#     ) %>% rename(
#       "plink2" = "bim_file_chroms",
#       "vcf" = "chrom_names"
#     )
#     
#     return(vcf_plink2_dictionary_table)
#   } else {
#     return(NULL) # NULL is a substitute for TRUE here
#   }
#   
# }
# 
# apply_dict <- function(dictionary, table){
#   
#   # The function assumes that -
#   # The dictionary table has two fields - 
#   # 1. plink2, the chromosomes in plink2 -
#   # 2. vcf, the chromosomes in vcf
#   
#   # The table should have the following fields -
#   # CHROM
#   
#   fixed_table <- table %>% 
#     left_join(dictionary, by = c("CHROM"="plink2")) %>% # fix the CHROM's by matching them to plink2 names
#     select(-"CHROM") %>%
#     rename("CHROM" = "vcf") # rename vcf column to CHROMs
#   
#   return(fixed_table)
# }
# 
# # Improved tbi file check function that works in both Shiny and CLI modes
# proper_tbi <- function(vcf_file, shiny_mode = FALSE) {
#   # get the extension of the supplied vcf file
#   vcf_file_extenion <- extention_func(path_to_file = vcf_file)
#   
#   if (vcf_file_extenion != "gz") {
#     stop(
#       "The file you have supplied is not compressed using BGZIP and BCFtools cannot work with it. \n 
#       Please compress your file using BGZIP and try again.
#       "
#     )
#   }
#   
#   if (!file.exists(paste0(vcf_file, ".tbi"))) {
#     # Handle prompt differently based on environment
#     if (shiny_mode) {
#       # Return FALSE - in Shiny mode we'll trigger a modal dialog later
#       return(FALSE)
#     } else {
#       # CLI mode - use the terminal prompt
#       message("The .tbi file for the given vcf file does not exist - a .tbi file is needed for panvaR to function. Would you like to generate it? (y/n) ")
#       answer <- tolower(readline())
#       if (answer == "y") {
#         return(generate_tbi_file(vcf_file))
#       } else {
#         stop("The .tbi file for the given vcf file does not have a tabix index. Please generate it and try again.")
#       }
#     }
#   } else {
#     return(TRUE)
#   }
# }
# 
# # Extract the tbi generation functionality to a separate function
# generate_tbi_file <- function(vcf_file) {
#   message("Asking BCFtools to generate the .tbi file...")
#   
#   # First check if bcftools is available
#   if (!check_binary_available("bcftools")) {
#     stop("BCFtools is not available in your PATH. Please install it and try again.")
#   }
#   
#   binary_name <- "bcftools"
#   the_args <- c("index", "-t", vcf_file)
#   
#   tryCatch(
#     {
#       std_err_file <- tempfile()  # Create a temporary file to store standard error
#       result <- exec_wait(
#         binary_name,
#         args = the_args,
#         std_out = NULL,  # We don't need to capture standard output
#         std_err = std_err_file  # Redirect standard error to a file
#       )
#       
#       # After exec_wait, read the standard error from the file
#       std_err_content <- readLines(std_err_file)
#       if (length(std_err_content) > 0) {
#         message("Please read this error produced by tabix:")
#         message(paste(std_err_content, collapse = "\n"))
#       }
#       
#       # Clean up: Remove the temporary file
#       unlink(std_err_file)
#       
#       # Check if the result is an error (non-zero return code)
#       if (!is.numeric(result) || result != 0) {
#         stop("Failed to generate .tbi file. Please check the error message and try again.")
#       }
#       
#       return(TRUE)
#     },
#     error = function(e) {
#       # Custom error message
#       message("There was an error while generating the .tbi file.")
#       message(paste("Execution attempt produced this error:", conditionMessage(e)))
#       stop("Failed to generate .tbi file. Please check the error message and try again.")
#     }
#   )
# }
# 
# # Helper function to check if a binary is available
# check_binary_available <- function(binary_name) {
#   os_type <- Sys.info()["sysname"]
#   
#   if (os_type == "Windows") {
#     # On Windows, use where.exe to check if the binary exists
#     result <- suppressWarnings(system(paste("where", binary_name), intern = TRUE, ignore.stderr = TRUE))
#     return(!is.null(result) && length(result) > 0)
#   } else {
#     # On Unix-like systems, use which to check if the binary exists
#     result <- suppressWarnings(system(paste("which", binary_name), intern = TRUE, ignore.stderr = TRUE))
#     return(!is.null(result) && length(result) > 0)
#   }
# }
# 
# # Function to manage the temporary directory
# # Checks if the directory exists; if not, creates it. Deletes contents if specified.
# 
# temporary_directory <- function(delete_files = FALSE, working_directory = "panvar") {
#   
#   # Use provided working directory or default to "panvar"
#   dir_path <- ifelse(is.null(working_directory), "panvar", working_directory)
#   
#   # Delete the contents of the directory if it exists and delete_files is TRUE
#   if (dir.exists(dir_path) && delete_files) {
#     file.remove(list.files(dir_path, full.names = TRUE))
#   }
#   
#   # Create the directory if it does not exist
#   if (!dir.exists(dir_path)) {
#     dir.create(dir_path)
#   }
#   
#   return(dir_path)
# }
# 
# # Function to generate a temporary file path in the specified directory
# # Optionally creates the file
# 
# # `file_name`: Is the name of the intermediate file that is being given the -
# # the temporary name. For example, in a  run you will make many different -
# # temporary files, but more specifically different intermediate files serve -
# # different purposes. For example, you will have intermediate files 
# 
# temp_file <- function(create_file = FALSE, working_directory = "panvar", prefix = NULL) {
#   
#   # Use provided working directory or default to "panvar"
#   dir_path <- ifelse(is.null(working_directory), "panvar", working_directory)
#   
#   # Generate a random letters that you can use for the filename
#   random_letters <- paste0(sample(LETTERS, 5, replace = TRUE), collapse = "")
#   
#   if(!is.null(prefix)){
#     random_name <- paste0(prefix,"_",random_letters)
#   } else{
#     random_name <- paste0(random_letters)
#   }
#   
#   # Create the final path
#   final_path <- file.path(dir_path, random_name)
#   
#   # Create the file if specified
#   if (create_file) {
#     file.create(final_path)
#   }
#   
#   return(final_path)
# }
# 
# split_vcf_eff <- function(input_file = NULL, output_file = NULL) {
#   # Function to process a single line
#   process_line <- function(line) {
#     if (grepl("^#", line)) {
#       return(line)
#     }
#     
#     fields <- strsplit(trimws(line), "\t")[[1]]
#     info <- fields[INFO_FIELD_NUM]
#     info_fields <- strsplit(info, ";")[[1]]
#     
#     eff_field <- grep("^(EFF|ANN)=", info_fields, value = TRUE)
#     if (length(eff_field) == 0) {
#       return(line)
#     }
#     
#     field_name <- sub("=.*", "", eff_field)
#     effs <- strsplit(sub("^[^=]+=", "", eff_field), ",")[[1]]
#     
#     other_info <- info_fields[!grepl("^(EFF|ANN)=", info_fields)]
#     
#     pre <- paste(fields[1:(INFO_FIELD_NUM-1)], collapse = "\t")
#     post <- paste(fields[(INFO_FIELD_NUM+1):length(fields)], collapse = "\t")
#     
#     output_lines <- character(length(effs))
#     for (i in seq_along(effs)) {
#       new_info <- if (length(other_info) > 0) {
#         paste(c(other_info, paste0(field_name, "=", effs[i])), collapse = ";")
#       } else {
#         paste0(field_name, "=", effs[i])
#       }
#       output_lines[i] <- paste(pre, new_info, post, sep = "\t")
#     }
#     return(output_lines)
#   }
#   
#   INFO_FIELD_NUM <- 8  # R uses 1-based indexing
#   
#   # Main processing
#   if (is.null(input_file)) {
#     con_in <- stdin()
#   } else {
#     con_in <- file(input_file, "r")
#   }
#   
#   # Handle output
#   if (is.null(output_file)) {
#     output_file <- temp_file(prefix = "oneperline")
#   }
#   con_out <- file(output_file, "w")
#   
#   # Main processing
#   tryCatch({
#     while (TRUE) {
#       line <- readLines(con_in, n = 1)
#       if (length(line) == 0) break
#       processed_lines <- process_line(line)
#       writeLines(processed_lines, con_out)
#     }
#   }, finally = {
#     if (!is.null(input_file)) close(con_in)
#     close(con_out)
#   })
#   
#   # Return the output file name
#   return(output_file)
# }
# 
# bcf_chroms_func <- function(vcf_file_path){
#   
#   # make a temporary file to catch the standard output
#   
#   std_output_name <- temp_file(create_file=TRUE, prefix = "bcf_chroms_function") # Defined in the file R/temporary.r
#   # The shell call that I want to run is:-
#   # `bcftools index -s file.vcf.gz``
#   
#   binary_call <- "bcftools"
#   
#   binary_args <- c(
#     "index","-s",
#     vcf_file_path
#   )
#   
#   tryCatch(
#     {
#       try <- exec_wait(
#         binary_call,
#         args = binary_args,
#         std_out = std_output_name,
#         std_err = FALSE
#       )
#     },
#     error = function(e){
#       # Custom error message
#       print(paste("Execution attempt produced error:-", e$message))
#       1 # Return 1 on error
#     }
#   )
#   
#   if (try == 0){
#     
#     final_index_table = fread(std_output_name)
#     
#     # get the column names
#     bcf_chroms <- final_index_table %>%
#       pull(1)
#     
#     return(bcf_chroms)
#   } else {
#     
#     print("There were errors when calculating index for this set of inputs.")
#     print("Please read the error message and re-try")
#     return(readLines(std_output_name))
#   }
#   
# }
# 
execute_snpsift <- function(file_path) {
  message("~~~~~~~~~~~~~~~ Extracting SNP impacts ~~~~~~~~~~~~~~~")
  # Create a temporary file for the output
  snp_eff_table <- temp_file(prefix = "snp_eff_table")

  # Set up the path to the SnpSift jar file
  jar_path <- system.file("java", "snpSift.jar", package = "panvaR")

  # Set up the arguments for the Java call
  binary_args <- c(
    "-jar", jar_path, "extractFields", file_path,
    "CHROM", "POS", "ANN[*].FEATUREID",
    "REF", "ALT", "ANN[*].EFFECT",
    "ANN[*].AA", "ANN[*].IMPACT"
  )

  # Try to execute the command and catch any errors
  tryCatch(
    {
      error_message <- temp_file(prefix = "snpsift_error")  # Create a temporary file for error messages

      # Execute the system call to run SnpSift
      exec_wait(
        "java",  # Main binary call is Java
        args = binary_args,  # Pass the arguments
        std_out = snp_eff_table,  # Capture standard output in the table
        std_err = error_message  # Capture error messages
      )

      # snp_eff_table is the path -
      # let's return a data.table as well
      snp_eff_datatable <- fread(snp_eff_table)

      # Change the names of the table
      colnames(snp_eff_datatable) <- c(
        "CHROM",
        "BP",
        "GENE",
        "REF",
        "ALT",
        "EFFECT",
        "AA",
        "IMPACT"
      )

      # Make a list and return it

      return(
        list(
          path = snp_eff_table,
          table = snp_eff_datatable
        )
      )

    },
    error = function(e) {
      print(paste("Execution attempt produced error:", e$message))
      return(1)  # Return 1 on error
    }
  )
}
