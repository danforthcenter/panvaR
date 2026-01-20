# Improved tbi file check function that works in both Shiny and CLI modes
proper_tbi <- function(vcf_file, shiny_mode = FALSE) {
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
    # Handle prompt differently based on environment
    if (shiny_mode) {
      # Return FALSE - in Shiny mode we'll trigger a modal dialog later
      return(FALSE)
    } else {
      # CLI mode - use the terminal prompt
      message("The .tbi file for the given vcf file does not exist - a .tbi file is needed for panvaR to function. Would you like to generate it? (y/n) ")
      answer <- tolower(readline())
      if (answer == "y") {
        return(generate_tbi_file(vcf_file))
      } else {
        stop("The .tbi file for the given vcf file does not have a tabix index. Please generate it and try again.")
      }
    }
  } else {
    return(TRUE)
  }
}

# Extract the tbi generation functionality to a separate function
generate_tbi_file <- function(vcf_file) {
  message("Asking BCFtools to generate the .tbi file...")
  
  # First check if bcftools is available
  if (!check_binary_available("bcftools")) {
    stop("BCFtools is not available in your PATH. Please install it and try again.")
  }
  
  binary_name <- "bcftools"
  the_args <- c("index", "-t", vcf_file)
  
  tryCatch(
    {
      std_err_file <- tempfile()  # Create a temporary file to store standard error
      result <- exec_wait(
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
      
      # Check if the result is an error (non-zero return code)
      if (!is.numeric(result) || result != 0) {
        stop("Failed to generate .tbi file. Please check the error message and try again.")
      }
      
      return(TRUE)
    },
    error = function(e) {
      # Custom error message
      message("There was an error while generating the .tbi file.")
      message(paste("Execution attempt produced this error:", conditionMessage(e)))
      stop("Failed to generate .tbi file. Please check the error message and try again.")
    }
  )
}

# Helper function to check if a binary is available
check_binary_available <- function(binary_name) {
  os_type <- Sys.info()["sysname"]
  
  if (os_type == "Windows") {
    # On Windows, use where.exe to check if the binary exists
    result <- suppressWarnings(system(paste("where", binary_name), intern = TRUE, ignore.stderr = TRUE))
    return(!is.null(result) && length(result) > 0)
  } else {
    # On Unix-like systems, use which to check if the binary exists
    result <- suppressWarnings(system(paste("which", binary_name), intern = TRUE, ignore.stderr = TRUE))
    return(!is.null(result) && length(result) > 0)
  }
}