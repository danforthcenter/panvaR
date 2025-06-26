#' vcf_to_plink
#'
#' @description 
#' This is a function that takes a vcf file and converts it to Plink format.
#'
#' @param vcf_file_path Path to the input VCF file.
#' @param output_prefix The path where the Plink file should be written.
#' Defaults to writing the output file in a `tempdir`. If you want the file to be kept around longer supply an accessible path.
#' @param auto_generate_tbi (Optional) If TRUE, automatically generates the .tbi index file for the VCF if it is missing. Defaults to FALSE.
#' 
#' @return A list containing the paths to the new bed and bim files.
#'
#' @examples
#' # \dontrun{
#' # vcf_to_plink2("path/to/your_file.vcf.gz", auto_generate_tbi = TRUE)
#' # }
#' 
#' @import tidyverse
#' @import data.table
#' @import parallel
#' @importFrom bigstatsr big_randomSVD
#' @import modelr
#' 
#' @export

# A general function to convert vcf files to Plink

vcf_to_plink2 <- function(vcf_file_path, output_prefix = NA, auto_generate_tbi = FALSE){
  
  # Check for and handle missing .tbi file before proceeding
  proper_tbi(vcf_file_path, auto_generate_tbi = auto_generate_tbi)
  
  # make a prefix for output
  if(is.na(output_prefix)){
    
    # tempoarary directory to store the files
    
    plink2_tempdir <- temporary_directory()
    
    base_name_prefix <- 
      base_name_func(vcf_file_path)
    
    output_path = paste0(plink2_tempdir,"/",base_name_prefix)
    
  } else {
    output_path <- output_prefix
  }
  # TODO :- If the user supplies a name copy the bed file from scratch to  
  # create a temporary directory 
  
  
  # check if bed and bim files already exist
  bim_file_path = paste0(output_path,".bim")
  bed_file_path = paste0(output_path,".bed")
  
  if(file.exists(bim_file_path) && file.exists(bed_file_path)){
    # return the path to the bed file
    print("It looks like your bed and bim files already exist. This was done with heuristics and regex.")
    
    print("The heuristics could be wrong, if that's the case please delete them and try again.")
    
    return(list(bed = bed_file_path,bim = bim_file_path))
  }
  binary_call <- "plink2"
  
  binary_args <- c(
    "--allow-extra-chr",
    "--vcf",
    vcf_file_path, "--make-bed",
    "--set-all-var-ids", "Chr_@_BP_#",
    "--out", output_path
  )
  
  # Rijan: Given how legible Plink2's error messages are
  # Rijan: I think a simple STDOUT catch is good enough
  
  # Rijan: By the looks of it Plink2 can handle .gz and .vcf paths just fine
  
  tryCatch(
    {
      error_message <- tempfile()
      try <- exec_wait(
        binary_call,
        args = binary_args,
        std_out = TRUE,
        std_err = error_message
      )
    },
    error = function(e){
      # Custom error message
      print(paste("Execution attempt produced error:-", e$message))
      1 # Return 1 on error
    }
  )
  
  
  
  if(try == 0){
    
    bed_file_path = paste0(output_path,".bed")
    
    bim_file_path = paste0(output_path,".bim")
    
    print(
      paste("The Plink files are available in",bed_file_path)
    )
    return(list(bed = bed_file_path,bim = bim_file_path))
  } else{
    
    print("Your VCF file could not be converted to Plink2's native format.")
    print("Please read the error message and re-try")
    return(readLines(error_message))
  }
}