# unexported functions I am actually using


# ------------------------------------------------------------------------\
# re-used functions --------
# ------------------------------------------------------------------------\

# re-used from old file
# not copying these over in order to retain functionality
# will iterate in the current file then copy over when I fully transfer over to using this file

# temporary_directory 

# ------------------------------------------------------------------------\
# new functions --------
# ------------------------------------------------------------------------\

# these functions are new functions that I have created that will be part of the package going forward 


#' get an input that has a certain prefix and suffix
#' will match exactly ^prefixsuffix$ if excluding the ".*" in the middle
#' might have issues with unescaped "."'s 
#'
#' @param inputs.dir directory to look in
#' @param in.prefix file prefix
#' @param suffix file suffix
#'
#' @returns
#' full path to file
#' 
#' @keywords internal
get_an_input <- function(inputs.dir, in.prefix, suffix){
  # list.files(path = inputs.dir, pattern = paste0("^", in.prefix, suffix, "$"), full.names = T)
  list.files(path = inputs.dir, pattern = paste0("^", in.prefix, ".*", suffix, "$"), full.names = T)
}



#' get_geno_filetype 
#' Get the format of a genotype file
#'
#' @param genotype.path path to genotype file
#'
#' @returns
#' either "vcf", "bed", or "unsupported_filetype"
#'
#' @keywords internal
get_geno_filetype <- function(genotype.path){
  f <- basename(genotype.path)
  if(str_detect(f, "bed$")){
    return("bed")
  } else if(str_detect(f, "vcf.gz$") | str_detect(f, "vcf$")){
    return("vcf")
  } else {
    return("unsupported_filetype")
  }
}



#' use plink to calculate ld
#'
#' @param plink.path path to plink2 executable
#' @param snp.name name of snp
#' @param window total window size in KB, all variants within .5 * window are calculated
#' @param bedfile bedfile prefix, no .bed
#' @param in.dir path to location of bed file
#' @param out.dir path to output temp file 
#'
#' @returns
#' ld values in file named 'ld_out_temp.vcor' 
#' @export
#'
#' @examples
#' # work in progress
make_ld <- function(plink.path,
                    snp.name,
                    window,
                    bedfile,
                    in.dir,
                    out.dir) {
  system(
    paste0(
      plink.path,
      " --silent --bfile ",
      in.dir,
      bedfile,
      " --r2-unphased --ld-snp ",
      snp.name,
      " --ld-window-kb ",
      window,
      " --ld-window 99999 --ld-window-r2 0 --out ",
      out.dir, 
      "/ld_out_temp"
    )
  )
}


# extract position from marker.ID in the form "CHR-POS"
get_bp_from_id <- function(marker.ID){
  as.numeric(stringr::str_extract(marker.ID, "-(.*)", group = 1))
}

# extract chromsome from marker.ID in the form "CHR-POS"
get_chrom_from_id <- function(marker.ID){
  as.numeric(stringr::str_extract(marker.ID, "(^.*)-(.*)", group = 1))
}


execute_snpsift2 <- function(file_path) {
  message("~~~~~~~~~~~~~~~ Extracting SNP impacts ~~~~~~~~~~~~~~~")
  file_path <- path.expand(file_path)

  # Create a temporary file for the output
  snp_eff_table <- temp_file(prefix = "snp_eff_table")

  # Set up the path to the SnpSift jar file
  jar_path <- system.file("java", "snpSift.jar", package = "panvaR")

  # Set up the arguments for the Java call
  binary_args <- c(
    "-jar", jar_path, "extractFields", file_path,
    "-e", ".", "-s", ",",
    "CHROM", "POS", "ANN[*].FEATUREID", "ANN[*].GENE",
    "REF", "ALT", "ANN[*].EFFECT", "ANN[*].BIOTYPE",
    "ANN[*].AA", "ANN[*].IMPACT"
  )

  # Try to execute the command and catch any errors
  try <- tryCatch(
    {
      error_message <- temp_file(prefix = "snpsift_error")  # Create a temporary file for error messages

      # Execute the system call to run SnpSift
      exec_wait(
        "java",  # Main binary call is Java
        args = binary_args,  # Pass the arguments
        std_out = snp_eff_table,  # Capture standard output in the table
        std_err = error_message  # Capture error messages
      )


    },
    error = function(e) {
      print(paste("Execution attempt produced error:", e$message))
      return(1)  # Return 1 on error
    }
  )
  if(try == 1){
    return(readLines(error_message))
  } else {
    return(snp_eff_table)
  }
}

temp_file <- function(create_file = FALSE, working_directory = NULL, prefix = NULL) {
  
  if(!is.null(options()$panvar_outdir)){
    dir_path <- options()$panvar_outdir
  } else if(is.null(working_directory)){
    dir_path <- tempdir()
  } else {
    dir_path <- working_directory
  }
  # dir_path <- ifelse(is.null(working_directory), tempdir(), working_directory)
  
  # Generate a random letters that you can use for the filename
  random_letters <- paste0(sample(LETTERS, 5, replace = TRUE), collapse = "")
  
  if(!is.null(prefix)){
    random_name <- paste0(prefix,"_",random_letters)
  } else{
    random_name <- paste0(random_letters)
  }
  
  # Create the final path
  final_path <- file.path(dir_path, random_name)
  
  # Create the file if specified
  if (create_file) {
    file.create(final_path)
  }
  
  return(final_path)
}
