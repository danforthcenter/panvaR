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

