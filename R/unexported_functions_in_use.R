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
  as.numeric(stringr::str_extract(x, "(^.*)-(.*)", group = 1))
}