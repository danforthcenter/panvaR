#' Extract position from marker.ID in the form "CHR-POS"
#'
#' @param marker.ID character, marker.ID with chromosome and position separated by some seperator
#' @param sep character, character that separates chromosome and position in marker.ID
#'
#' @returns
#' position for this marker as a number
#' @export
#'
#' @examples
#' my.marker <- "Chr01-123"
#' get_bp_from_id(my.marker, sep = "-")
get_bp_from_id <- function(marker.ID, sep = "-"){
  as.numeric(stringr::str_extract(marker.ID, paste0(sep, "(.*)"), group = 1))
}


#' Extract chromsome from marker.ID in the form "CHR-POS"
#'
#' @param marker.ID character, marker.ID with chromosome and position separated by some seperator
#' @param sep character, character that separates chromosome and position in marker.ID
#' @param tonumeric boolean, should the function attempt to convert the chromosome to numeric?
#'
#' @returns
#' chromsome for this marker as a number if you want
#' @export
#'
#' @examples
#' my.marker <- "Chr01-123"
#' get_chrom_from_id(my.marker, tonumeric = T)
#' get_chrom_from_id(my.marker, tonumeric = F)
get_chrom_from_id <- function(marker.ID, sep = "-", tonumeric = T){
  if(tonumeric){
    out <- stringr::str_extract(marker.ID, paste0("(^.*)", sep, "(.*)"), group = 1)
    
    if(any(letters_only(out))){
      warning("tonumeric selected but chromosome is only characters.")
    } else if(all(numbers_only(out))){
      out <- as.numeric(out)
    } else {
      out <- as.numeric(apply(str_extract_all(out, "\\d", simplify = T), 1, paste, collapse = ""))
    }
  
  } else {
    out <- stringr::str_extract(marker.ID, paste0("(^.*)", sep, "(.*)"), group = 1)
  }
  return(out)
}


# ------------------------------------------------------------------------\
# unexported --------
# ------------------------------------------------------------------------\

# Check that it doesn't match any non-letter
letters_only <- function(x) !grepl("[^A-Za-z]", x)

# Check that it doesn't match any non-number
numbers_only <- function(x) !grepl("\\D", x)

