#' sets global option "panvar_prefix" for output directory of intermediate panvar files. If not supplied will be sent to output of 'tempdir()'.
#' 
#' 
#' @param prefix character, a prefix to designate files from individual runs of the program. 
#'
#' @export
#'
#' @examples
#' # set_panvar_prefix("Phenotype1")
set_panvar_prefix <- function(prefix){
  options(panvar_prefix=prefix)
}