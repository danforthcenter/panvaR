#' Set panvar output file prefix
#' sets global option "panvar_prefix". Appended to outputs for organization and delineating different runs of the program. 
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