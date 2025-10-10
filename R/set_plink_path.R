#' Set path to plink2 executable
#'
#' @param plink_path character, path to plink2 executable
#'
#' @returns
#' sets global option that for path to global executable for plink2. If plink2 is on path and option is not supplied, will use the executable on the path. 
#' 
#' plink V2.00a6 minimum recommended 
#' @export
#'
#' @examples
set_plink_path <- function(plink_path){
  options(plink_path=plink_path)
}