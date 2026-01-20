#' Set path for output files
#'
#' @param out_dir character, path to directory
#'
#' @returns
#' sets global option "panvar_outdir" for output directory of intermediate panvar files. If not supplied will be sent to output of 'tempdir()'.
#' 
#' 
#' @export
#'
#' @examples
#' # set_out_dir("/home/user/output_directory")
set_out_dir <- function(out_dir){
  options(panvar_outdir=out_dir)
}