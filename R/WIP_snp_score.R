#' Title
#'
#' @param gwas.res 
#' @param tag.snp 
#' @param ld.list 
#' @param window 
#' @param score.cols 
#' @param score.dirs 
#' @param score.weights 
#'
#' @returns
#' @export
#'
#' @examples
make_snp_score <- function(gwas.res,
                           tag.snp,
                           ld.list,
                           window,
                           default.score.cols,
                           custom.score.cols,
                           custom.score.dirs,
                           score.weights){
  
  # default cols:
  # - LD
  # - distance
  # - P-value
  gwas.sub <- gwas.res %>%
    as.data.frame() %>% 
    mutate(marker.ID = paste(.data$CHR, .data$POS, sep = "-")) %>%
    left_join(ld.list$table, by = "marker.ID") %>%
    filter(!is.na(.data$R2)) %>% 
    # rename this here to work with special case "LD" given to point variable
    rename("LD" = "R2")
  
  
}