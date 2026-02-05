#' Make scores
#'
#' @param input.df data.frame, table with only snps that are to be included.
#' Must have "marker.ID" column
#' @param cols character, vector of column names indicating which variables to included in the score
#' @param directions numeric, a vector indicating which direction is to be considered more indicative
#' of an association. 1 indicates higher is better, -1 indicates lower is better. The order should correspond 
#' with the order in cols. 
#' @param weights numeric, a vector indicating weights for the variables. These must add up to 1. 
#'
#' @returns
#' A table with a score for each marker based on aggregating the variables indicated. 
#' 
#' 1) Variables are first made negative if indicated by `directions` vector. 
#' 2) Min/max normalization is applied to put all variables on same scale. 
#' 3) A weighted mean based on the `weights` vector is taken of the normalized variables and output as the final score. 
#' 
#' @export
#'
#' @examples
#' # work in progress
make_scores <- function(input.df, cols, directions, weights = NULL){
  # give equal weight if none provided 
  if(is.null(weights)){
    weights <- rep.int(1, length(cols)) / length(cols)
  }
  # make sure weights are valid
  if(sum(weights) != 1){
    stop("Weights must sum to 1.")
  }
  
  # correspond variables, weights and directions
  # based on order of vectors
  dirs.df <- data.frame(var = cols, direction = directions, weight = weights)
  
  # do the calculations
  out <- input.df %>% 
    select("marker.ID", all_of(cols)) %>% 
    tidyr::pivot_longer(cols = -.data$marker.ID, names_to = "var", values_to = "value") %>% 
    left_join(dirs.df, by = "var") %>% 
    rowwise() %>% 
    mutate(value.direction = .data$value * .data$direction) %>% 
    group_by(.data$var) %>% 
    mutate(value.normalized = minmaxnormal(.data$value.direction)) %>% 
    group_by(.data$marker.ID) %>% 
    mutate(snp.score = weighted.mean(.data$value.normalized, .data$weight)) %>%
    # stop pipe here to test weighting
    ungroup() %>% 
    select("marker.ID", "snp.score") %>% 
    distinct()
  
  return(out)
}

