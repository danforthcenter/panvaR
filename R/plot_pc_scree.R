#' Calculate and plot pcs
#'
#' @param bedfile.path character, path to bedfile (include .bed extension)
#' @param num.pcs numeric, total number of pcs to plot
#' @param cumulative boolean, if T will plot the cumulative proportion of variance
#'
#' @returns
#' GGplot of the given plot. 
#' 
#' The output is sensitive to the number of PC's indicated. The program used to calculate
#' PCs returns singular values only for components indicated. So the proportion can only calculate using 
#' the variance contained in these components. Consequently, shouldn't set the num.pcs too low, this could lead to strange results. 
#' @export
#'
#' @examples
#' # work in progress
plot_pc_scree <- function(bedfile.path,
                          num.pcs,
                          cumulative = F){
  
  bedfile.path <- normalizePath(bedfile.path)
  bedfile.path <- read.pcadapt(bedfile.path,
                               type = "bed")
  x <- pcadapt(input = bedfile.path, K = num.pcs) 
  
  plot.df <- data.frame(PC = 1:num.pcs, 
                        prop.var = x$singular.values^2 / sum(x$singular.values ^ 2), 
                        eigenvalues = x$singular.values^2) %>% 
    mutate(cumulative = cumsum(prop.var))
  
  if(!cumulative){
    p <-
      ggplot(aes(x = .data$PC, y = .data$prop.var), data = plot.df) +
      geom_point() +
      geom_line() +
      theme_bw(13) +
      labs(title = "Scree Plot",
           y = "Proportion of Variance Explained (approximate)")
  } else {
    p <- 
      ggplot(aes(x = .data$PC, y = .data$cumulative), data = plot.df) +
      geom_point() +
      geom_line() +
      theme_bw(13) +
      labs(title = "Scree Plot - Cumulative",
           y = "Cumulative Proportion of Variance Explained (approximate)")
    
  }
  return(p)
}