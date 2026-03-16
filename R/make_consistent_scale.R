#' Generate discrete scale
#' 
#' Make a ggplot discrete scale that has consistent correspondence between aesthetics and variables
#'
#' @param values vector, values to supply to given scale (e.g. colors or shapes)
#' @param vars character, variables to link to values. Will have one to one correspondence to values based on order
#' @param type character, either "fill", "color" or "shape"
#' @param show.example boolean, print a plot to show what values are tied to which variables
#' @param name optional, string to name scale, will appear as legend heading when used in a ggplot
#'
#' @returns ggplot scale object for use in downstream plotting
#' @export
#'
#' @examples
#' # work in progress
#' 
#' # fill
#' my.scale <- 
#' make_consistent_scale(
#' values =  c("red", "blue", "gold"),
#' vars = c("A", "B", "C"),
#' show.example = TRUE)
#' 
#' # shape 
#' my.scale <- 
#' make_consistent_scale(
#' values =  c(21, 22, 23),
#' vars = c("A", "B", "C"),
#' type = "shape",
#' show.example = TRUE)
make_consistent_scale <- function(values,
                                  vars,
                                  type = c("fill", "color", "shape"),
                                  show.example = F,
                                  name = NULL) {
  if(length(values) != length(vars)){
    stop("Values and variables should be same length.")
  }
  
  if(is.null(name)){
    name <- ggplot2::waiver()
  }
  
  names(values) <- levels(factor(vars, levels = vars))
  type <- match.arg(type)
  if (type == "fill") {
    out <- scale_fill_manual(name = name, values = values)
    if(show.example){
      khroma::plot_scheme(values, names = T)
    }
  } else if (type == "color") {
    out <- scale_color_manual(name = name, values = values)
    if(show.example){
      khroma::plot_scheme(values, names = T)
    }
  } else if (type == "shape") {
    out <- scale_shape_manual(name = name, values = values)
    if(show.example){
      x <- data.frame(shape = values, label = names(values), x.coord = seq(1, length(values), by = 1))
      # p <-
      #   ggplot(aes(x = .data$x.coord, y = 1), data = x) +
      #   geom_point(aes(shape = factor(.data$shape, levels=.data$shape)), show.legend =  F, size = 20, fill = "orange") +
      #   scale_shape_manual(values = x$shape) +
      #   geom_text(aes(label = .data$label, y = 2), size = 8) +
      #   theme_void() +
      #   theme(plot.margin = unit(rep(2,4), "cm")) +
      #   coord_cartesian(clip = "off") +
      #   ylim(0,3)
      p <-
        ggplot(aes(x = .data$x.coord, y = 1), data = x) +
        geom_point(aes(shape = .data$label), show.legend =  F, size = 20, fill = "orange") +
        # scale_shape_manual(values = x$shape) +
        out +
        geom_text(aes(label = .data$label, y = 2), size = 8) +
        theme_void() +
        theme(plot.margin = unit(rep(2,4), "cm")) +
        coord_cartesian(clip = "off") +
        ylim(0,3)
      print(p)
    }
  } else {
    stop("Unsupported type argument")
  }
  return(out)
}