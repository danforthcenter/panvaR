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
#' # These are some variables you'll be working with.
#' # We want to retain the correspondence between a color/shape in plots
#' # and these variables. 
#' my.vars <- c("A", "B", "C")
#' 
#' # fill
#' my.scale <- 
#' make_consistent_scale(
#' values =  c("red", "blue", "gold"),
#' vars = my.vars,
#' show.example = TRUE)
#' 
#' # shape 
#' 
#' my.scale <- 
#' make_consistent_scale(
#' values =  c(21, 22, 23),
#' vars = my.vars,
#' type = "shape",
#' show.example = TRUE)
#' 
#' # will retain correspondence when a level is dropped
#' plot.df <- data.frame(
#'   vars = my.vars[-2],
#'   value = c(2,1))
#'   
#' ggplot(data = plot.df, aes(x = value, y = 1)) +
#'   geom_point(aes(shape = vars),  fill = "orange", size = 10) +
#'   geom_text(aes(y = 1.1, label = vars)) +
#'   my.scale +
#'   ggplot2::theme_void() +
#'   ylim(.75, 1.25) +
#'   xlim(.5, 2.5) +
#'   theme(legend.position = "bottom")
#' 
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