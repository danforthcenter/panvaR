
#' panvar_gui
#' Run the panvaR shiny implementation. Will launch in a default browser or can navigate to IP address shown. 
#' 
#' @param display.mode display mode of shiny app, see documentation for `shiny::runApp()`
#'
#' @param launch.browser should browser be launched automatically or not
#' @param ... arguments passed to `shiny::runApp()`
#'
#' @export

panvar_gui <- function(display.mode = "auto", launch.browser = TRUE, ...) {
  gui_code <- system.file("shiny", "panvar_shiny_main.r", package = "panvaR")
  shiny::runApp(gui_code, display.mode = display.mode, launch.browser = launch.browser, ...)
}