# panvar_gui

#' @export

panvar_gui <- function(display.mode = "auto", launch.browser = TRUE, ...) {
  gui_code <- system.file("shiny", "panvar_shiny_main.r", package = "panvaR")
  shiny::runApp(gui_code, display.mode = display.mode, launch.browser = launch.browser, ...)
}