# Holds code that will make calling shiny app easier.

#' @export
panvar_gui <- function(display.mode = "normal", launch.browser = TRUE, ...) {
    gui_code <- system.file("shiny", "frontend.R", package = "panvaR")
    shiny::runApp(gui_code, display.mode = display.mode, launch.browser = launch.browser, ...)
}
