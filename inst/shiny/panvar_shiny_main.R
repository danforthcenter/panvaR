# app.R
library(shiny)
library(shinyFiles)
library(tidyverse)
library(data.table)
library(plotly)
library(DT)
library(shinyBS)
library(shinyjs)

# Dynamically get the path to the Input_dashboard_panvar.R and Output_dashboard_panvar.R files
input_dashboard_path <- system.file("shiny", "Input_dashboard_panvar.R", package = "panvaR")
output_dashboard_path <- system.file("shiny", "Output_dashboard_panvar.R", package = "panvaR")

source(input_dashboard_path)
source(output_dashboard_path)

ui <- fluidPage(
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
  ),
  useShinyjs(),
  tabsetPanel(
    id = "mainTabs",
    tabPanel("Input", input_dashboard_UI("module1")),
    tabPanel("Output", output_dashboard_UI("module2"))
  )
)

server <- function(input, output, session) {

  # We do not want the shiny code to just automatically run in the `app.R` dir -
  # .So, let's set it to work in the working directory of the R interpreter -
  # that initiates the GUI

  setwd(getwd())
  
  # Create shared reactive values
  shared <- reactiveValues(
    analysis_results = NULL
  )
  
  # Initialize modules
  input_results <- input_dashboard_Server("module1", shared)
  output_dashboard_Server("module2", shared)
  
  # Observe successful analysis completion and switch tabs
  observeEvent(shared$analysis_results, {
    if (!is.null(shared$analysis_results)) {
      updateTabsetPanel(session, "mainTabs", selected = "Output")
    }
  })
}

shinyApp(ui, server)
