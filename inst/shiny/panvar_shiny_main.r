# panvar_shiny_main.r
library(shiny)
library(shinyFiles)
library(tidyverse)
library(data.table)
library(plotly)
library(DT)
library(shinyBS)
library(shinyjs) # Ensure shinyjs is loaded

# --- Source all module files ---
# A more robust helper function to source module files
source_module <- function(file_name) {
  # Path for development environment (e.g., running app from project root)
  dev_path <- file.path("inst", "shiny", file_name)
  
  # Path for installed package environment
  pkg_path <- system.file("shiny", file_name, package = "panvaR")
  
  if (file.exists(dev_path)) {
    # Use development path if it exists
    source(dev_path)
    message("Sourced from dev path: ", dev_path)
  } else if (nzchar(pkg_path)) { # nzchar checks if the string is non-empty
    # Use package path if found
    source(pkg_path)
    message("Sourced from package path: ", pkg_path)
  } else {
    # If neither path works, stop with an informative error
    stop(paste("Could not find module file:", file_name,
               ". Looked in:", dev_path, "and in the installed package."))
  }
}


source_module("Input_dashboard_panvar.r")
source_module("Output_dashboard_panvar.r")
source_module("Gwas_input_dashboard_panvar.r")


ui <- fluidPage(
  tags$head(
    # Link to styles.css if you have one
    # tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
  ),
  useShinyjs(), # Initialize shinyjs
  div(
    class = "app-header",
    style = "background-color: #f8f9fa; padding: 10px 15px; border-bottom: 1px solid #dee2e6; margin-bottom: 20px; display: flex; justify-content: space-between; align-items: center;",
    div(
      style = "display: flex; align-items: center;",
      tags$h2("PanvaR Analysis Tool", style = "margin: 0;")
    ),
    div(
      actionButton("restart_app", "Restart Analysis", class = "btn-success", style = "margin-right: 10px;"),
      actionButton("exit_app", "Exit Application", class = "btn-danger")
    )
  ),
  tabsetPanel(
    id = "mainTabs",
    tabPanel("De Novo Analysis", input_dashboard_UI("module1")),
    tabPanel("Analysis from GWAS", Gwas_input_dashboard_UI("module3")), # New Tab
    tabPanel("Results", output_dashboard_UI("module2"))
  )
)

server <- function(input, output, session) {
  
  # Exit and Restart logic
  observeEvent(input$restart_app, {
    # Reset inputs in both input modules
    shinyjs::reset("module1-R2_threshold")
    shinyjs::reset("module1-tagSnps")
    # Add resets for other inputs in module1 as needed
    
    shinyjs::reset("module3-gwas_table_path") 
    shinyjs::reset("module3-vcf_file_path")
    shinyjs::reset("module3-tag_snps")
    
    shared$analysis_results <- NULL
    updateTabsetPanel(session, "mainTabs", selected = "De Novo Analysis")
    showNotification("Inputs and results have been cleared.", type = "message")
  })
  
  observeEvent(input$exit_app, {
    showModal(modalDialog(
      title = "Confirm Exit",
      "Are you sure you want to exit the application?",
      footer = tagList(
        actionButton("confirm_exit", "Yes, Exit", class = "btn-danger"),
        modalButton("Cancel")
      ),
      easyClose = TRUE
    ))
  })
  
  observeEvent(input$confirm_exit, {
    stopApp()
  })
  
  # Create shared reactive values
  shared <- reactiveValues(
    analysis_results = NULL
  )
  
  # Initialize all three modules
  input_dashboard_Server("module1", shared) # De Novo Analysis module
  output_dashboard_Server("module2", shared) # Results module
  Gwas_input_dashboard_Server("module3", shared) # Analysis from GWAS module
  
  # Observe successful analysis completion and switch to the Results tab
  observeEvent(shared$analysis_results, {
    if (!is.null(shared$analysis_results)) {
      updateTabsetPanel(session, "mainTabs", selected = "Results")
    }
  }, ignoreNULL = TRUE)
  
}

shinyApp(ui, server)