# panvar_shiny_main.r
library(shiny)
library(shinyFiles)
library(tidyverse)
library(data.table)
library(plotly)
library(DT)
library(shinyBS)
library(shinyjs) # Ensure shinyjs is loaded

# Dynamically get the path to the Input_dashboard_panvar.R and Output_dashboard_panvar.R files
# Adjust these paths if they are not correctly resolved by system.file in your environment
# For example, if running directly from the script's location:
# input_dashboard_path <- file.path("..", "shiny", "Input_dashboard_panvar.r") # Example adjustment
# output_dashboard_path <- file.path("..", "shiny", "Output_dashboard_panvar.r") # Example adjustment
input_dashboard_path <- tryCatch(system.file("shiny", "Input_dashboard_panvar.r", package = "panvaR"), error = function(e) "Input_dashboard_panvar.r") # Fallback if package structure isn't standard
output_dashboard_path <- tryCatch(system.file("shiny", "Output_dashboard_panvar.r", package = "panvaR"), error = function(e) "Output_dashboard_panvar.r") # Fallback

# Check if files exist before sourcing
if (file.exists(input_dashboard_path)) {
  source(input_dashboard_path)
} else {
  stop("Could not find Input_dashboard_panvar.r at path: ", input_dashboard_path)
}
if (file.exists(output_dashboard_path)) {
  source(output_dashboard_path)
} else {
  stop("Could not find Output_dashboard_panvar.r at path: ", output_dashboard_path)
}


ui <- fluidPage(
  tags$head(
    # Link to styles.css if you have one
    # tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
  ),
  useShinyjs(), # Initialize shinyjs
  # Add heading with exit and restart buttons
  div(
    class = "app-header",
    style = "background-color: #f8f9fa; padding: 10px 15px; border-bottom: 1px solid #dee2e6; margin-bottom: 20px; display: flex; justify-content: space-between; align-items: center;",
    div(
      style = "display: flex; align-items: center;",
      tags$h2("PanvaR Analysis Tool", style = "margin: 0;")
    ),
    # --- START MODIFICATION: Add Restart Button ---
    div(
      # Restart Button (Green)
      actionButton("restart_app", "Restart Analysis",
                   class = "btn-success", # Green color
                   style = "margin-right: 10px;"), # Add space between buttons
      
      # Exit Button (Red)
      actionButton("exit_app", "Exit Application",
                   class = "btn-danger")
    )
    # --- END MODIFICATION ---
  ),
  tabsetPanel(
    id = "mainTabs",
    tabPanel("Input", input_dashboard_UI("module1")),
    tabPanel("Output", output_dashboard_UI("module2"))
  )
)

server <- function(input, output, session) {
  
  # Set working directory (optional, consider if needed)
  # setwd(getwd())
  
  # --- START MODIFICATION: Add Restart Logic ---
  observeEvent(input$restart_app, {
    # 1. Reset inputs in the Input module (module1)
    #    Use shinyjs::reset for inputs with default values.
    #    Note: Resetting shinyFilesButton might require different handling
    #    if you need to clear the selection display text. For now, we reset standard inputs.
    shinyjs::reset("module1-R2_threshold")
    shinyjs::reset("module1-tagSnps")
    shinyjs::reset("module1-maf")
    shinyjs::reset("module1-missing_rate")
    shinyjs::reset("module1-window_span")
    shinyjs::reset("module1-PC_min")
    shinyjs::reset("module1-PC_max")
    shinyjs::reset("module1-use_specific_pcs")
    shinyjs::reset("module1-specific_pcs")
    shinyjs::reset("module1-dynamic_correlation")
    shinyjs::reset("module1-all_impacts")
    # Resetting file inputs visually is harder; this clears the reactive value holding the path
    # You might need custom JS or observeEvent logic in the module if full visual reset is needed.
    # For now, clearing the results and switching tab is the main goal.
    
    
    # 2. Clear the shared analysis results
    shared$analysis_results <- NULL
    
    # 3. Switch back to the Input tab
    updateTabsetPanel(session, "mainTabs", selected = "Input")
    
    # 4. Optionally, show a notification
    showNotification("Inputs reset. Ready for a new analysis.", type = "message")
  })
  # --- END MODIFICATION ---
  
  
  # Exit button handler (Confirmation Modal)
  observeEvent(input$exit_app, {
    showModal(modalDialog(
      title = "Confirm Exit",
      "Are you sure you want to exit the application? Any unsaved work will be lost.",
      footer = tagList(
        actionButton("confirm_exit", "Yes, Exit", class = "btn-danger"),
        modalButton("Cancel")
      ),
      easyClose = TRUE
    ))
  })
  
  # Handle confirmed exit
  observeEvent(input$confirm_exit, {
    stopApp()
  })
  
  # Create shared reactive values
  shared <- reactiveValues(
    analysis_results = NULL
  )
  
  # Initialize modules
  # Pass the shared reactive values to the modules
  input_dashboard_Server("module1", shared) # Input module server
  output_dashboard_Server("module2", shared) # Output module server
  
  # Observe successful analysis completion and switch tabs
  observeEvent(shared$analysis_results, {
    # Check if results are not NULL (analysis finished)
    if (!is.null(shared$analysis_results)) {
      # Check if results indicate an error (e.g., if panvar_func returns specific error structure)
      # Example: if (!inherits(shared$analysis_results, "try-error")) { ... }
      # Assuming successful completion for now:
      updateTabsetPanel(session, "mainTabs", selected = "Output")
    }
  }, ignoreNULL = TRUE) # ignoreNULL = TRUE prevents triggering when results are cleared by restart
  
}

shinyApp(ui, server)