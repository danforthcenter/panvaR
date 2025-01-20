# Load required libraries
library(shiny)
library(shinyFiles)
library(plotly)

# Store the main R session environment and working directory
.GlobalEnv$.main_session <- environment()
.GlobalEnv$.main_wd <- getwd()  # Store initial working directory

# Define UI for the app
ui <- fluidPage(
  titlePanel("Panvar Analysis Dashboard"),
  
  sidebarLayout(
    sidebarPanel(
      # Required inputs section
      div(style = "border: 1px solid #ddd; padding: 10px; margin-bottom: 10px;",
          h4("Required Inputs"),
          shinyFilesButton("phenotypeDataPath", "Select Phenotype Data File", 
                          "Please select a file", multiple = FALSE),
          verbatimTextOutput("phenotypePathDisplay"),
          
          shinyFilesButton("vcfFilePath", "Select VCF File", 
                          "Please select a file", multiple = FALSE),
          verbatimTextOutput("vcfPathDisplay")
      ),
      
      # Optional inputs section
      div(style = "border: 1px solid #ddd; padding: 10px; margin-bottom: 10px;",
          h4("Optional Inputs"),
          numericInput("r2Threshold", "R2 Threshold:", value = 0.6),
          textAreaInput("tagSnps", "Tag SNPs (one per line):", ""),  # Changed to textArea
          numericInput("maf", "Minor Allele Frequency:", value = 0.05),
          numericInput("missingRate", "Missing Rate:", value = 0.10),
          numericInput("window", "Window Size:", value = 500000),
          numericInput("pcMin", "Minimum PCs:", value = 5),
          numericInput("pcMax", "Maximum PCs:", value = 5),
          checkboxInput("dynamicCorrelation", "Dynamic Correlation:", FALSE),
          checkboxInput("allImpacts", "All Impacts:", FALSE)
      ),
      
      # Command preview and execution
      div(style = "border: 1px solid #ddd; padding: 10px; margin-bottom: 10px;",
          h4("Command Preview"),
          verbatimTextOutput("commandPreview"),
          actionButton("submit", "Execute in R Session", 
                      style = "margin-top: 10px; width: 100%;"),
          actionButton("copy", "Copy Command to Clipboard",
                      style = "margin-top: 10px; width: 100%;")
      )
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Command Output",
                 verbatimTextOutput("commandOutput")
        ),
        tabPanel("Results Preview",
                 plotlyOutput("resultPlot"),  # Changed to plotlyOutput
                 tableOutput("resultTable")
        )
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  # Get the main R session environment and working directory
  main_env <- .GlobalEnv$.main_session
  main_wd <- .GlobalEnv$.main_wd
  
  # Set working directory
  setwd(main_wd)
  
  # Set up file system roots
  volumes <- c(Home = fs::path_home(), 
              "Working Directory" = main_wd,
              shinyFiles::getVolumes()())
  
  # Set up file choosers
  shinyFiles::shinyFileChoose(input, "phenotypeDataPath", roots = volumes, session = session)
  shinyFiles::shinyFileChoose(input, "vcfFilePath", roots = volumes, session = session)
  
  # Reactive values for storing file paths
  selected_files <- reactiveValues(
    phenotype = NULL,
    vcf = NULL
  )
  
  # File selection handlers with relative path conversion
  observeEvent(input$phenotypeDataPath, {
    full_path <- shinyFiles::parseFilePaths(volumes, input$phenotypeDataPath)$datapath
    if (length(full_path) > 0) {
      # Convert to relative path if possible
      selected_files$phenotype <- fs::path_rel(full_path, start = main_wd)
    }
  })
  
  observeEvent(input$vcfFilePath, {
    full_path <- shinyFiles::parseFilePaths(volumes, input$vcfFilePath)$datapath
    if (length(full_path) > 0) {
      # Convert to relative path if possible
      selected_files$vcf <- fs::path_rel(full_path, start = main_wd)
    }
  })
  
  # Display selected file paths
  output$phenotypePathDisplay <- renderText({
    if (length(selected_files$phenotype) > 0) selected_files$phenotype else "No file selected"
  })
  
  output$vcfPathDisplay <- renderText({
    if (length(selected_files$vcf) > 0) selected_files$vcf else "No file selected"
  })
  
  # Process tag SNPs input
  clean_tag_snps <- reactive({
    if (nchar(input$tagSnps) == 0) return(NULL)
    # Split by newline and clean up
    snps <- strsplit(input$tagSnps, "\n")[[1]]
    snps <- trimws(snps)  # Remove whitespace
    snps <- snps[nchar(snps) > 0]  # Remove empty lines
    snps <- gsub('"', '', snps)  # Remove any existing quotes
    snps
  })
  
  # Generate command string
  command_string <- reactive({
    # Validate required inputs
    req(selected_files$phenotype, selected_files$vcf)
    
    # Start building args list
    args <- list(
      phenotype_data_path = sprintf('"%s"', selected_files$phenotype),
      vcf_file_path = sprintf('"%s"', selected_files$vcf),
      r2_threshold = input$r2Threshold,
      maf = input$maf,
      missing_rate = input$missingRate,
      window = input$window,
      pc_min = input$pcMin,
      pc_max = input$pcMax,
      dynamic_correlation = input$dynamicCorrelation,
      all.impacts = input$allImpacts
    )
    
    # Add tag SNPs if provided
    snps <- clean_tag_snps()
    if (!is.null(snps)) {
      args$tag_snps <- sprintf('c("%s")', paste(snps, collapse = '", "'))
    }
    
    # Convert to command string
    args_str <- paste(names(args), args, sep = " = ", collapse = ",\n    ")
    sprintf("result <- panvar_func(\n    %s\n)", args_str)
  })
  
  # Display command preview
  output$commandPreview <- renderText({
    tryCatch(command_string(), error = function(e) "Waiting for required inputs...")
  })
  
  # Copy command to clipboard
  observeEvent(input$copy, {
    if (requireNamespace("clipr", quietly = TRUE)) {
      clipr::write_clip(command_string())
      showNotification("Command copied to clipboard!", type = "message")
    } else {
      showNotification("Please install 'clipr' package to use clipboard", type = "warning")
    }
  })
  
  # Execute command in main R session
  observeEvent(input$submit, {
    tryCatch({
      # Set working directory before execution
      setwd(main_wd)
      
      # Execute in main R environment
      cmd <- command_string()
      result <- eval(parse(text = cmd), envir = main_env)
      
      # Store result in main environment
      assign("result", result, envir = main_env)
      
      # Show preview in Shiny
      output$resultPlot <- renderPlotly({ 
        # Convert ggplot to plotly if necessary
        if ("ggplot" %in% class(result$plot)) {
          ggplotly(result$plot)
        } else if ("plotly" %in% class(result$plot)) {
          result$plot
        } else {
          # Create a basic plotly object if neither
          plot_ly() %>% add_annotations(
            text = "Plot format not recognized",
            showarrow = FALSE
          )
        }
      })
      
      output$resultTable <- renderTable({ result$table })
      
      # Display success message
      output$commandOutput <- renderText({
        paste("Command executed successfully in R session.",
              "The result is stored in 'result' variable in your R environment.",
              "You can access it from your R console.",
              sep = "\n")
      })
      
    }, error = function(e) {
      output$commandOutput <- renderText({
        paste("Error executing command:", e$message)
      })
    })
  })
}

# Run the application
shinyApp(ui = ui, server = server)