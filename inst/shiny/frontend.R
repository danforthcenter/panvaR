# Load required libraries
library(shiny)
library(shinyFiles)
library(plotly)

# Store the main R session environment and working directory
.GlobalEnv$.main_session <- environment()
.GlobalEnv$.main_wd <- getwd()

# Custom CSS for styling
css <- tags$style(HTML("
  .required-label::after {
    content: ' *';
    color: red;
  }
  
  .required-field {
    border-color: #ff6b6b !important;
  }
  
  .default-hint {
    font-size: 0.8em;
    color: #28a745;
    margin-left: 5px;
  }
  
  .optional-hint {
    font-size: 0.8em;
    color: #007bff;
    margin-left: 5px;
  }
  
  .input-section {
    border: 1px solid #ddd;
    border-radius: 5px;
    padding: 15px;
    margin-bottom: 15px;
    background-color: #f8f9fa;
  }
  
  .kill-button {
    background-color: #dc3545 !important;
    color: white !important;
    margin-top: 10px;
  }
  
  .action-button {
    margin-top: 10px;
    width: 100%;
  }
"))

# Define UI for the app
ui <- fluidPage(
  css,
  titlePanel("Panvar Analysis Dashboard"),
  
  sidebarLayout(
    sidebarPanel(
      # Required inputs section
      div(class = "input-section",
          h4("Required Inputs"),
          
          div(
            tags$label("Phenotype Data File", class = "required-label"),
            shinyFilesButton("phenotypeDataPath", "Select Phenotype Data File", 
                           "Please select a file", multiple = FALSE),
            verbatimTextOutput("phenotypePathDisplay")
          ),
          
          div(style = "margin-top: 15px;",
              tags$label("VCF File", class = "required-label"),
              shinyFilesButton("vcfFilePath", "Select VCF File", 
                             "Please select a file", multiple = FALSE),
              verbatimTextOutput("vcfPathDisplay")
          )
      ),
      
      # Optional inputs section
      div(class = "input-section",
          h4("Optional Inputs"),
          
          div(
            numericInput("r2Threshold", "R2 Threshold:", value = 0.6, step = 0.1),
            tags$span("Default: '0.6'", class = "default-hint")
          ),
          
          div(
            textAreaInput("tagSnps", "Tag SNPs (one per line):", ""),
            tags$span("Optional", class = "optional-hint")
          ),
          
          div(
            numericInput("maf", "Minor Allele Frequency:", value = 0.05, step = 0.01),
            tags$span("Default: '0.05'", class = "default-hint")
          ),
          
          div(
            numericInput("missingRate", "Missing Rate:", value = 0.10, step = 0.01),
            tags$span("Default: '0.10'", class = "default-hint")
          ),
          
          div(
            numericInput("window", "Window Size:", value = 500000, step = 1000),
            tags$span("Default: '500000'", class = "default-hint")
          ),
          
          div(
            numericInput("pcMin", "Minimum PCs:", value = 5),
            tags$span("Default: '5'", class = "default-hint")
          ),
          
          div(
            numericInput("pcMax", "Maximum PCs:", value = 5),
            tags$span("Default: '5'", class = "default-hint")
          ),
          
          div(
            checkboxInput("dynamicCorrelation", "Dynamic Correlation:", FALSE),
            tags$span("Optional", class = "optional-hint")
          ),
          
          div(
            checkboxInput("allImpacts", "All Impacts:", FALSE),
            tags$span("Default: 'FALSE'", class = "default-hint")
          )
      ),
      
      # Command preview and execution
      div(class = "input-section",
          h4("Command Preview"),
          verbatimTextOutput("commandPreview"),
          actionButton("submit", "Execute in R Session", 
                      class = "action-button btn-primary"),
          actionButton("copy", "Copy Command to Clipboard",
                      class = "action-button btn-info"),
          actionButton("kill", "Kill Process", 
                      class = "action-button kill-button")
      )
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Command Output",
                 verbatimTextOutput("commandOutput")
        ),
        tabPanel("Results Preview",
                 plotlyOutput("resultPlot"),
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
  
  # Reactive values for storing file paths, process ID, and results
  values <- reactiveValues(
    phenotype = NULL,
    vcf = NULL,
    process_id = NULL,
    result = NULL,
    error = NULL
  )
  
  # File selection handlers with relative path conversion
  observeEvent(input$phenotypeDataPath, {
    full_path <- shinyFiles::parseFilePaths(volumes, input$phenotypeDataPath)$datapath
    if (length(full_path) > 0) {
      values$phenotype <- fs::path_rel(full_path, start = main_wd)
    }
  })
  
  observeEvent(input$vcfFilePath, {
    full_path <- shinyFiles::parseFilePaths(volumes, input$vcfFilePath)$datapath
    if (length(full_path) > 0) {
      values$vcf <- fs::path_rel(full_path, start = main_wd)
    }
  })
  
  # Display selected file paths
  output$phenotypePathDisplay <- renderText({
    if (length(values$phenotype) > 0) values$phenotype else "No file selected"
  })
  
  output$vcfPathDisplay <- renderText({
    if (length(values$vcf) > 0) values$vcf else "No file selected"
  })
  
  # Process tag SNPs input
  clean_tag_snps <- reactive({
    if (nchar(input$tagSnps) == 0) return(NULL)
    snps <- strsplit(input$tagSnps, "\n")[[1]]
    snps <- trimws(snps)
    snps <- snps[nchar(snps) > 0]
    snps <- gsub('"', '', snps)
    snps
  })
  
  # Generate command string
  command_string <- reactive({
    req(values$phenotype, values$vcf)
    
    args <- list(
      phenotype_data_path = sprintf('"%s"', values$phenotype),
      vcf_file_path = sprintf('"%s"', values$vcf),
      r2_threshold = input$r2Threshold,
      maf = input$maf,
      missing_rate = input$missingRate,
      window = input$window,
      pc_min = input$pcMin,
      pc_max = input$pcMax,
      dynamic_correlation = input$dynamicCorrelation,
      all.impacts = input$allImpacts
    )
    
    snps <- clean_tag_snps()
    if (!is.null(snps)) {
      args$tag_snps <- sprintf('c("%s")', paste(snps, collapse = '", "'))
    }
    
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
    values$error <- NULL
    values$result <- NULL
    
    tryCatch({
      setwd(main_wd)
      cmd <- command_string()
      
      # Execute in main R environment
      result <- eval(parse(text = cmd), envir = main_env)
      
      # Store result in both environments
      values$result <- result
      assign("result", result, envir = main_env)
      
      output$commandOutput <- renderText({
        paste("Command executed successfully in R session.",
              "The result is stored in 'result' variable in your R environment.",
              "You can access it from your R console.",
              sep = "\n")
      })
      
    }, error = function(e) {
      values$error <- e$message
      output$commandOutput <- renderText({
        paste("Error executing command:", e$message)
      })
    })
  })
  
  # Render plot only when result is available
  output$resultPlot <- renderPlotly({
    req(values$result)
    if (is.null(values$error)) {
      if ("ggplot" %in% class(values$result$plot)) {
        ggplotly(values$result$plot)
      } else if ("plotly" %in% class(values$result$plot)) {
        values$result$plot
      } else {
        plot_ly() %>% add_annotations(
          text = "Plot format not recognized",
          showarrow = FALSE
        )
      }
    }
  })
  
  # Render table only when result is available
  output$resultTable <- renderTable({
    req(values$result)
    if (is.null(values$error)) {
      values$result$table
    }
  })
  
  # Kill process handler
  observeEvent(input$kill, {
    if (!is.null(values$process_id)) {
      if (requireNamespace("tools", quietly = TRUE)) {
        tools::pskill(values$process_id)
        showNotification("Process killed", type = "message")
        values$process_id <- NULL
      } else {
        showNotification("Could not kill process - 'tools' package not available", 
                        type = "error")
      }
    } else {
      showNotification("No active process to kill", type = "warning")
    }
  })
}

# Run the application
shinyApp(ui = ui, server = server)
