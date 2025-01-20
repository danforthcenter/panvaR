# Load Shiny library
library(shiny)

# Define UI for the app
ui <- fluidPage(
  # Custom CSS for styling
  tags$head(
    tags$style(HTML("
      .input-group {
        display: flex;
        align-items: center;
        margin-bottom: 5px; /* Adjust this value to reduce space between input groups */
        flex-wrap: wrap;
      }
      .input-label {
        margin-right: 10px;
        font-weight: bold;
      }
      .required-tag {
        background-color: red;
        color: white;
        border-radius: 4px;
        padding: 2px 6px;
        margin-left: 10px;
      }
      .optional-tag {
        background-color: blue;
        color: white;
        border-radius: 4px;
        padding: 2px 6px;
        margin-left: 10px;
      }
      .default-tag {
        background-color: green;
        color: white;
        border-radius: 4px;
        padding: 2px 6px;
        margin-left: 10px;
      }
      .input-box {
        border: 1px solid #ccc;
        padding: 10px;
        margin-bottom: 10px; /* Adjust this value to reduce space between input boxes */
        border-radius: 5px;
      }
      .shiny-input-container {
        width: 100%;
      }
    "))
  ),
  
  # Title of the app
  titlePanel("Panvar Analysis Dashboard"),
  
  # Sidebar layout with input and output definitions
  sidebarLayout(
    sidebarPanel(
      # Section for required inputs
      # Required inputs section
      div(class = "input-box",
          h4("Required Inputs"),
          div(class = "input-group",
              shinyFilesButton("phenotypeDataPath", "Select Phenotype Data File", "Please select a file", multiple = FALSE),
              textInput("displayPhenotypeDataPath", "Selected Phenotype Data File Path:", value = "", placeholder = "File path will appear here"),
              tags$span("Required", class = "required-tag")
          ),
          div(class = "input-group",
              shinyFilesButton("vcfFilePath", "Select VCF File", "Please select a file", multiple = FALSE),
              textInput("displayVcfFilePath", "Selected VCF File Path:", value = "", placeholder = "File path will appear here"),
              tags$span("Required", class = "required-tag")
          )
      ),
      
      # Section for optional inputs with default values
      div(class = "input-box",
          h4("Optional Inputs with Default Values"),
          div(class = "input-group",
              numericInput("r2Threshold", "R2 Threshold:", value = 0.6),
              tags$span("Default: '0.6'", class = "default-tag")
          ),
          div(class = "input-group",
              textInput("tagSnps", "Tag SNPs (comma-separated):", placeholder = "Optional"),
              tags$span("Optional", class = "optional-tag")
          ),
          div(class = "input-group",
              numericInput("maf", "Minor Allele Frequency:", value = 0.05),
              tags$span("Default: '0.05'", class = "default-tag")
          ),
          div(class = "input-group",
              numericInput("missingRate", "Missing Rate:", value = 0.10),
              tags$span("Default: '0.10'", class = "default-tag")
          ),
          div(class = "input-group",
              numericInput("window", "Window Size:", value = 500000),
              tags$span("Default: '500000'", class = "default-tag")
          ),
          div(class = "input-group",
              numericInput("pcMin", "Minimum PCs:", value = 5),
              tags$span("Default: '5'", class = "default-tag")
          ),
          div(class = "input-group",
              numericInput("pcMax", "Maximum PCs:", value = 5),
              tags$span("Default: '5'", class = "default-tag")
          ),
          div(class = "input-group",
              checkboxInput("dynamicCorrelation", "Dynamic Correlation:", value = FALSE),
              tags$span("Optional", class = "optional-tag")
          ),
          div(class = "input-group",
              checkboxInput("allImpacts", "All Impacts:", value = FALSE),
              tags$span("Default: 'FALSE'", class = "default-tag")
          ),
          div(class = "input-group",
              actionButton("kill", "Kill Process") # Add Kill button
          )
      ),
      
      # Submit button
      actionButton("submit", "Submit")
    ),
    
    # Main panel for displaying outputs
    mainPanel(
      plotOutput("resultPlot"),
      tableOutput("resultTable"),
      verbatimTextOutput("output")
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  options(shiny.maxRequestSize = 2000 * 1024^2) # Set limit to 2000 MB

  # Set up file systems
  volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), shinyFiles::getVolumes()())

  shinyFiles::shinyFileChoose(input, "phenotypeDataPath", roots = volumes, session = session, restrictions = system.file(package = "base"))
  shinyFiles::shinyFileChoose(input, "vcfFilePath", roots = volumes, session = session, restrictions = system.file(package = "base"))

  # Reactive value to track if the process should be killed
  stopProcess <- reactiveVal(FALSE)

  observeEvent(input$kill, {
    stopProcess(TRUE)
  })

  observeEvent(input$phenotypeDataPath, {
    selectedFile <- shinyFiles::parseFilePaths(volumes, input$phenotypeDataPath)
    updateTextInput(session, "displayPhenotypeDataPath", value = as.character(selectedFile$datapath))
  })

  observeEvent(input$vcfFilePath, {
    selectedFile <- shinyFiles::parseFilePaths(volumes, input$vcfFilePath)
    updateTextInput(session, "displayVcfFilePath", value = as.character(selectedFile$datapath))
  })

  observeEvent(input$submit, {
    # Reset stopProcess to FALSE when a new process starts
    stopProcess(FALSE)

    # Ensure required inputs are provided
    req(input$displayPhenotypeDataPath, input$displayVcfFilePath)

    # Get the paths to the files
    phenotypeDataPath <- input$displayPhenotypeDataPath
    vcfFilePath <- input$displayVcfFilePath

    # Check if the files exist
    if (!file.exists(phenotypeDataPath) || !file.exists(vcfFilePath)) {
      showNotification("Please provide valid file paths for both the Phenotype Data and VCF File.", type = "error")
      return(NULL)
    }

    # Convert tag SNPs input to a vector
    tagSnps <- ifelse(nchar(input$tagSnps) == 0, NULL, unlist(strsplit(input$tagSnps, ",")))

    # Run the panvar function with process tracking
    result <- tryCatch({
      # Actual Panvar Analysis Function
      panvar_func(
        phenotype_data_path = phenotypeDataPath,
        vcf_file_path = vcfFilePath,
        tag_snps = tagSnps,
        r2_threshold = input$r2Threshold,
        maf = input$maf,
        missing_rate = input$missingRate,
        window = input$window,
        pc_min = input$pcMin,
        pc_max = input$pcMax,
        dynamic_correlation = input$dynamicCorrelation,
        all.impacts = input$allImpacts
      )
    }, error = function(e) {
      showNotification(paste("Error: ", e$message), type = "error")
      return(NULL)
    })

    # Render outputs if the result is successful
    if (!is.null(result)) {
      output$resultPlot <- renderPlot({ result$plot })
      output$resultTable <- renderTable({ result$table })
      output$resultText <- renderText({ result$text })  # Assume 'result$text' holds textual output
    }
  })
}



# Run the application
shinyApp(ui = ui, server = server)
