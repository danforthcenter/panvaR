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
      div(class = "input-box",
          h4("Required Inputs"),
          div(class = "input-group",
              fileInput("phenotypeData", "Upload Phenotype Data:", accept = c(".csv", ".tsv", ".txt")),
              tags$span("Required", class = "required-tag")
          ),
          div(class = "input-group",
              fileInput("vcfFile", "Upload VCF File:", accept = c(".vcf", ".vcf.gz")),
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
server <- function(input, output) {

  # So, turns out shiny limits file upload size to 
  # please see:
  # https://stackoverflow.com/questions/18037737/how-to-change-maximum-upload-size-exceeded-restriction-in-shiny-and-save-user and
  # https://groups.google.com/g/shiny-discuss/c/rU3vwGMZexQ/m/zeKhiYXrtEQJ?pli=1
  # So, for now I will manually set this to 2000 and fix it later.
  # TODO: options(shiny.maxRequestSize = 2000 * 1024^2) is ducktape and needs a better solution.

  options(shiny.maxRequestSize = 2000 * 1024^2) # Set limit to 2000 MB

  observeEvent(input$submit, {
    # Ensure required inputs are provided
    req(input$phenotypeData, input$vcfFile)
    
    # Get the paths to the uploaded files
    phenotypeDataPath <- input$phenotypeData$datapath
    vcfFilePath <- input$vcfFile$datapath
    
    # Check if the files have been uploaded correctly
    if (is.null(phenotypeDataPath) || is.null(vcfFilePath)) {
      showNotification("Please upload both the Phenotype Data and VCF File.", type = "error")
      return(NULL)
    }
    
    # Convert tag SNPs input to a vector
    tagSnps <- ifelse(nchar(input$tagSnps) == 0, NULL, unlist(strsplit(input$tagSnps, ",")))
    
    # Run the panvar function
    result <- tryCatch({
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
    
    # Render plot and table if the result is successful
    if (!is.null(result)) {
      output$resultPlot <- renderPlot({ result$plot })
      output$resultTable <- renderTable({ result$table })
    }
  })
}

# Run the application
shinyApp(ui = ui, server = server)
