# Rijan: Useful functions that will be used throughout the code.

# Define a function to clean the SNP tags
clean_snp_tags <- function(snp_input) {
  # Split the input string into lines
  snps <- strsplit(snp_input, "\n")[[1]]
  
  # Trim whitespace from each line
  snps <- trimws(snps)
  
  # Remove empty lines
  snps <- snps[nchar(snps) > 0]
  
  # Remove double-quote characters
  snps <- gsub('"', '', snps)
  
  return(snps)
}

parsePath <- function(input_path) {
  if (is.null(input_path)) return(NULL)
  parsed <- shinyFiles::parseFilePaths(c(Home = fs::path_home()), input_path)
  if (nrow(parsed) > 0) return(parsed$datapath[1])
  return(NULL)
}

# Check if a tbi file exists for the vcf file
check_tbi_exists <- function(vcf_file_path) {
  if (is.null(vcf_file_path)) return(FALSE)
  tbi_file_path <- paste0(vcf_file_path, ".tbi")
  return(file.exists(tbi_file_path))
}

# ---

# module1.R
# Module UI modifications for input_dashboard_UI function

input_dashboard_UI <- function(id) {
  ns <- NS(id)
  tagList(
    shinyjs::useShinyjs(),
    tabPanel(
      "PanvaR input  Dashboard"
    ),
    
    # Rijan: Split the panel into a side panel and a main panel.
    sidebarLayout(
      sidebarPanel(
        width = 3,
        div(class = "input-section",
            h4("Required Inputs"),
            # Rijan: Get the path to the genotype data.
            div(
              style = "display: flex; align-items: center; gap: 10px;",
              shinyFilesButton(
                ns("Genotype_data_path"),
                "Please select a Genotype data file",
                "Please select a file",
                multiple = FALSE
              ),
              tags$span(
                id = ns("Genotype_data_path_tooltip"),
                icon("question-circle"),
                style = "color: green;"
              )
            ),
            bsTooltip(
              id = ns("Genotype_data_path_tooltip"),
              title = "Please provide the path to a file with Genotype data.",
              placement = "right",
              trigger = "hover"
            ),
            div(
              textOutput(ns("Genotype_data_path_results"))
            ),
            div(
              style = "display: flex; align-items: center; gap: 10px;",
              shinyFilesButton(
                ns("phenotype_data"),
                "Please select a phenotype data file",
                "Please select a file",
                multiple = FALSE
              ),
              tags$span(
                id = ns("phenotype_data_tooltip"),
                icon("question-circle"),
                style = "color: green;"
              )
            ),
            bsTooltip(
              id = ns("phenotype_data_tooltip"),
              title = "Please provide the path to a Phenotype file.",
              placement = "right",
              trigger = "hover"
            ),  
            textOutput(ns("phenotype_data_results"))
        ),
        div(class = "input-section",
            h4("Default Inputs"),
            div(
              style = "display: flex; align-items: center; gap: 10px;",
              numericInput(
                ns("R2_threshold"),
                "R2 threshold:",
                value = 0.8,
                min = 0,
                max = 1
              ),
              span(
                icon("question-circle", style = "color: green;"),
                id = ns("R2_threshold_tooltip")
              )
            ),
            bsTooltip(
              id = ns("R2_threshold_tooltip"),
              title = "What is the R2 threshold that should be used? This is the value that be used to filter the LD between tag SNPs and other SNPs.",
              placement = "right",
              trigger = "hover"
            ),
            # ---
            # Rijan: The input section for the tag SNPs - Modified to make it optional
            div(
              style = "display: flex; align-items: center; gap: 10px;",
              textAreaInput(
                ns("tagSnps"),
                "Tag SNPs (optional - one per line):",
                ""
              ),
              span(
                icon("question-circle", style = "color: green;"),
                id = ns("tagSnps_tooltip")
              )
            ),
            bsTooltip(
              id = ns("tagSnps_tooltip"),
              title = "Optional: Specify tag SNPs in the format CHR:BP, one per line. If not provided, tag SNPs will be automatically inferred from GWAS results.",
              placement = "right",
              trigger = "hover"
            ),
            # ---
            # Rijan: The input section for the MAF
            div(
              style = "display: flex; align-items: center; gap: 10px;",
              numericInput(
                ns("maf"),
                "Minor Allele Frequency:",
                value = 0.05,
                min = 0,
                max = 1
              ),
              span(
                icon("question-circle", style = "color: green;"),
                id = "maf_tooltip"
              )
            ),
            bsTooltip(
              id = "maf_tooltip",
              title = "What is the Minor Allele Frequency that should be used? This is the value that be used to filter the LD between tag SNPs and other SNPs.",
              placement = "right",
              trigger = "hover"
            ),
            # ---
            # Rijan: The input section for missing rate.
            div(
              style = "display: flex; align-items: center; gap: 10px;",
              numericInput(
                ns("missing_rate"),
                "Missing rate:",
                value = 0.05,
                min = 0,
                max = 1
              ),
              span(
                icon("question-circle", style = "color: green;"),
                id = "missing_rate_tooltip"
              )
            ),
            bsTooltip(
              id = "missing_rate_tooltip",
              title = "What is the Missing rate that should be used? This is the value that be used to filter the LD between tag SNPs and other SNPs.",
              placement = "right",
              trigger = "hover"
            ),
            # ---
            # Rijan: The input section for the window span around the tag SNPs.
            div(
              style = "display: flex; align-items: center; gap: 10px;",
              numericInput(
                ns("window_span"),
                "Window span around the tag_SNP:",
                value = 500000,
                min = 0
              ),
              span(
                icon("question-circle", style = "color: green;"),
                id = "window_span_tooltip"
              )
            ),
            bsTooltip(
              id = "window_span_tooltip",
              title = "What is the window that should be used? This is the value that be used to filter the LD between tag SNPs and other SNPs.",
              placement = "right",
              trigger = "hover"
            ),
            #--- 
            # Rijan: The PC min and PC max inputs should actually be a fluidRow with two columns
            # and a numericInput for each column.
            div(
              style = "display: flex; align-items: center; gap: 10px;",
              div(
                style = "flex: 1;",
                # LD range as numeric inputs
                # Modify the PC range fluidRow
                div(
                  style = "display: flex; align-items: center; gap: 10px;",
                  div(
                    style = "flex: 1;",
                    fluidRow(
                      column(6,
                             numericInput(ns("PC_min"),
                                          "PC Min:",
                                          value = 5,
                                          min = 0,
                                          step = 1)),
                      column(6,
                             numericInput(ns("PC_max"),
                                          "PC Max:",
                                          value = 6,
                                          min = 1,
                                          step = 1))
                    )
                  ),
                  span(
                    icon("question-circle", style = "color: green;"),
                    id = "PC_range_tooltip"
                  )
                )
              )
            ),
            bsTooltip(
              id = "PC_range_tooltip",
              title = "What PC range should be used for GWAS? ",
              placement = "right",
              trigger = "hover"
            ),
            #---
            # Rijan: The input section for boolean for dynamic correlation among PCs.
            div(
              style = "display: flex; align-items: center; gap: 10px;",
              checkboxInput(
                ns("dynamic_correlation"),
                "Dynamic Correlation:",
                value = FALSE
              ),
              span(
                icon("question-circle", style = "color: green;"),
                id = "dynamic_correlation_tooltip"
              )
            ),
            # Add this after the PC range inputs but before the dynamic correlation checkbox
            div(
              style = "display: flex; align-items: center; gap: 10px;",
              div(
                style = "flex: 1;",
                checkboxInput(
                  ns("use_specific_pcs"),
                  "Use specific PCs instead of range",
                  value = FALSE
                ),
                conditionalPanel(
                  condition = sprintf("input['%s'] == true", ns("use_specific_pcs")),
                  textAreaInput(
                    ns("specific_pcs"),
                    "Specific PCs (comma-separated):",
                    value = "",
                    placeholder = "e.g., 1,3,5,7"
                  )
                )
              ),
              span(
                style = "flex-shrink: 0;",
                icon("question-circle", style = "color: green;"),
                id = "specific_pcs_tooltip"
              )
            ),
            bsTooltip(
              id = "specific_pcs_tooltip",
              title = "Choose whether to use a PC range or specify exact PCs to use. If specifying PCs, enter them as comma-separated numbers.",
              placement = "right",
              trigger = "hover"
            ),
            bsTooltip(
              id = "dynamic_correlation_tooltip",
              title = "Should the correlation among PCs be dynamic?",
              placement = "right",
              trigger = "hover"
            ),
            # ---
            # Rijan: The input section for the boolean for all.inputs
            div(
              style = "display: flex; align-items: center; gap: 10px;",
              checkboxInput(
                ns("all_impacts"),
                "All Inputs:",
                value = FALSE
              ),
              span(
                icon("question-circle", style = "color: green;"),
                id = "all_impacts_tooltip"
              )
            ),
            bsTooltip(
              id = "all_impacts_tooltip",
              title = "Should all inputs be used?",
              placement = "right",
              trigger = "hover"
            )
        )
      ),
      
      mainPanel ( # TODO: Can shiny capture the error from the underlying interpreter?
        div(
          uiOutput(ns("input_status_update"))
        ),
        actionButton(
          ns("run_analysis"),
          "Run Analysis",
          class = "btn-primary",
          style = "margin-top: 20px;"
        )
      )
    )
  )
  
}

# Module server modifications for input_dashboard_Server function

input_dashboard_Server <- function(id, shared) {
  moduleServer(
    id,
    function(input, output, session) {
      
      rootDir = c(Home = fs::path_home())
      
      # Define reactive values without req() to allow initial NULL states
      Genotype_data_path <- reactive({
        if (is.null(input$Genotype_data_path)) return(NULL)
        parsePath(input$Genotype_data_path)
      })
      
      # Add reactive value to track TBI status
      tbi_status <- reactiveVal(FALSE)
      
      # Function to generate TBI file
      generate_tbi_file <- function(vcf_path) {
        # Create a modal to show progress
        showModal(modalDialog(
          title = "Generating Index File",
          "This may take a few minutes depending on file size. Please wait...",
          footer = NULL,
          easyClose = FALSE
        ))
        
        # Run the tabix command in the background
        result <- tryCatch({
          # Create a system call to tabix
          system2("tabix", args = c("-p", "vcf", vcf_path), stdout = TRUE, stderr = TRUE)
          return(TRUE)
        }, error = function(e) {
          showNotification(
            paste("Error generating index file:", e$message),
            type = "error"
          )
          return(FALSE)
        }, finally = {
          # Ensure modal is removed in all cases
          removeModal()
        })
        
        return(result)
      }
      
      # When genotype path changes, check for TBI file
      observeEvent(Genotype_data_path(), {
        req(Genotype_data_path())
        path <- Genotype_data_path()
        
        # Check if TBI file exists
        has_tbi <- check_tbi_exists(path)
        tbi_status(has_tbi)
        
        # If TBI doesn't exist, prompt the user
        if (!has_tbi) {
          showModal(modalDialog(
            title = "Index File (.tbi) Required",
            "The selected genotype file does not have a corresponding index file (.tbi), which is required for the analysis. Would you like to generate it now?",
            footer = tagList(
              actionButton(session$ns("generate_tbi_yes"), "Yes, generate index file", class = "btn-primary"),
              actionButton(session$ns("exit_program"), "Don't generate - exit program", class = "btn-danger")
            ),
            easyClose = FALSE
          ))
        }
      })
      
      # Handle the generate TBI button click
      observeEvent(input$generate_tbi_yes, {
        req(Genotype_data_path())
        vcf_path <- Genotype_data_path()
        
        # Close the confirmation modal first
        removeModal()
        
        # Generate the TBI file
        success <- generate_tbi_file(vcf_path)
        
        # Update the status
        if (success) {
          tbi_status(TRUE)
          showNotification(
            "Index file (.tbi) has been successfully generated!",
            type = "message"
          )
        }
      })
      
      # Handle exit program button click
      observeEvent(input$exit_program, {
        # Exit the Shiny app
        stopApp()
      })
      
      phenotype_data <- reactive({
        if (is.null(input$phenotype_data)) return(NULL)
        parsePath(input$phenotype_data)
      })
      
      # For numeric inputs, return the input value directly since they have defaults
      R2_threshold <- reactive({
        input$R2_threshold
      })
      
      maf <- reactive({
        input$maf
      })
      
      missing_rate <- reactive({
        input$missing_rate
      })
      
      PC_min <- reactive({
        input$PC_min
      })
      
      PC_max <- reactive({
        input$PC_max
      })
      
      # For boolean inputs, they'll always have a value due to defaults
      dynamic_correlation <- reactive({
        input$dynamic_correlation
      })
      
      all_impacts <- reactive({
        input$all_impacts
      })
      
      window_span <- reactive({
        input$window_span
      })
      
      # File chooser setup remains the same
      shinyFileChoose(
        input,
        "Genotype_data_path",
        roots = rootDir,
        session = session,
        filetypes = c("gz","psam")
      )
      
      output$Genotype_data_path_results <- renderText({
        path <- Genotype_data_path()
        if (!is.null(path)) {
          paste("You have entered", path)
        }
      })
      
      shinyFileChoose(
        input,
        "phenotype_data",
        roots = rootDir,
        session = session,
        filetypes = c("csv","txt","tsv")
      )
      
      output$phenotype_data_results <- renderText({
        path <- phenotype_data()
        if (!is.null(path)) {
          paste("You have entered", path)
        }
      })
      
      clean_tag_snps <- reactive({
        if (is.null(input$tagSnps) || input$tagSnps == "") return(NULL)
        clean_snp_tags(input$tagSnps)
      })
      
      # In your server function
      observe({
        if (input$use_specific_pcs) {
          shinyjs::disable("PC_min")
          shinyjs::disable("PC_max")
        } else {
          shinyjs::enable("PC_min")
          shinyjs::enable("PC_max")
        }
      })
      
      # Add these new reactive values
      use_specific_pcs <- reactive({
        input$use_specific_pcs
      })
      
      specific_pcs <- reactive({
        if (!use_specific_pcs()) return(NULL)
        if (is.null(input$specific_pcs) || input$specific_pcs == "") return(NULL)
        
        # Parse comma-separated string into numeric vector
        pcs <- tryCatch({
          as.numeric(unlist(strsplit(gsub("\\s+", "", input$specific_pcs), ",")))
        }, error = function(e) NULL)
        
        # Validate that all values are positive integers
        if (!is.null(pcs) && all(!is.na(pcs)) && all(pcs > 0) && all(pcs == floor(pcs))) {
          return(sort(unique(pcs)))
        }
        return(NULL)
      })
      
      output$input_status_update <- renderUI({
        # Create status lists for both current values and missing inputs
        missing_inputs <- character(0)
        current_values <- list()
        
        # Check and store file paths
        if (!is.null(Genotype_data_path())) {
          current_values$`Genotype Data File` <- Genotype_data_path()
          
          # Check TBI status
          if (!tbi_status()) {
            missing_inputs <- c(missing_inputs, "Index file (.tbi) for the genotype data")
          }
        } else {
          missing_inputs <- c(missing_inputs, "Genotype Data File")
        }
        
        if (!is.null(phenotype_data())) {
          current_values$`Phenotype Data File` <- phenotype_data()
        } else {
          missing_inputs <- c(missing_inputs, "Phenotype Data File")
        }
        
        # Store numeric inputs with their current values
        current_values$`R² Threshold` <- sprintf("%.3f", R2_threshold())
        if (R2_threshold() < 0 || R2_threshold() > 1) {
          missing_inputs <- c(missing_inputs, "R² threshold (must be between 0 and 1)")
        }
        
        current_values$`Minor Allele Frequency` <- sprintf("%.3f", maf())
        if (maf() < 0 || maf() > 1) {
          missing_inputs <- c(missing_inputs, "MAF (must be between 0 and 1)")
        }
        
        current_values$`Missing Rate` <- sprintf("%.3f", missing_rate())
        current_values$`Window Span` <- format(window_span(), big.mark = ",")
        
        # Only add PC-related information based on the selected method
        if (use_specific_pcs()) {
          current_values$`PC Selection Method` <- "Specific PCs"
          if (!is.null(specific_pcs())) {
            current_values$`Selected PCs` <- paste(specific_pcs(), collapse = ", ")
          } else {
            missing_inputs <- c(missing_inputs, "Valid specific PCs (must be positive integers)")
          }
        } else {
          current_values$`PC Selection Method` <- "PC Range"
          current_values$`PC Range` <- sprintf("%d - %d", PC_min(), PC_max())
          if (PC_min() > PC_max()) { # > is wrong but = is fine.
            missing_inputs <- c(missing_inputs, "PC range (min must be less than or equal to max)")
          }
        }
        
        # Store boolean values
        current_values$`Dynamic Correlation` <- if(dynamic_correlation()) "Yes" else "No"
        current_values$`All Impacts` <- if(all_impacts()) "Yes" else "No"
        
        # Store tag SNPs - CHANGED: now tag SNPs are optional
        if (!is.null(clean_tag_snps())) {
          current_values$`Tag SNPs` <- paste(clean_tag_snps(), collapse = ", ")
        } else {
          current_values$`Tag SNPs` <- "Auto-infer from GWAS results"
        }
        
        # Create the UI
        tagList(
          div(
            class = "panel panel-default",
            style = "margin-top: 20px;",
            div(
              class = "panel-heading",
              h4(class = "panel-title", "Current Input Values")
            ),
            div(
              class = "panel-body",
              tags$table(
                class = "table table-striped table-bordered",
                tags$thead(
                  tags$tr(
                    tags$th("Parameter"),
                    tags$th("Value")
                  )
                ),
                tags$tbody(
                  lapply(names(current_values), function(param) {
                    tags$tr(
                      tags$td(strong(param)),
                      tags$td(current_values[[param]])
                    )
                  })
                )
              )
            )
          ),
          
          div(
            class = if(length(missing_inputs) > 0) "alert alert-warning" else "alert alert-success",
            style = "margin-top: 20px;",
            if(length(missing_inputs) > 0) {
              tagList(
                h4("Please provide or correct the following:"),
                tags$ul(
                  lapply(missing_inputs, function(x) {
                    tags$li(x)
                  })
                )
              )
            } else {
              tagList(
                h4("All inputs are valid!"),
                p("You can proceed with the analysis.")
              )
            }
          )
        )
      })
      
      # Create reactive values to store analysis results
      analysis_results <- reactiveVal(NULL)
      
      # Observe the run analysis button
      observeEvent(input$run_analysis, {
        # Validate inputs before running - CHANGED: Added TBI validation
        all_inputs_valid <- reactive({
          # Return TRUE only if all conditions are met
          !is.null(Genotype_data_path()) &&
            tbi_status() && # Check TBI status
            !is.null(phenotype_data()) &&
            R2_threshold() >= 0 && R2_threshold() <= 1 &&
            maf() >= 0 && maf() <= 1 &&
            missing_rate() >= 0 && missing_rate() <= 1 &&
            window_span() > 0 &&
            (use_specific_pcs() || PC_min() <= PC_max())
        })
        
        
        # If validation passes, run the analysis
        if (all_inputs_valid()) {
          withProgress(message = 'Running analysis...', value = 0, {
            
            # Prepare inputs for panvar_func
            pc_params <- if(use_specific_pcs()) {
              list(pc_min = min(specific_pcs()), 
                   pc_max = max(specific_pcs()))
            } else {
              list(pc_min = PC_min(), 
                   pc_max = PC_max())
            }
            
            # Run panvar_func with user inputs - pass NULL for tag_snps if not provided
            results <- try({
              panvar_func(
                phenotype_data = phenotype_data(),
                vcf_file_path = Genotype_data_path(),
                tag_snps = clean_tag_snps(),  # This will be NULL if no tag SNPs are provided
                r2_threshold = R2_threshold(),
                maf = maf(),
                missing_rate = missing_rate(),
                window = window_span(),
                pc_min = pc_params$pc_min,
                pc_max = pc_params$pc_max,
                dynamic_correlation = dynamic_correlation(),
                all.impacts = all_impacts()
              )
            })
            
            if (inherits(results, "try-error")) {
              showNotification(
                "Error in analysis. Please check your inputs and try again.",
                type = "error"
              )
            } else {
              # Store results in shared reactive values for access by output dashboard
              shared$analysis_results <- results
              
              showNotification(
                "Analysis completed successfully!",
                type = "message"  # Changed from "success" to "message"
              )
            }
          })
        } else {
          # Show error if TBI file is missing
          if (!tbi_status()) {
            # Check if we have a genotype file but no TBI
            if (!is.null(Genotype_data_path())) {
              # Prompt again to generate TBI
              showModal(modalDialog(
                title = "Index File (.tbi) Required",
                "You need to generate or provide an index file (.tbi) before running the analysis. Would you like to generate it now?",
                footer = tagList(
                  actionButton(session$ns("generate_tbi_yes"), "Yes, generate index file", class = "btn-primary"),
                  actionButton(session$ns("exit_program"), "Don't generate - exit program", class = "btn-danger")
                ),
                easyClose = FALSE
              ))
            } else {
              showNotification(
                "Please select a genotype file first.",
                type = "error"
              )
            }
          } else {
            showNotification(
              "Please check all inputs before running the analysis.",
              type = "error"
            )
          }
        }
      })
      
      # Return reactive values for use in output dashboard
      return(list(
        analysis_results = analysis_results
      ))
      
    }
  )
}