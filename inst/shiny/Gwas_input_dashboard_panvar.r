# Gwas_input_dashboard_panvar.r
# Module for a direct analysis run starting from a pre-computed GWAS table.

# Helper function to clean SNP tags (can be shared or defined here)
clean_snp_tags_gwas <- function(snp_input) {
  if (is.null(snp_input) || snp_input == "") return(NULL)
  snps <- strsplit(snp_input, "\n")[[1]]
  snps <- trimws(snps)
  snps <- snps[nchar(snps) > 0]
  snps <- gsub('"', '', snps)
  return(snps)
}

# Helper to parse file paths
parsePath_gwas <- function(input_path) {
  if (is.null(input_path)) return(NULL)
  # Ensure the input is a list-like structure as expected by parseFilePaths
  if (!is.list(input_path)) {
    # This handles cases where the input might not be in the expected format
    return(NULL)
  }
  parsed <- shinyFiles::parseFilePaths(c(Home = fs::path_home()), input_path)
  if (nrow(parsed) > 0) return(parsed$datapath[1])
  return(NULL)
}

# --- Module UI ---
Gwas_input_dashboard_UI <- function(id) {
  ns <- NS(id)
  tagList(
    shinyjs::useShinyjs(),
    sidebarLayout(
      sidebarPanel(
        width = 3,
        div(class = "input-section",
            h4("Required Inputs"),
            # VCF File Input
            div(
              style = "display: flex; align-items: center; gap: 10px;",
              shinyFilesButton(
                ns("vcf_file_path"),
                "Select VCF Genotype File (.vcf.gz)",
                "Please select a file",
                multiple = FALSE
              ),
              tags$span(id = ns("vcf_file_tooltip"), icon("question-circle"), style = "color: green;")
            ),
            bsTooltip(ns("vcf_file_tooltip"), "Select the SnpEff-annotated and bgzip-compressed VCF file.", "right", "hover"),
            textOutput(ns("vcf_file_path_results")),
            
            # GWAS Table Input
            div(
              style = "display: flex; align-items: center; gap: 10px;",
              shinyFilesButton(
                ns("gwas_table_path"),
                "Select GWAS Results File",
                "Please select a file",
                multiple = FALSE
              ),
              tags$span(id = ns("gwas_table_tooltip"), icon("question-circle"), style = "color: green;")
            ),
            bsTooltip(ns("gwas_table_tooltip"), "Select a file with GWAS results. Must contain CHROM, BP, and Pvalues columns.", "right", "hover"),
            textOutput(ns("gwas_table_path_results"))
        ),
        
        div(class = "input-section",
            h4("Analysis Parameters"),
            # Tag SNPs Input
            div(
              style = "display: flex; align-items: center; gap: 10px;",
              textAreaInput(ns("tag_snps"), "Tag SNPs (one per line):", ""),
              tags$span(id = ns("tag_snps_tooltip"), icon("question-circle"), style = "color: green;")
            ),
            bsTooltip(ns("tag_snps_tooltip"), "Specify tag SNPs in CHR:BP format. If empty, the top hit from the GWAS table will be used.", "right", "hover"),
            
            # R2 Threshold
            div(
              style = "display: flex; align-items: center; gap: 10px;",
              numericInput(ns("r2_threshold"), "R² threshold:", value = 0.6, min = 0.01, max = 0.99, step = 0.05),
              tags$span(id = ns("r2_tooltip"), icon("question-circle"), style = "color: green;")
            ),
            bsTooltip(ns("r2_tooltip"), "LD threshold (R²) to find correlated SNPs.", "right", "hover"),
            
            # MAF
            div(
              style = "display: flex; align-items: center; gap: 10px;",
              numericInput(ns("maf"), "Minor Allele Frequency:", value = 0.05, min = 0, max = 1),
              tags$span(id = ns("maf_tooltip"), icon("question-circle"), style = "color: green;")
            ),
            bsTooltip(ns("maf_tooltip"), "Minimum Minor Allele Frequency for filtering.", "right", "hover"),
            
            # Missing Rate
            div(
              style = "display: flex; align-items: center; gap: 10px;",
              numericInput(ns("missing_rate"), "Missing rate:", value = 0.1, min = 0, max = 1),
              tags$span(id = ns("missing_rate_tooltip"), icon("question-circle"), style = "color: green;")
            ),
            bsTooltip(ns("missing_rate_tooltip"), "Maximum allowed missing rate for variants.", "right", "hover"),
            
            # Window Span
            div(
              style = "display: flex; align-items: center; gap: 10px;",
              numericInput(ns("window_span"), "Window around tag SNP (bp):", value = 500000, min = 0),
              tags$span(id = ns("window_span_tooltip"), icon("question-circle"), style = "color: green;")
            ),
            bsTooltip(ns("window_span_tooltip"), "The genomic window (in base pairs) to analyze around each tag SNP.", "right", "hover"),
            
            # All Impacts Checkbox
            div(
              style = "display: flex; align-items: center; gap: 10px;",
              checkboxInput(ns("all_impacts"), "Include all SNP impacts", value = FALSE),
              tags$span(id = ns("all_impacts_tooltip"), icon("question-circle"), style = "color: green;")
            ),
            bsTooltip(ns("all_impacts_tooltip"), "If checked, includes LOW and MODIFIER impacts. If unchecked, only HIGH and MODERATE impacts are considered.", "right", "hover")
        )
      ), # End sidebarPanel
      mainPanel(
        uiOutput(ns("input_status_update")),
        actionButton(ns("run_analysis"), "Run Analysis", class = "btn-primary", style = "margin-top: 20px;")
      ) # End mainPanel
    ) # End sidebarLayout
  ) # End tagList
}

# --- Module Server ---
Gwas_input_dashboard_Server <- function(id, shared) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    rootDir <- c(Home = fs::path_home())
    
    # --- Reactives for File Inputs ---
    vcf_file_path <- reactive({ parsePath_gwas(input$vcf_file_path) })
    gwas_table_path <- reactive({ parsePath_gwas(input$gwas_table_path) })
    
    # --- TBI status reactive for the VCF file ---
    tbi_status <- reactiveVal(FALSE)
    
    # --- File Chooser Setup ---
    shinyFileChoose(input, "vcf_file_path", roots = rootDir, session = session, filetypes = c("gz"))
    shinyFileChoose(input, "gwas_table_path", roots = rootDir, session = session, filetypes = c("csv", "txt", "tsv"))
    
    # --- Display Selected File Paths ---
    output$vcf_file_path_results <- renderText({
      path <- vcf_file_path()
      if (!is.null(path)) paste("VCF selected:", basename(path))
    })
    output$gwas_table_path_results <- renderText({
      path <- gwas_table_path()
      if (!is.null(path)) paste("GWAS file selected:", basename(path))
    })
    
    # --- Check for TBI file when VCF path changes ---
    observeEvent(vcf_file_path(), {
      path <- vcf_file_path()
      if (is.null(path) || !file.exists(path)) {
        tbi_status(FALSE)
        return()
      }
      if (!endsWith(path, ".gz")) {
        showNotification("Warning: The selected VCF file is not gzipped (.gz). Indexing is required and only works on bgzip-compressed files.", type = "warning", duration=8)
        tbi_status(FALSE)
        return()
      }
      
      has_tbi <- file.exists(paste0(path, ".tbi"))
      tbi_status(has_tbi)
      if (!has_tbi) {
        showModal(modalDialog(
          title = "VCF Index File (.tbi) Missing",
          "The selected VCF file is missing its index (.tbi), which is required. Would you like to create it now? This may take a moment.",
          footer = tagList(
            actionButton(ns("generate_tbi_yes"), "Yes, Create Index", class = "btn-primary"),
            modalButton("Cancel")
          )
        ))
      }
    })
    
    # --- Handle TBI Generation ---
    observeEvent(input$generate_tbi_yes, {
      removeModal()
      vcf_path <- vcf_file_path()
      req(vcf_path)
      
      withProgress(message = 'Generating VCF index...', value = 0.5, {
        tryCatch({
          generate_tbi_file(vcf_path) # Assumes generate_tbi_file is available from the package
          tbi_status(TRUE)
          showNotification("VCF index file (.tbi) created successfully!", type = "message")
        }, error = function(e) {
          showNotification(paste("Error creating index:", e$message, ". Please ensure bcftools/tabix is installed and in your PATH."), type = "error", duration = 10)
          tbi_status(FALSE)
        })
      })
    })
    
    # --- Reactive for Input Validation ---
    all_inputs_valid <- reactive({
      !is.null(vcf_file_path()) && !is.null(gwas_table_path()) &&
        file.exists(vcf_file_path()) && file.exists(gwas_table_path()) &&
        tbi_status() == TRUE &&
        is.numeric(input$r2_threshold) && input$r2_threshold >= 0 && input$r2_threshold <= 1 &&
        is.numeric(input$maf) && input$maf >= 0 && input$maf <= 1 &&
        is.numeric(input$missing_rate) && input$missing_rate >= 0 && input$missing_rate <= 1 &&
        is.numeric(input$window_span) && input$window_span >= 0
    })
    
    # --- UI for Input Status ---
    output$input_status_update <- renderUI({
      missing_inputs <- character(0)
      current_values <- list()
      
      if (!is.null(vcf_file_path()) && file.exists(vcf_file_path())) {
        current_values$`VCF File` <- basename(vcf_file_path())
        if(!tbi_status()) {
          missing_inputs <- c(missing_inputs, "VCF file index (.tbi) is missing")
        }
      } else {
        missing_inputs <- c(missing_inputs, "A valid VCF file (.vcf.gz)")
      }
      
      if (!is.null(gwas_table_path()) && file.exists(gwas_table_path())) {
        current_values$`GWAS Table` <- basename(gwas_table_path())
      } else {
        missing_inputs <- c(missing_inputs, "A valid GWAS results file")
      }
      
      tag_snps <- clean_snp_tags_gwas(input$tag_snps)
      current_values$`Tag SNPs` <- if (!is.null(tag_snps)) paste(tag_snps, collapse = ", ") else "Auto-infer from GWAS results"
      
      current_values$`R² Threshold` <- input$r2_threshold
      current_values$`Minor Allele Frequency` <- input$maf
      current_values$`Missing Rate` <- input$missing_rate
      current_values$`Window Span (bp)` <- format(input$window_span, big.mark = ",")
      current_values$`Include All Impacts` <- if(input$all_impacts) "Yes" else "No"
      
      tagList(
        div(
          class = "panel panel-default", style = "margin-top: 20px;",
          div(class = "panel-heading", h4(class = "panel-title", "Current Settings")),
          div(
            class = "panel-body",
            tags$table(
              class = "table table-striped table-bordered",
              tags$thead(tags$tr(tags$th("Parameter"), tags$th("Value"))),
              tags$tbody(
                lapply(names(current_values), function(param) {
                  tags$tr(tags$td(strong(param)), tags$td(current_values[[param]]))
                })
              )
            )
          )
        ),
        div(
          class = if(length(missing_inputs) > 0) "alert alert-warning" else "alert alert-success",
          style = "margin-top: 20px;",
          if(length(missing_inputs) > 0) {
            tagList(h4("Please provide or correct the following:"), tags$ul(lapply(missing_inputs, tags$li)))
          } else {
            tagList(h4("All inputs are valid!"), p("You can proceed with the analysis."))
          }
        )
      )
    })
    
    # --- Run Analysis Observer ---
    observeEvent(input$run_analysis, {
      on.exit(shinyjs::enable("run_analysis"))
      
      if (all_inputs_valid()) {
        shinyjs::disable("run_analysis")
        
        withProgress(message = 'Running analysis from GWAS table...', value = 0, {
          results <- tryCatch({
            incProgress(0.1, detail = "Initializing analysis...")
            
            # Call the modified panvar_func
            panvar_func(
                vcf_file_path = vcf_file_path(),
                gwas_table_path = gwas_table_path(),
                phenotype_data = NULL, # Explicitly pass NULL for phenotype
                tag_snps = clean_snp_tags_gwas(input$tag_snps),
                r2_threshold = input$r2_threshold,
                maf = input$maf,
                missing_rate = input$missing_rate,
                window = input$window_span,
                all.impacts = input$all_impacts
            )
            
          }, error = function(e) {
            showNotification(paste("Error during analysis:", e$message), type = "error", duration = 10)
            NULL
          })
          
          if (!is.null(results)) {
            shared$analysis_results <- results
            showNotification("Analysis completed successfully!", type = "message")
          } else {
            shared$analysis_results <- NULL
          }
        })
      } else {
        showNotification("Please provide all required and valid inputs before running the analysis.", type = "error")
      }
    })
  })
}