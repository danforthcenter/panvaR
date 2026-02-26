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
  # bsplus::use_bs_tooltip()
  tagList(
    shinyjs::useShinyjs(),
    sidebarLayout(
      sidebarPanel(
        width = 3,
        div(class = "input-section",
            h4("Input Files"),
            # bed File Input
            div(
              style = "display: flex; align-items: center; gap: 10px;",
              shinyFilesButton(
                ns("bed_file_path"),
                "Select Genotype File (.bed)",
                "Select a genotype file in bed format",
                multiple = FALSE
              ),
              tags$span(id = ns("bed_file_tooltip"), icon("question-circle"), style = "color: green;")
            ),
            shinyBS::bsTooltip(ns("bed_file_tooltip"), "Select the bedfile produced from panvaR::make_panvar_inputs()", "right", "hover"),
            textOutput(ns("bed_file_path_results")),
            
            # GWAS Table Input
            div(
              style = "display: flex; align-items: center; gap: 10px;",
              shinyFilesButton(
                ns("gwas_table_path"),
                "Select GWAS Results File",
                "Select a table of GWAS results",
                multiple = FALSE
              ),
              tags$span(id = ns("gwas_table_tooltip"), icon("question-circle"), style = "color: green;")
            ),
            shinyBS::bsTooltip(ns("gwas_table_tooltip"), "Select a file with GWAS results. Must contain CHR, POS, PVAL columns.", "right", "hover"),
            textOutput(ns("gwas_table_path_results")),
            
            # annotation table
            div(
              style = "display: flex; align-items: center; gap: 10px;",
              shinyFilesButton(
                ns("annotation_table_path"),
                "Select Annotation Table File",
                "Select a table of gene annotations",
                multiple = FALSE
              ),
              tags$span(id = ns("anno_table_tooltip"), icon("question-circle"), style = "color: green;")
            ),
            bsTooltip(ns("anno_table_tooltip"), "Select an annotation table. Must contain columns geneID, CHR, start, end, annotation.", "right", "hover"),
            textOutput(ns("anno_table_path_results"))
        ),
        
        div(class = "input-section",
            h4("Plink path"),
            
            # Plink path
            div(
              style = "display: flex; align-items: center; gap: 10px;",
              shinyFilesButton(
                ns("plink_path"),
                "Select plink2 executable",
                "Select plink2 executable",
                multiple = FALSE
              ),
              tags$span(id = ns("plink_path_tooltip"), icon("question-circle"), style = "color: green;")
            ),
            bsTooltip(ns("plink_path_tooltip"), "Select the bedfile produced from panvaR::make_panvar_inputs()", "right", "hover"),
            textOutput(ns("plink_path_results"))
        ),
        
        div(class = "input-section",
            h4("Required inputs"),
            
            # Tag SNPs Input
            div(
              style = "display: flex; align-items: center; gap: 10px;",
              textInput(ns("tag_snp"), "Tag SNP:", value = "Chr_05-6857045"),
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
            
            # Pvals in log 
            div(
              style = "display: flex; align-items: center; gap: 10px;",
              checkboxInput(ns("pvals_in_log"), "P-values in Log?", value = T),
              tags$span(id = ns("pvals_in_log_tooltip"), icon("question-circle"), style = "color: green;")
            ),
            bsTooltip(ns("pvals_in_log_tooltip"), "Are values in PVAL column already in log10?", "right", "hover"),
            
            # Window 
            div(
              style = "display: flex; align-items: center; gap: 10px;",
              numericInput(ns("window_span"), "Window around tag SNP in kilobases:", value = 100, min = 0),
              tags$span(id = ns("window_span_tooltip"), icon("question-circle"), style = "color: green;")
            ),
            bsTooltip(ns("window_span_tooltip"), "The physical distance (in kilobases) on either side of tag SNP.", "right", "hover")
        ),
        
        div(class = "input-section",
            h4("Required inputs"),
            
            # compute scores
            div(
              style = "display: flex; align-items: center; gap: 10px;",
              checkboxInput(ns("compute_scores"), "Compute scores?", value = F),
              tags$span(id = ns("compute_scores_tooltip"), icon("question-circle"), style = "color: green;")
            ),
            bsTooltip(ns("compute_scores_tooltip"), "Compute snp scores, an aggregate of some variables?", "right", "hover"),
            
            # score.vars
            
            # score.dirs
            
            # score.weights
            
            # snp.to.gene.vars
            
            # snp.to.gene.buffer
            div(
              style = "display: flex; align-items: center; gap: 10px;",
              numericInput(ns("snp.to.gene.buffer"), "snp to gene buffer", value = 5, min = 0),
              tags$span(id = ns("snp.to.gene.buffer_tooltip"), icon("question-circle"), style = "color: green;")
            ),
            bsTooltip(ns("snp.to.gene.buffer_tooltip"), "Buffer around genes in KB to associate a snp with a gene", "right", "hover")
            
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
    
    # --- Reactives for File Inputs ---
    bed_file_path <- reactive({ parsePath_gwas(input$bed_file_path) })
    gwas_table_path <- reactive({ parsePath_gwas(input$gwas_table_path) })
    annotation_table_path <- reactive({ parsePath_gwas(input$annotation_table_path) })
    plink_path <- reactive({ parsePath_gwas(input$plink_path) })
    
    # --- TBI status reactive for the VCF file ---
    tbi_status <- reactiveVal(FALSE)
    
    # --- File Chooser Setup ---
    rootDir <- c(Home = fs::path_home())
    # make this example path to reduce click
    # rootDir <- "inst/extdata"
    shinyFileChoose(input, "bed_file_path", roots = rootDir, defaultPath = "Projects/panvaR/inst/extdata", session = session)
    shinyFileChoose(input, "gwas_table_path", roots = rootDir, defaultPath = "Projects/panvaR/inst/extdata", session = session)
    shinyFileChoose(input, "annotation_table_path", roots = rootDir, defaultPath = "Projects/panvaR/inst/extdata", session = session)
    shinyFileChoose(input, "plink_path", roots = rootDir, defaultPath = "Projects/panvaR/inst/extdata", session = session)
    
    # --- Display Selected File Paths ---
    output$bed_file_path_results <- renderText({
      path <- bed_file_path()
      if (!is.null(path)) paste("Bed file:", basename(path))
    })
    output$gwas_table_path_results <- renderText({
      path <- gwas_table_path()
      if (!is.null(path)) paste("GWAS results file:", basename(path))
    })
    output$anno_table_path_results <- renderText({
      path <- annotation_table_path()
      if (!is.null(path)) paste("Annotation table file:", basename(path))
    })
    
    
    # --- Check for TBI file when VCF path changes ---
    # observeEvent(vcf_file_path(), {
    #   path <- vcf_file_path()
    #   if (is.null(path) || !file.exists(path)) {
    #     tbi_status(FALSE)
    #     return()
    #   }
    #   if (!endsWith(path, ".gz")) {
    #     showNotification("Warning: The selected VCF file is not gzipped (.gz). Indexing is required and only works on bgzip-compressed files.", type = "warning", duration=8)
    #     tbi_status(FALSE)
    #     return()
    #   }
    #   
    #   has_tbi <- file.exists(paste0(path, ".tbi"))
    #   tbi_status(has_tbi)
    #   if (!has_tbi) {
    #     showModal(modalDialog(
    #       title = "VCF Index File (.tbi) Missing",
    #       "The selected VCF file is missing its index (.tbi), which is required. Would you like to create it now? This may take a moment.",
    #       footer = tagList(
    #         actionButton(ns("generate_tbi_yes"), "Yes, Create Index", class = "btn-primary"),
    #         modalButton("Cancel")
    #       )
    #     ))
    #   }
    # })
    
    # --- Handle TBI Generation ---
    # observeEvent(input$generate_tbi_yes, {
    #   removeModal()
    #   vcf_path <- vcf_file_path()
    #   req(vcf_path)
    #   
    #   withProgress(message = 'Generating VCF index...', value = 0.5, {
    #     tryCatch({
    #       generate_tbi_file(vcf_path) # Assumes generate_tbi_file is available from the package
    #       tbi_status(TRUE)
    #       showNotification("VCF index file (.tbi) created successfully!", type = "message")
    #     }, error = function(e) {
    #       showNotification(paste("Error creating index:", e$message, ". Please ensure bcftools/tabix is installed and in your PATH."), type = "error", duration = 10)
    #       tbi_status(FALSE)
    #     })
    #   })
    # })
    
    # --- Reactive for Input Validation ---
    all_inputs_valid <- reactive({
      !is.null(bed_file_path()) && !is.null(gwas_table_path()) &&
        file.exists(bed_file_path()) && file.exists(gwas_table_path()) &&
        is.numeric(input$r2_threshold) && input$r2_threshold >= 0 && input$r2_threshold <= 1 &&
        is.numeric(input$window_span) && input$window_span >= 0
    })
    
    # --- UI for Table of inputs ---
    output$input_status_update <- renderUI({
      missing_inputs <- character(0)
      current_values <- list()
      
      if (!is.null(bed_file_path()) && file.exists(bed_file_path())) {
        current_values$`Bed File` <- basename(bed_file_path())
      } else {
        missing_inputs <- c(missing_inputs, "A valid bed file")
      }
      
      if (!is.null(gwas_table_path()) && file.exists(gwas_table_path())) {
        current_values$`GWAS Table` <- basename(gwas_table_path())
      } else {
        missing_inputs <- c(missing_inputs, "A valid GWAS results file")
      }
      
      if (!is.null(annotation_table_path()) && file.exists(annotation_table_path())) {
        current_values$`Annotation Table` <- basename(annotation_table_path())
      } else {
        missing_inputs <- c(missing_inputs, "An annotation table")
      }
      
      tag_snps <- clean_snp_tags_gwas(input$tag_snps)
      current_values$`Tag SNPs` <- if (!is.null(tag_snps)) paste(tag_snps, collapse = ", ") else "Auto-infer from GWAS results"
      current_values$`R² Threshold` <- input$r2_threshold
      current_values$`Window (kb)` <- format(input$window_span, big.mark = ",")
      
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
        
        withProgress(message = 'Running make_panvar_tables()', value = 0, {
            # get bedfile path and filename 
            bedfullpath <- bed_file_path()
            bedfile <- str_replace(basename(bedfullpath), "\\.bed", "")
            beddir <- dirname(bedfullpath)
          
            # read in files 
            gwas.df <- data.table::fread(gwas_table_path(), data.table = F)
            anno.df <- data.table::fread(annotation_table_path(), data.table = F)
            
            # run the function
            results <- 
            make_panvar_tables(gwas.res = gwas.df,
                               annotation.table = anno.df,
                               geno.bed.filename = bedfile,
                               geno.bed.directory = beddir,
                               tag.snp = input$tag_snp,
                               pvals.in.log = input$pvals_in_log,
                               window = input$window_span,
                               snp.to.gene.buffer = input$snp.to.gene.buffer,
                               compute.scores = input$compute_scores)
            
            # Call the modified panvar_func
            # panvar_func(
            #     vcf_file_path = vcf_file_path(),
            #     gwas_data = gwas_table_path(),
            #     phenotype_data = NULL, # Explicitly pass NULL for phenotype
            #     tag_snps = clean_snp_tags_gwas(input$tag_snps),
            #     r2_threshold = input$r2_threshold,
            #     maf = input$maf,
            #     missing_rate = input$missing_rate,
            #     window = input$window_span,
            #     all.impacts = input$all_impacts
            # )
            
          })
          
          if (!is.null(results)) {
            shared$analysis_results <- results
            showNotification("Analysis completed successfully!", type = "message")
          } else {
            shared$analysis_results <- NULL
          }
        
        # save window span for plotting 
        shared$window_span <- input$window_span
        
      } else {
        showNotification("Please provide all required and valid inputs before running the analysis.", type = "error")
      }
    })
  })
}