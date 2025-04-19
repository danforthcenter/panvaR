# Output_dashboard_panvar.r
# Ensure necessary libraries are loaded in the main app script (panvar_shiny_main.r)
# library(shiny)
# library(shinyFiles)
# library(tidyverse) # Loaded via data.table usually
# library(data.table)
# library(plotly)
# library(DT)
# library(shinyBS)
# library(shinyjs)

# --- Helper functions ---

# Enhanced load and filter module with NULL handling and specific tag SNP BP logic
load_and_filter_module <- function(my_data, input) {
  # Check if data is NULL or empty
  if (is.null(my_data) || nrow(my_data) == 0) {
    return(NULL)
  }
  
  # Safely get tag_snp BP from the *current* data subset
  tag_snp_row <- my_data %>% filter(Type == "tag_snp") %>% head(1) # Get the first tag snp row for this subset
  tag_snp_bp <- if(nrow(tag_snp_row) > 0) as.numeric(tag_snp_row$BP) else NULL
  
  # If no tag_snp found in this specific dataset, warn and disable BP filter for it
  if (is.null(tag_snp_bp)) {
    # This warning might be noisy if the tag SNP itself was filtered out
    # warning("No tag_snp BP found in the current data subset, BP range filtering disabled for this tab.")
  }
  
  # Safely apply filters with error handling
  filtered_data <- tryCatch({
    result <- my_data
    
    # Apply numeric filters if input exists
    if (!is.null(input$LD_min) && !is.null(input$LD_max)) {
      # Handle potential NA in LD column during filtering
      result <- result %>%
        filter(is.na(LD) | (LD >= as.numeric(input$LD_min) & LD <= as.numeric(input$LD_max)))
    }
    
    # Apply BP filter only if tag_snp_bp was successfully found FOR THIS SUBSET
    if (!is.null(tag_snp_bp) && !is.null(input$BP_LHS) && !is.null(input$BP_RHS)) {
      # Ensure BP is numeric before filtering
      if(is.numeric(result$BP)){
        result <- result %>%
          filter(BP >= as.numeric(tag_snp_bp - input$BP_LHS) &
                   BP <= as.numeric(tag_snp_bp + input$BP_RHS))
      } else {
        warning("BP column is not numeric, cannot apply BP range filter.")
      }
    }
    
    
    if (!is.null(input$pvalue_min) && !is.null(input$pvalue_max)) {
      # Ensure Pvalues is numeric before filtering
      if(is.numeric(result$Pvalues)){
        result <- result %>%
          filter(Pvalues >= as.numeric(input$pvalue_min) &
                   Pvalues <= as.numeric(input$pvalue_max))
      } else {
        warning("Pvalues column is not numeric, cannot apply P-value range filter.")
      }
    }
    
    # Apply categorical filters only if selections exist
    if (!is.null(input$selected_genes) && length(input$selected_genes) > 0) {
      result <- result %>% filter(GENE %in% input$selected_genes)
    }
    
    if (!is.null(input$selected_effect_types) && length(input$selected_effect_types) > 0) {
      result <- result %>% filter(EFFECT %in% input$selected_effect_types)
    }
    
    if (!is.null(input$selected_amino_acid) && length(input$selected_amino_acid) > 0) {
      result <- result %>% filter(AA %in% input$selected_amino_acid)
    }
    
    if (!is.null(input$selected_REF_types) && length(input$selected_REF_types) > 0) {
      result <- result %>% filter(REF %in% input$selected_REF_types)
    }
    
    if (!is.null(input$selected_ALT_types) && length(input$selected_ALT_types) > 0) {
      result <- result %>% filter(ALT %in% input$selected_ALT_types)
    }
    
    result
    
  }, error = function(e) {
    warning("Error filtering data: ", e$message)
    return(NULL) # Return NULL if filtering fails
  })
  
  return(filtered_data)
}


# Rijan: Plotly object generator - handles NA LD, uses correct color scale application
panvar_plotly_function <- function(panvar_results_table, nrows_in_gwas = NULL, pvalue_threshold = 0.05, point_size = 3, alpha_base = 0.7){
  
  # Handle NULL or empty input table gracefully
  if (is.null(panvar_results_table) || nrow(panvar_results_table) == 0) {
    return(
      plot_ly() %>%
        layout(title = "No Data to Display", xaxis = list(visible = FALSE), yaxis = list(visible = FALSE)) %>%
        add_annotations(
          text = "No data available based on current filters for this tag SNP.",
          x = 0.5, y = 0.5, xref = "paper", yref = "paper", showarrow = FALSE, font = list(size=14)
        )
    )
  }
  
  # Check if essential columns exist before plotting
  required_plot_cols <- c("BP", "Pvalues", "IMPACT", "LD", "GENE", "Type")
  missing_plot_cols <- setdiff(required_plot_cols, names(panvar_results_table))
  if (length(missing_plot_cols) > 0) {
    warning(paste("Plotting function is missing required columns:", paste(missing_plot_cols, collapse=", ")))
    return(
      plot_ly() %>%
        layout(title = "Plotting Error", xaxis = list(visible = FALSE), yaxis = list(visible = FALSE)) %>%
        add_annotations(
          text = paste("Cannot generate plot. Missing columns:", paste(missing_plot_cols, collapse=", ")),
          x = 0.5, y = 0.5, xref = "paper", yref = "paper", showarrow = FALSE, font = list(size=14)
        )
    )
  }
  
  # Ensure required numeric columns exist and are numeric
  if (!is.numeric(panvar_results_table$BP) || !is.numeric(panvar_results_table$Pvalues) || !("LD" %in% names(panvar_results_table)) || !is.numeric(panvar_results_table$LD)) {
    warning("One or more required numeric columns (BP, Pvalues, LD) are missing or not numeric.")
    return(
      plot_ly() %>%
        layout(title = "Plotting Error", xaxis = list(visible = FALSE), yaxis = list(visible = FALSE)) %>%
        add_annotations(
          text = "Cannot generate plot. BP, Pvalues, or LD column is missing or not numeric.",
          x = 0.5, y = 0.5, xref = "paper", yref = "paper", showarrow = FALSE, font = list(size=14)
        )
    )
  }
  
  
  # --- Check for valid LD values for color mapping ---
  has_valid_ld <- any(!is.na(panvar_results_table$LD))
  # --- End Check ---
  
  # Rijan: What is the tag_snp
  tag_df <- panvar_results_table %>%
    filter(Type == 'tag_snp')
  tag_snp <- if(nrow(tag_df) > 0) tag_df %>% pull(BP) %>% unique() %>% head(1) else NULL
  
  
  # Rijan: Bonferroni correction setup
  if(is.null(nrows_in_gwas)){
    print("You did not supply a value for the number of tests that were in the GWAS - a place holder value will be used. This is not ideal.")
    nrows_in_gwas <- 5e6 # Use a default if not provided
  }
  bonf.cor <- -log10(pvalue_threshold / nrows_in_gwas)
  
  
  # --- Create Plotly object conditionally based on LD values ---
  if (has_valid_ld) {
    # CASE 1: Valid LD values exist, map color to LD and use colorscale in marker
    panvar_plotly <- plot_ly(
      panvar_results_table,
      x = ~BP,
      y = ~Pvalues,
      color = ~LD,        # Map color variable
      # colors argument removed from here
      type = 'scatter',
      mode = 'markers',
      symbol = ~IMPACT,
      marker = list(
        size = 10,
        colorscale = "Viridis", # <<< MODIFIED: Specify colorscale within marker
        colorbar = list(title = "LD (R²)") # Include color bar
      ),
      hoverinfo = 'text',
      text = ~paste('BP:', BP, '<br>P-value:', round(Pvalues, 2), '<br>LD:', round(LD, 2), '<br>Impact:', IMPACT, '<br>Gene:', GENE),
      showlegend = TRUE,
      name = ~IMPACT
    )
  } else {
    # CASE 2: No valid LD values (all NA or empty), do not map color to LD
    panvar_plotly <- plot_ly(
      panvar_results_table,
      x = ~BP,
      y = ~Pvalues,
      # NO 'color' mapping argument here
      type = 'scatter',
      mode = 'markers',
      symbol = ~IMPACT,
      marker = list(
        size = 10,
        color = 'gray' # Assign a default fixed color
        # NO 'colorscale' or 'colorbar' here
      ),
      hoverinfo = 'text',
      text = ~paste('BP:', BP, '<br>P-value:', round(Pvalues, 2), '<br>LD: NA', '<br>Impact:', IMPACT, '<br>Gene:', GENE), # Indicate LD is NA
      showlegend = TRUE,
      name = ~IMPACT
    )
  }
  # --- End Conditional Plot Creation ---
  
  
  # --- Add common layout elements and shapes ---
  panvar_plotly <- layout(panvar_plotly, title = "Interactive PanvaR Plot",
                          xaxis = list(title = "Position (BP)", titlefont = list(size = 14)),
                          yaxis = list(title = "-log<sub>10</sub>(P-value)", titlefont = list(size = 14)),
                          legend = list(title = list(text = '<b>Impact</b>'), orientation = "h", y = -0.2)
  )
  
  # Functions for lines
  hline <- function(y = 0, color = "blue") list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = y, y1 = y, yref = "y", line = list(color = color, dash = "dash"))
  vline <- function(x = 0, color = "red") list(type = "line", y0 = 0, y1 = 1, yref = "paper", x0 = x, x1 = x, xref = "x", line = list(color = color, dash = "dot"))
  
  # Add shapes
  shapes_list <- list(hline(bonf.cor))
  if (!is.null(tag_snp) && is.finite(tag_snp)) {
    shapes_list <- c(shapes_list, list(vline(tag_snp)))
  }
  panvar_plotly <- layout(panvar_plotly, shapes = shapes_list)
  # --- End Common Layout ---
  
  return(panvar_plotly)
  
}

# --- UI Modification ---
output_dashboard_UI <- function(id) {
  ns <- NS(id)
  tagList(
    # Sidebar remains the same structure, but IDs are namespaced
    sidebarLayout(
      sidebarPanel(
        width = 3,
        # --- Add a placeholder for dynamic filter updates message ---
        uiOutput(ns("filter_update_status")), # Placeholder for messages
        # --- End placeholder ---
        
        # Data source selection with radio buttons
        radioButtons(
          ns("data_source"),
          "Select Data Source:",
          choices = c(
            "Load from dynamic analysis" = "dynamic",
            "Load from file" = "file"
          ),
          selected = "dynamic"
        ),
        
        # Conditional panel for file selection - only shows when "file" is selected
        conditionalPanel(
          condition = sprintf("input['%s'] === 'file'", ns("data_source")),
          div(
            style = "display: flex; align-items: center; gap: 10px;",
            shinyFilesButton(
              ns("Pre_existing_panvaR_results_path"),
              "Select Pre-existing panvaR results",
              "Please select a file",
              multiple = FALSE # Keep as FALSE, expecting a single RDS/TSV file containing the list or single result
            ),
            tags$span(
              id = ns("Pre_existing_panvaR_results_path_tooltip"),
              icon("question-circle"),
              style = "color: green;"
            )
          ),
          bsTooltip(
            id = ns("Pre_existing_panvaR_results_path_tooltip"),
            title = "Provide path to a pre-existing panvaR result (can be an .rds file containing a list for multiple tag SNPs, or a .tsv/.csv for a single result).",
            placement = "right",
            trigger = "hover"
          ),
          textOutput(ns("Pre_existing_panvaR_results"))
        ),
        
        # Rest of the filtering controls - conditionally enabled
        # Wrap filters in a div to easily enable/disable with shinyjs
        div(id = ns("filter_controls"), # Add an ID for shinyjs targeting
            # Select Genes
            div(
              style = "display: flex; align-items: center; gap: 10px;",
              selectizeInput(
                ns("selected_genes"),
                "Select Genes:",
                choices = NULL,
                multiple = TRUE
              ),
              tags$span(
                id = ns("selected_genes_tooltip"),
                icon("question-circle"),
                style = "color: green;"
              )
            ),
            bsTooltip(
              id = ns("selected_genes_tooltip"),
              title = "Filter results for the selected tag SNP by Gene.",
              placement = "right",
              trigger = "hover"
            ),
            
            # Select Effect types
            div(
              style = "display: flex; align-items: center; gap: 10px;",
              selectizeInput(
                ns("selected_effect_types"),
                "Select Effect types:",
                choices = NULL,
                multiple = TRUE
              ),
              span(
                icon("question-circle", style = "color: green;"),
                id = ns("selected_effect_types_tooltip") # Unique ID for tooltip
              )
            ),
            bsTooltip(
              id = ns("selected_effect_types_tooltip"), # Match span ID
              title = "Filter results for the selected tag SNP by Effect type.",
              placement = "right",
              trigger = "hover"
            ),
            
            
            # Select Amino Acid Changes
            div(
              style = "display: flex; align-items: center; gap: 10px;",
              selectizeInput(
                ns("selected_amino_acid"),
                "Select Amino Acid Changes:",
                choices = NULL,
                multiple = TRUE
              ),
              span(
                icon("question-circle", style = "color: green;"),
                id = ns("selected_amino_acid_tooltip") # Unique ID
              )
            ),
            bsTooltip(
              id = ns("selected_amino_acid_tooltip"), # Match span ID
              title = "Filter results for the selected tag SNP by Amino Acid change.",
              placement = "right",
              trigger = "hover"
            ),
            
            # Select ALT type
            div(
              style = "display: flex; align-items: center; gap: 10px;",
              selectizeInput(
                ns("selected_ALT_types"),
                "Select ALT type:",
                choices = NULL,
                multiple = TRUE
              ),
              span(
                icon("question-circle", style = "color: green;"),
                id = ns("selected_ALT_types_tooltip") # Unique ID
              )
            ),
            bsTooltip(
              id = ns("selected_ALT_types_tooltip"), # Match span ID
              title = "Filter results for the selected tag SNP by ALT allele.",
              placement = "right",
              trigger = "hover"
            ),
            
            # Select REF type
            div(
              style = "display: flex; align-items: center; gap: 10px;",
              selectizeInput(
                ns("selected_REF_types"),
                "Select REF type :",
                choices = NULL,
                multiple = TRUE
              ),
              span(
                icon("question-circle", style = "color: green;"),
                id = ns("selected_REF_types_tooltip") # Unique ID
              )
            ),
            bsTooltip(
              id = ns("selected_REF_types_tooltip"), # Match span ID
              title = "Filter results for the selected tag SNP by REF allele.",
              placement = "right",
              trigger = "hover"
            ),
            
            
            # LD range
            div(
              style = "display: flex; align-items: center; gap: 10px;",
              div(
                style = "flex: 1;",
                fluidRow(
                  column(6,
                         numericInput(ns("LD_min"),
                                      "LD Min:",
                                      value = 0,
                                      min = 0,
                                      max = 1,
                                      step = 0.01)),
                  column(6,
                         numericInput(ns("LD_max"),
                                      "LD Max:",
                                      value = 1,
                                      min = 0,
                                      max = 1,
                                      step = 0.01))
                )
              ),
              span(
                style = "flex-shrink: 0;",
                icon("question-circle", style = "color: green;"),
                id = ns("LD_range_tooltip") # Unique ID
              )
            ),
            bsTooltip(
              id = ns("LD_range_tooltip"), # Match span ID
              title = "Filter results for the selected tag SNP by LD (R²) range.",
              placement = "right",
              trigger = "hover"
            ),
            
            
            # Base position range
            div(
              style = "display: flex; align-items: center; gap: 10px;",
              div(
                style = "flex: 1;",
                fluidRow(
                  column(6,
                         numericInput(ns("BP_LHS"),
                                      "BP range LHS of tag_SNP:", # Label clarity
                                      value = 0,
                                      min = 0,
                                      # Max will be updated dynamically
                                      step = 1)
                  ),
                  column(6,
                         numericInput(ns("BP_RHS"),
                                      "BP range RHS of tag_snp:", # Label clarity
                                      value = 0,
                                      min = 0,
                                      # Max will be updated dynamically
                                      step = 1)
                  )
                )
              ),
              span(
                style = "flex-shrink: 0;",
                icon("question-circle", style = "color: green;"),
                id = ns("BP_range_tooltip") # Unique ID
              )
            ),
            bsTooltip(
              id = ns("BP_range_tooltip"), # Match span ID
              title = "Filter results by BP range relative to the *selected* tag SNP. LHS = Max distance left, RHS = Max distance right.",
              placement = "right",
              trigger = "hover"
            ),
            
            
            # P-value range
            div(
              style = "display: flex; align-items: center; gap: 10px;",
              div(
                style = "flex: 1;",
                fluidRow(
                  column(6,
                         numericInput(ns("pvalue_min"),
                                      "P-value Min (-log10):", # Label clarity
                                      value = 0,
                                      min = 0)), # Min P-value usually 0
                  column(6,
                         numericInput(ns("pvalue_max"),
                                      "P-value Max (-log10):", # Label clarity
                                      value = 30)) # Default max
                )
              ),
              span(
                style = "flex-shrink: 0;",
                icon("question-circle", style = "color: green;"),
                id = ns("Pvalue_range_tooltip") # Unique ID
              )
            ),
            bsTooltip(
              id = ns("Pvalue_range_tooltip"), # Match span ID
              title = "Filter by -log10(P-value) range for the selected tag SNP.",
              placement = "right",
              trigger = "hover"
            ),
            
            
            # Bonferroni correction inputs
            div(
              style = "display: flex; align-items: center; gap: 10px;",
              div(
                style = "flex: 1;",
                fluidRow(
                  column(6,
                         numericInput(ns("pvalue_threshold"),
                                      "Sig. Threshold (alpha):", # Label clarity
                                      value = 0.05, min = 0, max = 1, step = 0.01)),
                  column(6,
                         numericInput(ns("total_rows"),
                                      "Total GWAS Tests:", # Label clarity
                                      value = 1e6, min = 1)) # Default, min 1 test
                )
              ),
              span(
                style = "flex-shrink: 0;",
                icon("question-circle", style = "color: green;"),
                id = ns("Bonferroni_correction_tooltip") # Unique ID
              )
            ),
            bsTooltip(
              id = ns("Bonferroni_correction_tooltip"), # Match span ID
              title = "Parameters for calculating the Bonferroni correction line on the plot (alpha / total tests). Applied to the plot for the selected tag SNP.",
              placement = "right",
              trigger = "hover"
            )
        ) # End div for filter_controls
        
      ), # End sidebarPanel
      
      # --- Modified mainPanel to use uiOutput ---
      mainPanel(
        # This will hold the dynamically generated tabsetPanel or direct output
        uiOutput(ns("results_output_area"))
      )
      # --- End Modification ---
    ) # End sidebarLayout
  ) # End tagList
}
# --- END UI ---


# --- START SERVER ---
output_dashboard_Server <- function(id, shared) {
  moduleServer(
    id,
    function(input, output, session) {
      ns <- session$ns # Ensure ns is available
      
      # Reactive value to store the raw results (list or single)
      raw_results <- reactiveVal(NULL)
      # Reactive value to store the names of the results (tag SNPs) if it's a list
      result_names <- reactiveVal(NULL)
      # Removed active_data_table - will derive active data within reactives/observers
      
      # --- File Chooser Setup ---
      shinyFileChoose(
        input,
        "Pre_existing_panvaR_results_path",
        roots = c(Home = fs::path_home()),
        session = session,
        filetypes = c("rds", "tsv", "csv") # Allow RDS for lists, tsv/csv for single tables
      )
      
      # Display chosen file path
      output$Pre_existing_panvaR_results <- renderText({
        req(input$Pre_existing_panvaR_results_path)
        file_info <- shinyFiles::parseFilePaths(c(Home = fs::path_home()), input$Pre_existing_panvaR_results_path)
        if (nrow(file_info) > 0) {
          paste("Selected file:", file_info$name[1])
        } else {
          "No file selected."
        }
      })
      
      # --- Observe Data Source and Load Data ---
      observe({
        source_type <- input$data_source
        results_to_store <- NULL
        names_to_store <- NULL
        
        if (source_type == "dynamic") {
          # Use results directly from the shared reactive value
          if (!is.null(shared$analysis_results)) {
            results_to_store <- shared$analysis_results
            # Try to get names if it's a list (might need adjustment based on how names are set in panvar_func)
            # Check if it's a list BUT not a data.frame (data.frames are lists) and not the single result structure
            if (is.list(results_to_store) && !is.data.frame(results_to_store) && !all(c("plot", "table") %in% names(results_to_store))) {
              names_to_store <- names(results_to_store)
              # If names are NULL, generate default names based on index
              if (is.null(names_to_store)) {
                names_to_store <- paste("Tag SNP", seq_along(results_to_store))
                # Optionally try to extract from the table if names are missing
                # This assumes the tag snp info is consistently available
                possible_names <- sapply(results_to_store, function(res) {
                  if (is.list(res) && "table" %in% names(res) && is.data.frame(res$table)) {
                    tag_row <- res$table %>% filter(Type == "tag_snp") %>% head(1)
                    if (nrow(tag_row) > 0) {
                      return(paste0(tag_row$CHROM, ":", tag_row$BP))
                    }
                  }
                  return(NA_character_) # Return NA if name cannot be extracted
                })
                # Use extracted names if available, otherwise stick to default index names
                valid_possible_names <- !is.na(possible_names)
                if(any(valid_possible_names)) names_to_store[valid_possible_names] <- possible_names[valid_possible_names]
                
              }
            } else if (is.list(results_to_store) && all(c("plot", "table") %in% names(results_to_store))) {
              # Handle the case where dynamic result is a single item but not in a list of length 1
              # We can try to extract the tag SNP name from the table
              if (is.data.frame(results_to_store$table)) {
                tag_row <- results_to_store$table %>% filter(Type == "tag_snp") %>% head(1)
                if (nrow(tag_row) > 0) {
                  names_to_store <- paste0(tag_row$CHROM, ":", tag_row$BP)
                } else {
                  names_to_store <- "Single Result"
                }
              } else {
                names_to_store <- "Single Result"
              }
              
            }
          } else {
            showNotification("No dynamic analysis results available yet.", type = "warning")
          }
        } else { # source_type == "file"
          req(input$Pre_existing_panvaR_results_path)
          file_info <- shinyFiles::parseFilePaths(c(Home = fs::path_home()), input$Pre_existing_panvaR_results_path)
          
          if (nrow(file_info) > 0) {
            file_path <- file_info$datapath[1]
            tryCatch({
              if (endsWith(tolower(file_path), ".rds")) {
                # Load RDS file, expecting a list or single result object
                loaded_object <- readRDS(file_path)
                # Check if it's a list of results (list of lists, where each inner list has plot/table)
                if (is.list(loaded_object) && !is.data.frame(loaded_object) && length(loaded_object)>0 && is.list(loaded_object[[1]]) && all(c("plot", "table") %in% names(loaded_object[[1]]))) {
                  # Assume it's a list of results
                  results_to_store <- loaded_object
                  names_to_store <- names(results_to_store)
                  if (is.null(names_to_store)) {
                    # Try to extract names from tables
                    possible_names <- sapply(results_to_store, function(res) {
                      if (is.list(res) && "table" %in% names(res) && is.data.frame(res$table)) {
                        tag_row <- res$table %>% filter(Type == "tag_snp") %>% head(1)
                        if (nrow(tag_row) > 0) return(paste0(tag_row$CHROM, ":", tag_row$BP))
                      }
                      return(NA_character_)
                    })
                    valid_possible_names <- !is.na(possible_names)
                    if(any(valid_possible_names)) {
                      names_to_store <- possible_names
                      # Fill remaining NAs with default names if needed
                      names_to_store[!valid_possible_names] <- paste("Result", which(!valid_possible_names))
                    } else {
                      names_to_store <- paste("Result", seq_along(results_to_store))
                    }
                  }
                } else if (is.list(loaded_object) && all(c("plot", "table") %in% names(loaded_object))) {
                  # Assume it's a single result structure
                  results_to_store <- loaded_object
                  # Try to get name for the single result
                  if (is.data.frame(results_to_store$table)) {
                    tag_row <- results_to_store$table %>% filter(Type == "tag_snp") %>% head(1)
                    if (nrow(tag_row) > 0) names_to_store <- paste0(tag_row$CHROM, ":", tag_row$BP) else names_to_store <- "Single Result"
                  } else {
                    names_to_store <- "Single Result"
                  }
                } else {
                  stop("RDS file does not contain a valid panvaR result list or single result structure.")
                }
                
              } else if (endsWith(tolower(file_path), ".tsv") || endsWith(tolower(file_path), ".csv")) {
                # Load delimited file, assuming it's a single result table
                single_table <- data.table::fread(file_path)
                # Basic validation for required columns in the single table
                required_cols_single <- c("CHROM", "BP", "Pvalues", "LD", "Type", "IMPACT", "GENE", "EFFECT", "AA", "REF", "ALT", "final_weight") # Add all expected columns
                missing_cols_single <- setdiff(required_cols_single, names(single_table))
                if(length(missing_cols_single) > 0) {
                  stop(paste("Loaded table file is missing required columns:", paste(missing_cols_single, collapse=", ")))
                }
                # We can't reconstruct the original plot object, so set plot to NULL
                results_to_store <- list(plot = NULL, table = single_table)
                # Try to get name for the single result
                tag_row <- single_table %>% filter(Type == "tag_snp") %>% head(1)
                if (nrow(tag_row) > 0) names_to_store <- paste0(tag_row$CHROM, ":", tag_row$BP) else names_to_store <- "Single Result File"
                showNotification("Loaded table from file. Plot generation from file is not supported.", type = "info", duration=7)
                
              } else {
                stop("Unsupported file type. Please select an .rds, .tsv, or .csv file.")
              }
              
              # Check if results are empty after loading
              is_empty_result <- FALSE
              if(is.null(results_to_store)) {
                is_empty_result <- TRUE
              } else if (is.list(results_to_store) && !is.data.frame(results_to_store) && length(results_to_store) == 0) {
                is_empty_result <- TRUE
              } else if (is.list(results_to_store) && all(c("plot", "table") %in% names(results_to_store)) && (is.null(results_to_store$table) || nrow(results_to_store$table) == 0) ){
                is_empty_result <- TRUE
              }
              
              if (is_empty_result) {
                showNotification("Loaded file contains no valid data.", type = "warning")
                results_to_store <- NULL
                names_to_store <- NULL
              }
              
            }, error = function(e) {
              showNotification(paste("Error loading or processing file:", e$message), type = "error")
              results_to_store <- NULL
              names_to_store <- NULL
            })
          } else {
            # This case should ideally not be reached if req() is used, but added for safety
            results_to_store <- NULL
            names_to_store <- NULL
          }
        }
        raw_results(results_to_store)
        result_names(names_to_store)
        
        # Enable/disable filter controls based on whether data loaded
        if (is.null(results_to_store)) {
          shinyjs::disable("filter_controls")
          output$filter_update_status <- renderUI({ # Clear status message
            p(em("Load data or run analysis to enable filters."), style="color:gray;")
          })
        } else {
          shinyjs::enable("filter_controls")
          # Trigger initial filter update explicitly after enabling
          # Use a slight delay to ensure UI is ready
          shinyjs::delay(100, {
            # This needs to trigger the observeEvent(active_tab_data()...)
            # One way is to invalidate the active_tab_data reactive, but safer
            # might be to manually call the update logic here or trigger its dependency.
            # For now, let's rely on the automatic triggering.
            # Resetting a dummy reactiveVal could also work if needed.
            output$filter_update_status <- renderUI({ # Clear status message
              p(em("Filters enabled. Select tab to update choices."), style="color:gray;")
            })
          })
        }
        
      })
      
      
      # --- Dynamically generate UI (Tabs or Direct Output) ---
      output$results_output_area <- renderUI({
        results <- raw_results()
        names_res <- result_names()
        
        if (is.null(results)) {
          return(div(style = "padding: 20px; text-align: center;", h4("No analysis results loaded or available.")))
        }
        
        # Check if it's a list representing multiple results
        is_multiple_results <- is.list(results) && !is.data.frame(results) && length(results) > 1 &&
          !all(c("plot", "table") %in% names(results)) # Differentiates from single result structure
        
        if (is_multiple_results) {
          # Ensure names are available and match length
          if(is.null(names_res) || length(names_res) != length(results)){
            names_res <- paste("Result", seq_along(results)) # Fallback names
            warning("Result names missing or incorrect length for tabs, using default names.")
          }
          # Create tabs
          tabs <- lapply(seq_along(results), function(i) {
            tabPanel(title = names_res[[i]], # Use tag SNP or default name as title
                     value = ns(paste0("tab_", i)), # Unique value for each tab using index
                     plotlyOutput(ns(paste0("plot_", i))), # Unique plot output ID
                     DTOutput(ns(paste0("table_", i)))      # Unique table output ID
            )
          })
          do.call(tabsetPanel, c(tabs, id = ns("resultTabs"))) # Use ns() for tabsetPanel ID
        } else {
          # Handle single result (could be direct structure or list of length 1)
          single_result_data <- NULL
          single_result_name <- "Single Result" # Default name
          
          if (is.list(results) && !is.data.frame(results) && length(results) == 1) {
            single_result_data <- results[[1]]
            if (!is.null(names_res) && length(names_res)==1) single_result_name <- names_res[[1]]
          } else if (is.list(results) && all(c("plot", "table") %in% names(results))) {
            single_result_data <- results
            if (!is.null(names_res) && length(names_res)==1) single_result_name <- names_res[[1]]
          }
          
          if (!is.null(single_result_data) && is.list(single_result_data) && "table" %in% names(single_result_data)) {
            # Display single result directly (or could wrap in a single tab)
            # Option 1: Direct display
            # tagList(
            #     h4(single_result_name), # Display the name/title
            #     plotlyOutput(ns("plot_single")),
            #     DTOutput(ns("table_single"))
            # )
            # Option 2: Wrap in a single-tab tabsetPanel for consistency
            do.call(tabsetPanel, list(
              id = ns("resultTabs"), # Still give it an ID
              tabPanel(title = single_result_name,
                       value = ns("tab_1"), # Use index 1 for consistency
                       plotlyOutput(ns("plot_1")), # Use indexed IDs
                       DTOutput(ns("table_1"))       # Use indexed IDs
              )
            ))
            
          } else {
            # Handle invalid structure or empty single result
            div(style = "padding: 20px; text-align: center;", h4("Loaded result has an unexpected format or is empty."))
          }
        }
      })
      
      # --- Reactive for Currently Selected Tab Data Table ---
      active_tab_data <- reactive({
        results <- raw_results()
        # Use input$resultTabs which now exists for both single and multiple tabs
        selected_tab_value <- input$resultTabs
        
        req(results) # Need results to proceed
        req(selected_tab_value) # Need a selected tab value
        
        # Extract index from tab value (e.g., "module2-tab_1" -> 1)
        match_val <- regmatches(selected_tab_value, regexpr("\\d+$", selected_tab_value))
        if(length(match_val) == 0) return(NULL)
        selected_tab_index <- as.integer(match_val)
        
        
        # Determine if results represent multiple entries or a single one
        is_multiple_results <- is.list(results) && !is.data.frame(results) && length(results) > 1 &&
          !all(c("plot", "table") %in% names(results))
        
        current_item <- NULL
        if (is_multiple_results) {
          # Select the item from the list based on index
          if (selected_tab_index > 0 && selected_tab_index <= length(results)) {
            current_item <- results[[selected_tab_index]]
          }
        } else {
          # It's a single result (either direct or list of 1)
          if (is.list(results) && !is.data.frame(results) && length(results) == 1) {
            current_item <- results[[1]]
          } else if (is.list(results) && all(c("plot", "table") %in% names(results))) {
            current_item <- results
          }
        }
        
        # Validate the structure and return the table
        if(is.list(current_item) && "table" %in% names(current_item) && is.data.frame(current_item$table)){
          return(current_item$table)
        } else {
          warning(paste("Result item for selected tab", selected_tab_value, "has invalid structure or no table."))
          return(NULL)
        }
      })
      
      # --- Reactive for Currently Selected Tab Plot Object ---
      # Needed if loading plot from file is not supported/desired
      active_tab_plot_object <- reactive({
        results <- raw_results()
        selected_tab_value <- input$resultTabs
        req(results, selected_tab_value)
        
        match_val <- regmatches(selected_tab_value, regexpr("\\d+$", selected_tab_value))
        if(length(match_val) == 0) return(NULL)
        selected_tab_index <- as.integer(match_val)
        
        is_multiple_results <- is.list(results) && !is.data.frame(results) && length(results) > 1 && !all(c("plot", "table") %in% names(results))
        
        current_item <- NULL
        if (is_multiple_results) {
          if (selected_tab_index > 0 && selected_tab_index <= length(results)) current_item <- results[[selected_tab_index]]
        } else {
          if (is.list(results) && !is.data.frame(results) && length(results) == 1) current_item <- results[[1]]
          else if (is.list(results) && all(c("plot", "table") %in% names(results))) current_item <- results
        }
        
        # Return the plot object if it exists, otherwise NULL
        if(is.list(current_item) && "plot" %in% names(current_item)){
          return(current_item$plot) # This might be NULL if loaded from TSV/CSV
        } else {
          return(NULL)
        }
      })
      
      
      # --- Update Filter Choices Based on Active Tab ---
      observeEvent(active_tab_data(), {
        current_data <- active_tab_data() # This is the data.frame for the active tab
        req(current_data) # Require data for the active tab
        
        # Show brief status message
        output$filter_update_status <- renderUI({
          tagList(
            hr(),
            p(em("Updating filter choices for the selected tab..."), style="color:gray; font-size: small;")
          )
        })
        
        # --- Start UI Update ---
        # Helper to safely update selectize, preserving selection if possible, handling NULLs
        safe_update_selectize <- function(inputId, choices, current_selection = NULL) {
          valid_choices <- NULL
          if (!is.null(choices) && length(choices) > 0) {
            valid_choices <- sort(unique(na.omit(choices)))
          }
          # Keep existing selection if it's still valid, otherwise clear
          selected_val <- intersect(current_selection, valid_choices)
          if (length(selected_val) == 0) selected_val <- character(0)
          
          updateSelectizeInput(session, inputId, choices = valid_choices, selected = selected_val, server = TRUE) # server=TRUE can improve performance for large lists
        }
        
        # Use existing input values as current_selection where appropriate
        safe_update_selectize("selected_genes", current_data$GENE, input$selected_genes)
        safe_update_selectize("selected_effect_types", current_data$EFFECT, input$selected_effect_types)
        safe_update_selectize("selected_amino_acid", current_data$AA, input$selected_amino_acid)
        safe_update_selectize("selected_REF_types", current_data$REF, input$selected_REF_types)
        safe_update_selectize("selected_ALT_types", current_data$ALT, input$selected_ALT_types)
        
        
        # --- Update Numeric Ranges based on Active Tab ---
        # BP range relative to tag SNP in the current tab
        tag_snp_bp_active <- tryCatch({
          current_data %>% filter(Type == "tag_snp") %>% pull(BP) %>% unique() %>% head(1) %>% as.numeric()
        }, error = function(e) NULL)
        
        # Update BP range sliders, resetting values to max range initially
        if(!is.null(tag_snp_bp_active) && !is.na(tag_snp_bp_active) && is.numeric(current_data$BP)) {
          min_BP_active <- min(current_data$BP, na.rm = TRUE)
          max_BP_active <- max(current_data$BP, na.rm = TRUE)
          
          if(is.finite(min_BP_active) && is.finite(max_BP_active)) {
            bp_LHS_max_active <- tag_snp_bp_active - min_BP_active
            bp_RHS_max_active <- max_BP_active - tag_snp_bp_active
            bp_LHS_max_active <- max(0, bp_LHS_max_active, na.rm = TRUE)
            bp_RHS_max_active <- max(0, bp_RHS_max_active, na.rm = TRUE)
            
            # Update max first, then value to avoid conflicts
            updateNumericInput(session, "BP_LHS", max = bp_LHS_max_active, value = bp_LHS_max_active)
            updateNumericInput(session, "BP_RHS", max = bp_RHS_max_active, value = bp_RHS_max_active)
          } else {
            # Handle cases where min/max BP aren't finite
            updateNumericInput(session, "BP_LHS", max = 0, value = 0)
            updateNumericInput(session, "BP_RHS", max = 0, value = 0)
          }
        } else {
          # Reset or disable BP range if no tag SNP found or BP not numeric
          updateNumericInput(session, "BP_LHS", value = 0, max = 0)
          updateNumericInput(session, "BP_RHS", value = 0, max = 0)
        }
        
        # P-value range for active tab - reset values initially
        if(is.numeric(current_data$Pvalues)) {
          min_pval_active <- min(current_data$Pvalues, na.rm = TRUE)
          max_pval_active <- max(current_data$Pvalues, na.rm = TRUE)
          if(is.finite(min_pval_active) && is.finite(max_pval_active)) {
            update_min_pval <- floor(min_pval_active)
            update_max_pval <- ceiling(max_pval_active)
            if(update_max_pval < update_min_pval) update_max_pval <- update_min_pval # Ensure max >= min
            updateNumericInput(session, "pvalue_min", min = update_min_pval, value = update_min_pval)
            updateNumericInput(session, "pvalue_max", max = update_max_pval + 1, value = update_max_pval) # Add 1 to max for range
          } else {
            updateNumericInput(session, "pvalue_min", value = 0, min = 0)
            updateNumericInput(session, "pvalue_max", value = 30, max = 100) # Reset to default
          }
        } else {
          updateNumericInput(session, "pvalue_min", value = 0, min = 0)
          updateNumericInput(session, "pvalue_max", value = 30, max = 100) # Reset to default
        }
        
        # LD Range (Min/Max = 0/1, reset values)
        if(is.numeric(current_data$LD)){
          updateNumericInput(session, "LD_min", value = min(current_data$LD, 0, na.rm=TRUE)) # Set initial min value based on data
          updateNumericInput(session, "LD_max", value = max(current_data$LD, 1, na.rm=TRUE)) # Set initial max value based on data
        } else {
          updateNumericInput(session, "LD_min", value = 0)
          updateNumericInput(session, "LD_max", value = 1)
        }
        
        
        # Clear the update message after a short delay
        shinyjs::delay(1000, { # Delay 1 second
          output$filter_update_status <- renderUI(NULL)
        })
        
      }, ignoreNULL = TRUE, ignoreInit = TRUE) # Important flags
      
      
      # --- Filter Data for Active Tab ---
      filtered_active_data <- reactive({
        # This reactive now depends on active_tab_data() AND the filter inputs
        current_active_data <- active_tab_data()
        req(current_active_data) # Require data for the active tab
        
        # Apply filtering using the helper function and current input values
        # We pass input directly, load_and_filter_module accesses namespaced inputs
        result <- load_and_filter_module(current_active_data, input)
        
        # Notification logic - now applies to the active tab's data
        if (!is.null(result) && nrow(result) > 0) {
          is_only_tag_after_filter <- all(result$Type == "tag_snp", na.rm = TRUE) || all(is.na(result$final_weight))
          
          if (is_only_tag_after_filter) {
            original_had_candidates <- any(current_active_data$Type != "tag_snp", na.rm = TRUE)
            if(original_had_candidates) {
              # Only show notification if candidates were present initially but filtered out
              # Debounce this? Could be annoying if filtering rapidly.
              # showNotification(
              #    "Applied filters resulted in only the Tag SNP remaining for the current tab.",
              #    type = "warning", duration = 5, id=ns("filter_notification")) # Use ID to prevent stacking
            }
            # Don't show notification if original data only had tag SNP
          }
        } else if (is.null(result) || nrow(result) == 0) {
          # Only show if the original data for the tab was not empty
          if (nrow(current_active_data)>0){
            # showNotification("No data matches the current filter criteria for this tag SNP.",
            # type = "warning", duration = 5, id=ns("filter_notification"))
          }
        }
        
        
        return(result)
      })
      
      
      # --- Render Outputs Dynamically ---
      # --- Render Outputs Dynamically ---
      observe({
        results <- raw_results()
        names_res <- result_names() # Get the potential names (tag SNPs)
        req(results) # Need results to proceed
        
        # Determine number of results/tabs
        num_results <- 0
        is_multiple <- is.list(results) && !is.data.frame(results) && !all(c("plot", "table") %in% names(results))
        is_single_list <- is.list(results) && all(c("plot", "table") %in% names(results))
        
        if(is_multiple) {
          num_results <- length(results)
        } else if (is_single_list){
          num_results <- 1 # Single result structure counts as 1
        }
        
        if (num_results == 0) return() # No valid results
        
        # Use lapply to create outputs for 1 to num_results
        lapply(1:num_results, function(i) {
          local({ # Ensure correct value of 'i' is captured
            current_index <- i
            plot_output_id <- paste0("plot_", current_index)
            table_output_id <- paste0("table_", current_index)
            
            # --- Determine a safe filename base for the current tab ---
            tab_name <- "single_result" # Default
            if (!is.null(names_res) && length(names_res) >= current_index) {
              tab_name <- names_res[[current_index]]
            }
            # Sanitize the name for use in filename (remove spaces, colons, etc.)
            safe_filename_base <- gsub("[^A-Za-z0-9_-]", "_", tab_name)
            dynamic_filename <- paste0("panvar_data_", safe_filename_base)
            # --- End filename determination ---
            
            
            # Render Plot for this index (code unchanged from previous version)
            output[[plot_output_id]] <- renderPlotly({
              req(input$resultTabs == ns(paste0("tab_", current_index)))
              data_to_plot <- filtered_active_data()
              plot_obj <- panvar_plotly_function(
                panvar_results_table = data_to_plot,
                nrows_in_gwas = input$total_rows,
                pvalue_threshold = input$pvalue_threshold
              )
              plot_obj
            }) # End renderPlotly
            
            
            # Render Table for this index with customized export filenames
            output[[table_output_id]] <- renderDT({
              req(input$resultTabs == ns(paste0("tab_", current_index)))
              data_to_display <- filtered_active_data()
              
              if (is.null(data_to_display) || nrow(data_to_display) == 0) {
                return(datatable(data.frame(Message = "No data available based on current filters for this tag SNP."),
                                 options = list(dom = 't', searching = FALSE, paging = FALSE),
                                 rownames = FALSE))
              }
              
              # Standard DT rendering with modified buttons option
              datatable(data_to_display,
                        options = list(
                          pageLength = 10, scrollX = TRUE, scrollY = "350px",
                          deferRender = TRUE, scrollCollapse = TRUE,
                          dom = 'Bfrtip', # Ensure Buttons (B) is included
                          # --- MODIFICATION START: Button configuration ---
                          buttons = list(
                            list(extend = 'copy', text = 'Copy'), # Standard copy
                            list(
                              extend = 'csv', # CSV export
                              filename = dynamic_filename # Use dynamic filename
                            ),
                            list(
                              extend = 'excel', # Excel export
                              filename = dynamic_filename, # Use dynamic filename
                              text = 'Excel' # Optional: Change button text
                            )
                          )
                          # --- MODIFICATION END ---
                        ),
                        filter = 'top', # Enable column filters
                        rownames = FALSE,
                        extensions = c('Buttons'), # Declare extension
                        selection = 'single' # Allow single row selection
              )
              
            }) # End renderDT
          }) # End local()
        }) # End lapply
      }) # End observe block for rendering outputs
      
    } # End function(input, output, session)
  ) # End moduleServer
} # End output_dashboard_Server
# --- END SERVER ---