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

# Rijan: Plotly object generator - Modified for DISCRETE LD Color Scale (Manual Palette)
# panvar_plotly_function_old <- function(panvar_results_table, nrows_in_gwas = NULL, pvalue_threshold = 0.05, point_size = 10, alpha_base = 0.7){
#   
#   # --- 1. Input Validation and Graceful Handling ---
#   if (is.null(panvar_results_table) || nrow(panvar_results_table) == 0) {
#     # Return message plot
#     return(
#       plot_ly() %>%
#         layout(title = "No Data to Display", xaxis = list(visible = FALSE), yaxis = list(visible = FALSE)) %>%
#         add_annotations(text = "No data available based on current filters.", x = 0.5, y = 0.5, xref = "paper", yref = "paper", showarrow = FALSE, font = list(size=14))
#     )
#   }
#   
#   required_plot_cols <- c("BP", "Pvalues", "IMPACT", "LD", "GENE", "Type") # LD is still required for binning
#   missing_plot_cols <- setdiff(required_plot_cols, names(panvar_results_table))
#   if (length(missing_plot_cols) > 0) {
#     warning(paste("Plotting function is missing required columns:", paste(missing_plot_cols, collapse=", ")))
#     # Return error plot
#     return(
#       plot_ly() %>%
#         layout(title = "Plotting Error", xaxis = list(visible = FALSE), yaxis = list(visible = FALSE)) %>%
#         add_annotations(text = paste("Cannot generate plot. Missing columns:", paste(missing_plot_cols, collapse=", ")), x = 0.5, y = 0.5, xref = "paper", yref = "paper", showarrow = FALSE, font = list(size=14))
#     )
#   }
#   
#   # --- 2. Data Preparation & Type Safety ---
#   plot_data <- tryCatch({
#     as.data.frame(panvar_results_table) %>%
#       dplyr::mutate(
#         BP = as.numeric(BP),
#         Pvalues = as.numeric(Pvalues),
#         LD = as.numeric(LD),
#         negLog10P = -log10(Pvalues)
#       ) %>%
#       dplyr::mutate(negLog10P = ifelse(is.infinite(negLog10P), max(negLog10P[is.finite(negLog10P)], 0, na.rm = TRUE) + 1, negLog10P)) %>%
#       # --- *** START: Bin LD into Categories *** ---
#       dplyr::mutate(
#         LD_Category = cut(
#           LD,
#           breaks = c(-Inf, 0.2, 0.4, 0.6, 0.8, 1.0), # Bins: (-Inf, 0.2], (0.2, 0.4], ..., (0.8, 1.0]
#           labels = c("0.0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0"),
#           right = TRUE, # Intervals are closed on the right (e.g., 0.2 falls in first bin)
#           include.lowest = TRUE # Include 0 in the first bin
#         ),
#         # Explicitly handle NAs - convert factor NA to a specific level
#         LD_Category = factor(ifelse(is.na(LD), "NA", as.character(LD_Category)),
#                              levels = c("0.0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0", "NA")) # Ensure correct level order + NA
#       )
#     # --- *** END: Bin LD into Categories *** ---
#   }, error = function(e) {
#     warning("Error during data preparation for plotting: ", e$message)
#     return(NULL)
#   })
#   
#   if (is.null(plot_data) || !is.numeric(plot_data$BP) || !is.numeric(plot_data$negLog10P)) {
#     warning("One or more required columns (BP, Pvalues) could not be coerced to numeric.")
#     # Return error plot
#     return(
#       plot_ly() %>%
#         layout(title = "Plotting Error", xaxis = list(visible = FALSE), yaxis = list(visible = FALSE)) %>%
#         add_annotations(text = "Cannot generate plot. BP or Pvalues column is not numeric.", x = 0.5, y = 0.5, xref = "paper", yref = "paper", showarrow = FALSE, font = list(size=14))
#     )
#   }
#   
#   
#   # --- 3. Plot Specific Calculations ---
#   tag_df <- plot_data %>% dplyr::filter(Type == 'tag_snp')
#   tag_snp_bp <- if(nrow(tag_df) > 0) tag_df %>% dplyr::pull(BP) %>% unique() %>% head(1) else NULL
#   
#   if(is.null(nrows_in_gwas)) nrows_in_gwas <- 1e6
#   if(is.null(pvalue_threshold)) pvalue_threshold <- 0.05
#   bonf.cor <- -log10(pvalue_threshold / nrows_in_gwas)
#   
#   # --- 4. Define Discrete Colors ---
#   # Define 5 colors for the bins + 1 color for NA
#   
#   # --- !!! MANUALLY DEFINED COLORS !!! ---
#   # Replace these 5 hex codes with your desired palette, ordered from low LD to high LD bin
#   ld_bin_colors <- c("#005f73", "#5e548e", "#9d0208", "#ae4a24", "#606c38") # Dark Teal, Dark Purple, Dark Red, Dark Orange/Brown, Dark Olive
#   # --- !!! END MANUAL DEFINITION !!! ---
#   
#   na_color <- "#808080" # Grey for NA (can also be changed)
#   # Combine colors in the order of the factor levels defined above
#   discrete_colors <- c(ld_bin_colors, na_color)
#   
#   
#   # --- 5. Build Plotly Object using Discrete Color Mapping ---
#   p <- plot_ly(
#     data = plot_data,
#     x = ~BP,
#     y = ~Pvalues,
#     # --- *** Map color to the CATEGORICAL variable *** ---
#     color = ~LD_Category,
#     # --- *** Provide the DISCRETE color palette *** ---
#     colors = discrete_colors, # Use the manually defined colors
#     type = 'scatter',
#     mode = 'markers',
#     symbol = ~IMPACT,     # Symbol still mapped to IMPACT
#     marker = list(
#       size = point_size,
#       opacity = alpha_base
#       # NO colorbar needed here
#     ),
#     hoverinfo = 'text',
#     # --- Update hover text to show LD category ---
#     text = ~paste(
#       'BP:', BP,
#       '<br>-log10(P):', round(Pvalues, 2),
#       '<br>LD Cat.:', LD_Category, # Show category
#       '<br>LD Val:', round(LD, 2),   # Optionally show original LD value too
#       '<br>Impact:', IMPACT,
#       '<br>Gene:', GENE
#     )
#     # Legend definition handled in layout
#   )
#   
#   # --- 6. Define and Apply Layout (Adjusted for Discrete Legends) ---
#   hline <- function(y = 0, color = "grey") list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = y, y1 = y, yref = "y", line = list(color = color, dash = "dash"))
#   vline <- function(x = 0, color = "red") list(type = "line", y0 = 0, y1 = 1, yref = "paper", x0 = x, x1 = x, xref = "x", line = list(color = color, dash = "dot", width=1.5))
#   
#   shapes_list <- list(hline(bonf.cor))
#   if (!is.null(tag_snp_bp) && is.finite(tag_snp_bp)) {
#     shapes_list <- c(shapes_list, list(vline(tag_snp_bp)))
#   }
#   
#   # Define the consolidated layout list
#   final_layout <- list(
#     title = "Interactive PanvaR Plot",
#     xaxis = list(title = "Position (BP)", titlefont = list(size = 14)),
#     yaxis = list(title = "-log<sub>10</sub>(P-value)", titlefont = list(size = 14)),
#     shapes = shapes_list,
#     # --- *** Define Legend(s) *** ---
#     # Plotly should automatically create legends for 'color' and 'symbol'
#     # We primarily control their appearance and position here.
#     legend = list(
#       title = list(text = '<b>Legend</b>'), # Generic title or specific if only one shows clearly
#       x = 1.02, # Position horizontally (right of plot)
#       y = 1.0,  # Position vertically (top of plot)
#       xanchor = 'left', # Anchor legend from its left
#       yanchor = 'top',  # Anchor legend from its top
#       traceorder = 'grouped', # Group items by trace (color first, then symbol)
#       tracegroupgap = 10, # Gap between the color group and symbol group
#       bgcolor = 'rgba(255,255,255,0.7)', # Semi-transparent background
#       borderwidth = 1,
#       bordercolor = '#CCCCCC'
#     ),
#     # --- *** REMOVE coloraxis *** ---
#     # coloraxis = NULL, # Not needed for discrete colors mapped via 'colors'
#     margin = list(r = 180) # Keep or adjust right margin as needed
#   )
#   
#   
#   # Apply the final layout
#   p <- layout(p,
#               title = final_layout$title,
#               xaxis = final_layout$xaxis,
#               yaxis = final_layout$yaxis,
#               shapes = final_layout$shapes,
#               legend = final_layout$legend,
#               margin = final_layout$margin
#               # No coloraxis argument here
#   )
#   
#   # --- 7. Return Plotly Object ---
#   return(p)
#   
# }

panvar_plotly_function <- function(panvar_results_table, nrows_in_gwas = NULL, pvalue_threshold = 0.05, point_size = 10, alpha_base = 0.7, window.size = 5e5){
  p <-
    panvar_plot(panvar_results_table, 
                nrows_in_gwas = nrows_in_gwas,
                pvalue_threshold = pvalue_threshold,
                point_size = point_size,
                alpha_base = alpha_base,
                window.size = window.size)
  
  gp <- 
    ggplotly(p)
  
  return(gp)
  
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
        
        # --- START: Add Bulk Download Button ---
        hr(), # Add a horizontal rule for separation
        h4("Downloads"),
        downloadButton(ns("bulk_download"), "Download All Results (.zip)",
                       class = "btn-primary", # Optional: Style the button
                       style="width:100%; margin-bottom: 15px;"), # Optional: Style
        bsTooltip(
          id = ns("bulk_download"),
          title = "Download all generated plots (as PNG) and tables (as CSV) in a single zip file.",
          placement = "right",
          trigger = "hover"
        ),
        hr(), # Add another separator
        h4("Filters (Apply to Selected Tab)"),
        # --- END: Add Bulk Download Button ---
        
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

# Output_dashboard_panvar.r
# --- START SERVER ---
output_dashboard_Server <- function(id, shared) {
  moduleServer(
    id,
    function(input, output, session) {
      ns <- session$ns # Ensure ns is available
      
      # Ensure ggplot2 is available for saving plots
      # require(ggplot2) # Uncomment if not loaded globally
      # require(zip) # Uncomment if not loaded globally
      
      
      # Reactive value to store the raw results (list or single)
      raw_results <- reactiveVal(NULL)
      # Reactive value to store the names of the results (tag SNPs) if it's a list
      result_names <- reactiveVal(NULL)
      # Reactive value to store the original ggplot objects (if available)
      # We need this because panvar_plotly_function creates plotly, but we want to save ggplot
      original_plots <- reactiveVal(NULL)
      
      
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
        # Use tryCatch to handle potential errors if input$Pre_existing_panvaR_results_path is invalid
        file_info <- tryCatch({
          shinyFiles::parseFilePaths(c(Home = fs::path_home()), input$Pre_existing_panvaR_results_path)
        }, error = function(e) NULL)
        
        if (!is.null(file_info) && nrow(file_info) > 0) {
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
        plots_to_store <- list() # Initialize list for original ggplots
        
        if (source_type == "dynamic") {
          # Use results directly from the shared reactive value
          if (!is.null(shared$analysis_results)) {
            results_to_store <- shared$analysis_results
            # Try to get names if it's a list (might need adjustment based on how names are set in panvar_func)
            # Check if it's a list BUT not a data.frame (data.frames are lists) and not the single result structure
            if (is.list(results_to_store) && !is.data.frame(results_to_store) && !all(c("plot", "table") %in% names(results_to_store))) {
              names_to_store <- names(results_to_store)
              # If names are NULL, generate default names based on index or try to extract from table
              if (is.null(names_to_store) || any(names_to_store == "")) {
                possible_names <- sapply(results_to_store, function(res) {
                  if (is.list(res) && "table" %in% names(res) && is.data.frame(res$table)) {
                    tag_row <- res$table %>% filter(Type == "tag_snp") %>% head(1)
                    if (nrow(tag_row) > 0) return(paste0(tag_row$CHROM, ":", tag_row$BP))
                  }
                  return(NA_character_)
                })
                default_names <- paste("Result", seq_along(results_to_store))
                names_to_store <- ifelse(is.na(possible_names), default_names, possible_names)
                
              }
              # Extract original plots
              plots_to_store <- lapply(results_to_store, function(res) if(is.list(res) && "plot" %in% names(res)) res$plot else NULL)
              
            } else if (is.list(results_to_store) && all(c("plot", "table") %in% names(results_to_store))) {
              # Handle the case where dynamic result is a single item but not in a list of length 1
              # We can try to extract the tag SNP name from the table
              single_name <- "Single Result"
              if (is.data.frame(results_to_store$table)) {
                tag_row <- results_to_store$table %>% filter(Type == "tag_snp") %>% head(1)
                if (nrow(tag_row) > 0) {
                  single_name <- paste0(tag_row$CHROM, ":", tag_row$BP)
                }
              }
              names_to_store <- single_name
              plots_to_store <- list(results_to_store$plot) # Store the single plot in a list
              
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
                is_list_of_results <- is.list(loaded_object) && !is.data.frame(loaded_object) && length(loaded_object)>0 && is.list(loaded_object[[1]]) && all(c("plot", "table") %in% names(loaded_object[[1]]))
                is_single_result <- is.list(loaded_object) && all(c("plot", "table") %in% names(loaded_object))
                
                if (is_list_of_results) {
                  # Assume it's a list of results
                  results_to_store <- loaded_object
                  names_from_rds <- names(results_to_store)
                  if (is.null(names_from_rds) || any(names_from_rds == "")) {
                    # Try to extract names from tables
                    possible_names <- sapply(results_to_store, function(res) {
                      if (is.list(res) && "table" %in% names(res) && is.data.frame(res$table)) {
                        tag_row <- res$table %>% filter(Type == "tag_snp") %>% head(1)
                        if (nrow(tag_row) > 0) return(paste0(tag_row$CHROM, ":", tag_row$BP))
                      }
                      return(NA_character_)
                    })
                    default_names <- paste("Result", seq_along(results_to_store))
                    names_to_store <- ifelse(is.na(possible_names), default_names, possible_names)
                  } else {
                    names_to_store <- names_from_rds
                  }
                  # Extract plots
                  plots_to_store <- lapply(results_to_store, function(res) if(is.list(res) && "plot" %in% names(res)) res$plot else NULL)
                  
                } else if (is_single_result) {
                  # Assume it's a single result structure
                  results_to_store <- loaded_object
                  single_name <- "Single Result"
                  # Try to get name for the single result
                  if (is.data.frame(results_to_store$table)) {
                    tag_row <- results_to_store$table %>% filter(Type == "tag_snp") %>% head(1)
                    if (nrow(tag_row) > 0) single_name <- paste0(tag_row$CHROM, ":", tag_row$BP)
                  }
                  names_to_store <- single_name
                  plots_to_store <- list(results_to_store$plot) # Store single plot
                  
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
                single_name <- "Single Result File"
                # Try to get name for the single result
                tag_row <- single_table %>% filter(Type == "tag_snp") %>% head(1)
                if (nrow(tag_row) > 0) single_name <- paste0(tag_row$CHROM, ":", tag_row$BP)
                names_to_store <- single_name
                plots_to_store <- list(NULL) # No plot available from table file
                showNotification("Loaded table from file. Plot generation/download from file is not supported.", type = "info", duration=7)
                
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
                plots_to_store <- NULL
              }
              
            }, error = function(e) {
              showNotification(paste("Error loading or processing file:", e$message), type = "error")
              results_to_store <- NULL
              names_to_store <- NULL
              plots_to_store <- NULL
            })
          } else {
            # This case should ideally not be reached if req() is used, but added for safety
            results_to_store <- NULL
            names_to_store <- NULL
            plots_to_store <- NULL
          }
        }
        raw_results(results_to_store)
        result_names(names_to_store)
        original_plots(plots_to_store) # Store the original plots
        
        # Enable/disable filter controls and download buttons based on whether data loaded
        if (is.null(results_to_store)) {
          shinyjs::disable("filter_controls")
          shinyjs::disable("bulk_download") # Disable bulk download if no results
          output$filter_update_status <- renderUI({ # Clear status message
            p(em("Load data or run analysis to enable filters and downloads."), style="color:gray;")
          })
        } else {
          shinyjs::enable("filter_controls")
          shinyjs::enable("bulk_download") # Enable bulk download
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
        is_multiple_results <- is.list(results) && !is.data.frame(results) &&
          !all(c("plot", "table") %in% names(results)) && length(results) > 0
        
        # Check if it's a single result structure
        is_single_result <- is.list(results) && all(c("plot", "table") %in% names(results))
        
        if (is_multiple_results) {
          # Ensure names are available and match length
          if(is.null(names_res) || length(names_res) != length(results)){
            names_res <- paste("Result", seq_along(results)) # Fallback names
            warning("Result names missing or incorrect length for tabs, using default names.")
          }
          # Create tabs
          tabs <- lapply(seq_along(results), function(i) {
            tab_title <- names_res[[i]]
            tab_value <- ns(paste0("tab_", i))
            plot_output_id <- ns(paste0("plot_", i))
            table_output_id <- ns(paste0("table_", i))
            plot_download_button_id <- ns(paste0("download_plot_", i)) # ID for the download button
            
            tabPanel(title = tab_title,
                     value = tab_value,
                     plotlyOutput(plot_output_id),
                     # --- Add Plot Download Button UI ---
                     downloadButton(plot_download_button_id, "Download Plot (.png)",
                                    class="btn-sm", style="margin-top: 5px; margin-bottom: 15px;"),
                     bsTooltip(
                       id = plot_download_button_id,
                       title = "Download the original plot as a PNG image. Filters applied in the UI will not affect the downloaded plot.",
                       placement = "right",
                       trigger = "hover"
                     ),
                     # --- End Plot Download Button UI ---
                     DTOutput(table_output_id)
            )
          })
          do.call(tabsetPanel, c(tabs, id = ns("resultTabs"))) # Use ns() for tabsetPanel ID
          
        } else if (is_single_result) {
          # Handle single result (could be direct structure or list of length 1 wrapped by mistake)
          single_result_name <- if(!is.null(names_res) && length(names_res)==1) names_res[[1]] else "Single Result"
          
          single_plot_output_id <- ns("plot_1")
          single_table_output_id <- ns("table_1")
          single_plot_download_id <- ns("download_plot_1") # Use index 1
          
          do.call(tabsetPanel, list(
            id = ns("resultTabs"), # Still give it an ID
            tabPanel(title = single_result_name,
                     value = ns("tab_1"), # Use index 1 for consistency
                     plotlyOutput(single_plot_output_id), # Use indexed IDs
                     # Add download button for single plot
                     downloadButton(single_plot_download_id, "Download Plot (.png)",
                                    class="btn-sm", style="margin-top: 5px; margin-bottom: 15px;"),
                     bsTooltip(
                       id = single_plot_download_id,
                       title = "Download the original plot as a PNG image. Filters applied in the UI will not affect the downloaded plot.",
                       placement = "right",
                       trigger = "hover"
                     ),
                     DTOutput(single_table_output_id)       # Use indexed IDs
            )
          ))
          
        } else {
          # Handle invalid structure or empty single result
          div(style = "padding: 20px; text-align: center;", h4("Loaded result has an unexpected format or is empty."))
        }
      })
      
      # --- Reactive for Currently Selected Tab Data Table ---
      active_tab_data <- reactive({
        results <- raw_results()
        selected_tab_value <- input$resultTabs # Value like ns("tab_1")
        
        # Directly require results, but handle NULL input$resultTabs later
        req(results)
        
        selected_tab_index <- 1 # Default to the first tab initially
        
        # If input$resultTabs is not NULL, parse the index from its value
        if (!is.null(selected_tab_value)) {
          # Extract index from tab value (e.g., "module2-tab_1" -> 1)
          match_val <- regmatches(selected_tab_value, regexpr("\\d+$", selected_tab_value))
          if(length(match_val) > 0) {
            selected_tab_index <- as.integer(match_val)
          } else {
            # Fallback if parsing fails, though unlikely with ns()
            warning("Could not parse index from tab value: ", selected_tab_value)
            return(NULL)
          }
        }
        # Now selected_tab_index is either 1 (default) or parsed value
        
        # Determine if results represent multiple entries or a single one
        is_multiple_results <- is.list(results) &&
          !is.data.frame(results) &&
          !all(c("plot", "table") %in% names(results)) &&
          length(results) > 0
        
        is_single_result <- is.list(results) && all(c("plot", "table") %in% names(results))
        
        current_item <- NULL
        if (is_multiple_results) {
          # Select the item from the list based on index
          if (selected_tab_index > 0 && selected_tab_index <= length(results)) {
            current_item <- results[[selected_tab_index]]
          } else {
            warning(paste("Selected tab index", selected_tab_index, "out of bounds for results list."))
            return(NULL) # Index out of bounds
          }
        } else if (is_single_result) {
          # It's potentially a single result (either direct or list of 1)
          # We still use index 1 because the UI renders it as ns("tab_1")
          if(selected_tab_index != 1) {
            # This shouldn't happen if it's truly a single result UI
            warning("Unexpected tab index for single result view.")
            return(NULL)
          }
          current_item <- results # The results object itself is the single item
        }
        
        # Validate the structure and return the table
        if(is.list(current_item) && "table" %in% names(current_item) && (is.data.frame(current_item$table) || is.data.table(current_item$table)) ){
          # Ensure it's a data.table for consistency downstream if needed
          return(as.data.table(current_item$table))
        } else {
          # This might be normal if a result item is invalid
          # warning(paste("Result item for selected tab index", selected_tab_index, "has invalid structure or no table."))
          return(NULL) # Return NULL if structure is wrong or table missing
        }
      }) # End active_tab_data reactive
      
      # --- Reactive for Currently Selected Tab's ORIGINAL Plot Object ---
      # This retrieves the ggplot object, not the plotly one
      active_tab_original_plot <- reactive({
        plots <- original_plots() # Get the list of original plot objects
        selected_tab_value <- input$resultTabs
        req(plots, selected_tab_value) # Require plots list and selected tab
        
        match_val <- regmatches(selected_tab_value, regexpr("\\d+$", selected_tab_value))
        if(length(match_val) == 0) return(NULL)
        selected_tab_index <- as.integer(match_val)
        
        # Check if index is valid
        if (selected_tab_index > 0 && selected_tab_index <= length(plots)) {
          plot_obj <- plots[[selected_tab_index]]
          # Check if it's a ggplot object
          if (ggplot2::is.ggplot(plot_obj)) {
            return(plot_obj)
          } else {
            # warning("Original plot object for this tab is not a ggplot object.")
            return(NULL)
          }
        } else {
          warning("Selected tab index out of bounds for original plots list.")
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
            # Use unique() directly on the vector, omit NAs first
            valid_choices <- tryCatch(
              sort(unique(na.omit(choices))),
              error = function(e) { # Handle potential errors during sort (e.g., mixed types)
                warning(paste("Error sorting choices for", inputId, ":", e$message))
                unique(na.omit(choices)) # Fallback to unsorted unique if sort fails
              }
            )
          }
          # Keep existing selection if it's still valid, otherwise clear
          selected_val <- intersect(current_selection, valid_choices)
          if (length(selected_val) == 0) selected_val <- character(0)
          
          updateSelectizeInput(session, inputId, choices = valid_choices, selected = selected_val, server = TRUE) # server=TRUE can improve performance for large lists
        }
        
        # Use existing input values as current_selection where appropriate
        # Add error handling for potential missing columns during filtering
        tryCatch({ safe_update_selectize("selected_genes", current_data$GENE, input$selected_genes) }, error=function(e) warning(paste("Column GENE not found:", e$message)))
        tryCatch({ safe_update_selectize("selected_effect_types", current_data$EFFECT, input$selected_effect_types) }, error=function(e) warning(paste("Column EFFECT not found:", e$message)))
        tryCatch({ safe_update_selectize("selected_amino_acid", current_data$AA, input$selected_amino_acid) }, error=function(e) warning(paste("Column AA not found:", e$message)))
        tryCatch({ safe_update_selectize("selected_REF_types", current_data$REF, input$selected_REF_types) }, error=function(e) warning(paste("Column REF not found:", e$message)))
        tryCatch({ safe_update_selectize("selected_ALT_types", current_data$ALT, input$selected_ALT_types) }, error=function(e) warning(paste("Column ALT not found:", e$message)))
        
        
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
          min_ld_active <- min(current_data$LD, 0, na.rm = TRUE) # Handle potential all NAs -> Inf
          max_ld_active <- max(current_data$LD, 1, na.rm = TRUE) # Handle potential all NAs -> -Inf
          
          updateNumericInput(session, "LD_min", value = if(is.finite(min_ld_active)) min_ld_active else 0)
          updateNumericInput(session, "LD_max", value = if(is.finite(max_ld_active)) max_ld_active else 1)
          
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
          # Check if only tag SNP(s) remain after filtering
          # Handle cases where 'final_weight' might not exist or be all NA
          is_only_tag_after_filter <- all(result$Type == "tag_snp", na.rm = TRUE) ||
            (!("final_weight" %in% names(result)) || all(is.na(result$final_weight)))
          
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
      observe({
        results <- raw_results()
        names_res <- result_names() # Get the potential names (tag SNPs)
        req(results) # Need results to proceed
        
        # Determine number of results/tabs
        num_results <- 0
        is_multiple <- is.list(results) && !is.data.frame(results) && !all(c("plot", "table") %in% names(results)) && length(results) > 0
        is_single <- is.list(results) && all(c("plot", "table") %in% names(results))
        
        if(is_multiple) {
          num_results <- length(results)
        } else if (is_single){
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
            tab_name <- "result" # Default
            if (!is.null(names_res) && length(names_res) >= current_index) {
              # Attempt to use provided name, sanitize it
              safe_tab_name <- gsub("[^A-Za-z0-9_.-]", "_", names_res[[current_index]]) # Allow dot and hyphen
              # Avoid names starting with non-alphanumeric or being empty after sanitization
              if (nchar(safe_tab_name) > 0 && grepl("^[A-Za-z0-9]", safe_tab_name)) {
                tab_name <- safe_tab_name
              } else {
                tab_name <- paste0("result_", current_index) # Fallback if sanitization fails
              }
            } else {
              tab_name <- paste0("result_", current_index) # Fallback if names_res is invalid
            }
            # --- End filename determination ---
            
            # Render Plotly Plot for this index
            output[[plot_output_id]] <- renderPlotly({
              # Important: Ensure this reactive dependency only triggers re-render
              # when the specific tab is active.
              req(input$resultTabs == ns(paste0("tab_", current_index)))
              
              data_to_plot <- filtered_active_data() # Use the filtered data for the active tab
              
              # Generate the plotly plot using the function
              plot_obj <- panvar_plotly_function(
                panvar_results_table = data_to_plot,
                nrows_in_gwas = input$total_rows,
                pvalue_threshold = input$pvalue_threshold)
              plot_obj # Return the plotly object
            }) # End renderPlotly
            
            
            # Render Table for this index with customized export filenames
            output[[table_output_id]] <- renderDT({
              # Ensure dependency on active tab
              req(input$resultTabs == ns(paste0("tab_", current_index)))
              data_to_display <- filtered_active_data() # Use filtered data
              
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
                              filename = paste0("panvar_table_", tab_name) # Use dynamic filename
                            ),
                            list(
                              extend = 'excel', # Excel export
                              filename = paste0("panvar_table_", tab_name), # Use dynamic filename
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
            
            # --- START: Individual Plot Download Handler ---
            plot_download_id <- paste0("download_plot_", current_index)
            output[[plot_download_id]] <- downloadHandler(
              filename = function() {
                paste0("panvar_plot_", tab_name, ".png") # Use sanitized name
              },
              content = function(file) {
                # Get the original plot for the current index
                plots <- original_plots()
                plot_obj <- NULL
                if (current_index > 0 && current_index <= length(plots)) {
                  plot_obj <- plots[[current_index]]
                }
                
                # Check if it's a ggplot object and save
                if (ggplot2::is.ggplot(plot_obj)) {
                  # Use ggsave; specify device if needed, adjust dpi, width, height
                  ggsave(file, plot = plot_obj, device = "png", width = 10, height = 6, dpi = 300)
                } else {
                  # Handle cases where plot is NULL or not ggplot
                  # Option 1: Show notification
                  showNotification(paste("Original plot for", tab_name, "is not available or not a ggplot object. Cannot download."), type = "warning")
                  # Option 2: Create a placeholder image (e.g., using base R plot)
                  png(file, width=400, height=200)
                  plot.new()
                  text(0.5, 0.5, "Plot not available for download", cex = 1.2)
                  dev.off()
                  # Option 3: Return error (download will fail client-side)
                  # stop("Plot not available")
                }
              }
            )
            # --- END: Individual Plot Download Handler ---
            
          }) # End local()
        }) # End lapply
      }) # End observe block for rendering outputs
      
      
      # --- START: Bulk Download Handler ---
      output$bulk_download <- downloadHandler(
        filename = function() {
          paste0("panvar_all_results_", Sys.Date(), ".zip")
        },
        content = function(file) {
          # Get all raw results and names
          all_results <- raw_results()
          all_names <- result_names()
          all_plots <- original_plots() # Get original ggplot objects
          
          if (is.null(all_results) || length(all_results) == 0) {
            showNotification("No results available to download.", type = "warning")
            # Create an empty zip file or handle error appropriately
            file.create(file) # Create an empty file to prevent Shiny error
            return()
          }
          
          # Create a temporary directory
          temp_dir <- tempdir()
          download_dir <- file.path(temp_dir, paste0("panvar_bulk_", format(Sys.time(), "%Y%m%d_%H%M%S")))
          dir.create(download_dir, recursive = TRUE)
          # Ensure cleanup of the temp directory
          on.exit(unlink(download_dir, recursive = TRUE), add = TRUE)
          
          # Determine if it's a list of results or a single result structure
          is_multiple <- is.list(all_results) && !is.data.frame(all_results) && !all(c("plot", "table") %in% names(all_results))
          is_single <- is.list(all_results) && all(c("plot", "table") %in% names(all_results))
          
          num_items <- 0
          if(is_multiple) {
            num_items <- length(all_results)
          } else if (is_single) {
            num_items <- 1
          }
          
          if (num_items == 0) { # Should not happen if check above passed, but safety first
            file.create(file)
            return()
          }
          
          
          # --- Loop through results and save ---
          files_to_zip <- c() # Keep track of files created
          
          withProgress(message = 'Preparing bulk download...', value = 0, {
            for (i in 1:num_items) {
              incProgress(1/num_items, detail = paste("Processing result", i))
              
              # Get current item, name, and plot
              current_item <- if(is_multiple) all_results[[i]] else all_results
              current_name <- if(!is.null(all_names) && length(all_names) >= i) all_names[[i]] else paste0("result_", i)
              current_plot <- if(!is.null(all_plots) && length(all_plots) >= i) all_plots[[i]] else NULL
              
              # Sanitize name for filename
              safe_name <- gsub("[^A-Za-z0-9_.-]", "_", current_name)
              if (nchar(safe_name) == 0 || !grepl("^[A-Za-z0-9]", safe_name)) {
                safe_name <- paste0("result_", i)
              }
              
              # --- Save Table ---
              if (is.list(current_item) && "table" %in% names(current_item) && (is.data.frame(current_item$table) || is.data.table(current_item$table))) {
                table_data <- as.data.table(current_item$table)
                if (nrow(table_data) > 0) {
                  csv_filename <- file.path(download_dir, paste0("panvar_table_", safe_name, ".csv"))
                  tryCatch({
                    fwrite(table_data, csv_filename)
                    files_to_zip <- c(files_to_zip, csv_filename)
                  }, error = function(e) {
                    warning(paste("Failed to write table for", safe_name, ":", e$message))
                  })
                }
              }
              
              # --- Save Plot (Original ggplot) ---
              if (ggplot2::is.ggplot(current_plot)) {
                plot_filename <- file.path(download_dir, paste0("panvar_plot_", safe_name, ".png"))
                tryCatch({
                  ggsave(plot_filename, plot = current_plot, device = "png", width = 10, height = 6, dpi = 300)
                  files_to_zip <- c(files_to_zip, plot_filename)
                }, error = function(e) {
                  warning(paste("Failed to save plot for", safe_name, ":", e$message))
                })
              } else {
                # Optionally create a placeholder or log that plot wasn't available/ggplot
                # print(paste("Plot for", safe_name, "is not a valid ggplot object - skipping plot save."))
              }
            } # End loop
          }) # End withProgress
          
          # --- Create Zip File ---
          if (length(files_to_zip) > 0) {
            # Use the zip library
            # Need to specify the root directory for zipping to avoid including full temp path
            # zip::zip(zipfile = file, files = files_to_zip, root = download_dir) # Old syntax?
            setwd(download_dir) # Temporarily change working directory
            zip::zip(zipfile = file, files = list.files(download_dir)) # Zip contents relative to download_dir
            setwd(temp_dir) # Change back
          } else {
            showNotification("No valid files were generated to include in the zip archive.", type = "warning")
            file.create(file) # Create empty file
          }
        },
        contentType = "application/zip"
      )
      # --- END: Bulk Download Handler ---
      
      
    } # End function(input, output, session)
  ) # End moduleServer
} # End output_dashboard_Server
# --- END SERVER ---