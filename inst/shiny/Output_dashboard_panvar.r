# Output_dashboard_panvar.r

# ---
# Helper functions
# Enhanced load and filter module with additional filtering options
# Enhanced load and filter module with NULL handling
load_and_filter_module <- function(my_data, input) {
  # Check if data is NULL or empty
  if (is.null(my_data) || nrow(my_data) == 0) {
    return(NULL)
  }
  
  # Safely get tag_snp with error handling
  tag_snp <- tryCatch({
    my_data %>%
      filter(Type == "tag_snp") %>%
      pull(BP) %>%
      unique() %>%
      as.numeric()
  }, error = function(e) {
    warning("Could not extract tag_snp: ", e$message)
    return(NULL)
  })
  
  # If no tag_snp found, return NULL (or handle differently if needed)
  if (is.null(tag_snp)) {
    # This case might need specific handling depending on requirements
    # If tag_snp is essential for filtering, returning NULL might be appropriate
    # Otherwise, filtering might proceed without BP range filter
    warning("No tag_snp found in data, BP range filtering disabled.")
    # return(NULL) # Option: Return NULL if tag_snp is essential
  }
  
  # Safely apply filters with error handling
  filtered_data <- tryCatch({
    result <- my_data
    
    # Apply numeric filters if input exists
    if (!is.null(input$LD_min) && !is.null(input$LD_max)) {
      result <- result %>%
        filter(LD >= as.numeric(input$LD_min) & LD <= as.numeric(input$LD_max))
    }
    
    # Apply BP filter only if tag_snp was successfully found
    if (!is.null(tag_snp) && !is.null(input$BP_LHS) && !is.null(input$BP_RHS)) {
      result <- result %>%
        filter(BP >= as.numeric(tag_snp - input$BP_LHS) &
                 BP <= as.numeric(tag_snp + input$BP_RHS))
    }
    
    
    if (!is.null(input$pvalue_min) && !is.null(input$pvalue_max)) {
      result <- result %>%
        filter(Pvalues >= as.numeric(input$pvalue_min) &
                 Pvalues <= as.numeric(input$pvalue_max))
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


# Rijan: Ploty object generator
panvar_plotly_function <- function(panvar_results_table, nrows_in_gwas = NULL, pvalue_threshold = 0.05, point_size = 3, alpha_base = 0.7){
  
  # Rijan: What is the tag_snp
  tag_df <- panvar_results_table %>%
    filter(Type == 'tag_snp')
  
  # Handle cases where tag_df might be empty after filtering
  if(nrow(tag_df) > 0) {
    tag_snp <- tag_df %>%
      pull(BP) %>%
      unique() %>%
      head(1) # Ensure only one value if multiple rows match
  } else {
    tag_snp <- NULL # Set to NULL if no tag SNP found
  }
  
  
  # Rijan: Did the user supply a default value for the bonforoni correction?
  if(is.null(nrows_in_gwas)){
    print("You did not supply a value for the number of tests that were in the GWAS - a place holder value will be used. This is not ideal.")
    nrows_in_gwas <- 5e6 # Use a default if not provided
  }
  
  # Rijan: What is the Bonferroni correction
  bonf.cor <- -log10(pvalue_threshold / nrows_in_gwas)
  
  # Rijan: The main Object that we will gradually add to.
  panvar_plotly <- plot_ly(
    panvar_results_table,
    x = ~BP,
    y = ~Pvalues,
    type = 'scatter',
    mode = 'markers',
    symbol = ~IMPACT,
    marker = list(
      size = 10 # Adjusted size for better visibility
    ),
    # Add hover text for more info
    hoverinfo = 'text',
    text = ~paste('BP:', BP, '<br>P-value:', round(Pvalues, 2), '<br>LD:', round(LD, 2), '<br>Impact:', IMPACT, '<br>Gene:', GENE),
    showlegend = TRUE,
    name = ~IMPACT # Group points by IMPACT in legend
  )
  
  # Add layout details
  panvar_plotly <- layout(panvar_plotly, title = "Interactive PanvaR Plot",
                          xaxis = list(title = "Position (BP)", titlefont = list(size = 14)),
                          yaxis = list(title = "-log<sub>10</sub>(P-value)", titlefont = list(size = 14)), # Use subscript
                          legend = list(title = list(text = '<b>Impact</b>'), orientation = "h", y = -0.2) # Bold title, horizontal legend below
  )
  
  
  # --- Functions for lines ---
  hline <- function(y = 0, color = "blue") {
    list(
      type = "line",
      x0 = 0, x1 = 1, xref = "paper", # Stretches across plot width
      y0 = y, y1 = y, yref = "y",    # Position based on y-axis value
      line = list(color = color, dash = "dash"),
      name = "Bonferroni Threshold" # Name for potential hover/legend
    )
  }
  
  vline <- function(x = 0, color = "red") {
    list(
      type = "line",
      y0 = 0, y1 = 1, yref = "paper", # Stretches across plot height
      x0 = x, x1 = x, xref = "x",    # Position based on x-axis value
      line = list(color = color, dash = "dot"),
      name = "Tag SNP" # Name for potential hover/legend
    )
  }
  # --- End Functions for lines ---
  
  # List to hold shapes
  shapes_list <- list(hline(bonf.cor)) # Start with Bonferroni line
  
  # Add vertical line for tag SNP only if tag_snp is not NULL
  if (!is.null(tag_snp)) {
    shapes_list <- c(shapes_list, list(vline(tag_snp)))
  }
  
  # Add shapes to the layout
  panvar_plotly <- layout(panvar_plotly, shapes = shapes_list)
  
  return(panvar_plotly)
  
}
# ---

# ---
# The UI logic (remains unchanged from your provided code)

output_dashboard_UI <- function(id) {
  ns <- NS(id)
  tagList(
    tabPanel(
      "PanvaR output exploration Dashboard"
    ),
    sidebarLayout(
      sidebarPanel(
        width = 3,
        
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
              multiple = FALSE
            ),
            tags$span(
              id = ns("Pre_existing_panvaR_results_path_tooltip"),
              icon("question-circle"),
              style = "color: green;"
            )
          ),
          bsTooltip(
            id = ns("Pre_existing_panvaR_results_path_tooltip"),
            title = "Please provide the path to a pre-existing panvaR results table.",
            placement = "right",
            trigger = "hover"
          ),
          textOutput(ns("Pre_existing_panvaR_results"))
        ),
        
        # Rest of the filtering controls - disabled when no data is loaded
        conditionalPanel(
          # Simplified condition: Enable if dynamic OR if file is chosen
          condition = sprintf("input['%s'] === 'dynamic' || input['%s'] === 'file'",
                              ns("data_source"), ns("data_source")),
          
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
            title = "What genes should the results table be filtered for?",
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
            title = "What effect types should the results be filtered for?",
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
            title = "What amino acids should the results be filtered for?",
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
            title = "What ALT types should the results be filtered for?",
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
            title = "What REF types should the results be filtered for?",
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
            title = "What LD range should the results be filtered for?",
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
            title = "Filter results by BP range relative to the tag SNP. LHS = Max distance left, RHS = Max distance right.",
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
            title = "Filter by -log10(P-value) range.",
            placement = "right",
            trigger = "hover"
          ),
          
          
          # Bonferroni correction inputs (moved to plot function, remove from UI if desired)
          # If keeping for dynamic plot adjustment:
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
            title = "Parameters for calculating the Bonferroni correction line on the plot (alpha / total tests).",
            placement = "right",
            trigger = "hover"
          )
        ) # End conditionalPanel
      ), # End sidebarPanel
      
      mainPanel(
        plotlyOutput(ns("data_plotly")),
        DTOutput(ns("dataTable"))
      )
    ) # End sidebarLayout
  ) # End tagList
}
# --- END UI ---


# --- START SERVER ---
output_dashboard_Server <- function(id, shared) {
  moduleServer(
    id,
    function(input, output, session) {
      
      # File chooser setup
      shinyFileChoose(
        input,
        "Pre_existing_panvaR_results_path",
        roots = c(Home = fs::path_home()),
        session = session
      )
      
      # Reactive value to store the data (either dynamic or from file)
      my_data <- reactive({
        data_source <- input$data_source
        data_to_use <- NULL # Initialize
        
        if (data_source == "dynamic") {
          if (!is.null(shared$analysis_results$table)) {
            data_to_use <- shared$analysis_results$table
            # Ensure it's a data.table for consistency
            if (!is.data.table(data_to_use)) {
              data_to_use <- as.data.table(data_to_use)
            }
          } else {
            showNotification("No dynamic analysis results available yet.", type = "warning")
            return(NULL) # Return NULL if no dynamic data
          }
        } else { # data_source == "file"
          req(input$Pre_existing_panvaR_results_path) # Require file input if source is 'file'
          file_info <- shinyFiles::parseFilePaths(c(Home = fs::path_home()), input$Pre_existing_panvaR_results_path)
          
          if (nrow(file_info) > 0) {
            file_path <- file_info$datapath[1]
            tryCatch({
              data_to_use <- data.table::fread(file_path)
              if (nrow(data_to_use) == 0) {
                showNotification("Loaded file contains no data.", type = "warning")
                return(NULL) # Return NULL if file is empty
              }
            }, error = function(e) {
              showNotification(paste("Error loading file:", e$message), type = "error")
              return(NULL) # Return NULL on error
            })
          } else {
            return(NULL) # Return NULL if file path parsing fails
          }
        }
        
        # --- Data Validation ---
        if (!is.null(data_to_use)) {
          required_cols <- c("CHROM", "BP", "Pvalues", "LD", "Type", "IMPACT", "GENE", "EFFECT", "AA", "REF", "ALT", "final_weight")
          missing_cols <- setdiff(required_cols, names(data_to_use))
          if(length(missing_cols) > 0) {
            showNotification(paste("Loaded data is missing required columns:", paste(missing_cols, collapse=", ")), type = "error")
            return(NULL)
          }
          # Add more specific type checks if necessary
        }
        # --- End Data Validation ---
        
        return(data_to_use)
      })
      
      
      # Update UI elements reactively based on loaded data
      observeEvent(my_data(), {
        current_data <- my_data() # Get the data
        req(current_data)       # Require data to proceed
        
        # --- Start UI Update ---
        # Safely update selection inputs, handling potential NULLs/NAs
        safe_update_selectize <- function(inputId, choices) {
          valid_choices <- sort(unique(na.omit(choices))) # Get unique, non-NA choices
          updateSelectizeInput(session, inputId, choices = valid_choices, selected = character(0))
        }
        
        safe_update_selectize("selected_genes", current_data$GENE)
        safe_update_selectize("selected_effect_types", current_data$EFFECT)
        safe_update_selectize("selected_amino_acid", current_data$AA)
        safe_update_selectize("selected_REF_types", current_data$REF)
        safe_update_selectize("selected_ALT_types", current_data$ALT)
        
        # --- Update Numeric Ranges ---
        # BP range relative to tag SNP
        tag_snp_bp <- tryCatch({
          current_data %>% filter(Type == "tag_snp") %>% pull(BP) %>% unique() %>% head(1) %>% as.numeric()
        }, error = function(e) NULL)
        
        if(!is.null(tag_snp_bp) && !is.na(tag_snp_bp)) {
          min_BP <- min(current_data$BP, na.rm = TRUE)
          max_BP <- max(current_data$BP, na.rm = TRUE)
          
          # Check for valid min/max BP
          if(is.finite(min_BP) && is.finite(max_BP)) {
            bp_LHS_max <- tag_snp_bp - min_BP
            bp_RHS_max <- max_BP - tag_snp_bp
            
            # Ensure max values are non-negative
            bp_LHS_max <- max(0, bp_LHS_max, na.rm = TRUE)
            bp_RHS_max <- max(0, bp_RHS_max, na.rm = TRUE)
            
            updateNumericInput(session, "BP_LHS", max = bp_LHS_max, value = bp_LHS_max) # Default to max range
            updateNumericInput(session, "BP_RHS", max = bp_RHS_max, value = bp_RHS_max) # Default to max range
          } else {
            # Handle cases where min/max BP aren't finite (e.g., all NA)
            updateNumericInput(session, "BP_LHS", max = 0, value = 0)
            updateNumericInput(session, "BP_RHS", max = 0, value = 0)
          }
        } else {
          # Reset or disable BP range if no tag SNP found
          updateNumericInput(session, "BP_LHS", value = 0, max = 0)
          updateNumericInput(session, "BP_RHS", value = 0, max = 0)
          # Optionally disable inputs here using shinyjs
        }
        
        
        # P-value range
        min_pval <- min(current_data$Pvalues, na.rm = TRUE)
        max_pval <- max(current_data$Pvalues, na.rm = TRUE)
        if(is.finite(min_pval) && is.finite(max_pval)) {
          updateNumericInput(session, "pvalue_min", min = floor(min_pval), value = floor(min_pval))
          updateNumericInput(session, "pvalue_max", max = ceiling(max_pval), value = ceiling(max_pval))
        } else {
          # Handle cases where min/max Pvalues aren't finite
          updateNumericInput(session, "pvalue_min", value = 0)
          updateNumericInput(session, "pvalue_max", value = 30) # Reset to default or suitable max
        }
        # --- End UI Update ---
        
      }, ignoreNULL = TRUE, ignoreInit = TRUE) # ignoreNULL=T prevents running when my_data() is NULL
      
      
      # Filter data reactively based on UI inputs
      filtered_data <- reactive({
        req(my_data()) # Require base data
        data_in <- my_data()
        
        # Apply filtering using the helper function
        result <- load_and_filter_module(data_in, input)
        
        # --- Add Notification for "Tag SNP Only" Scenario ---
        if (!is.null(result) && nrow(result) > 0) {
          # Check if only the tag SNP row(s) remain after filtering
          # This is true if all rows have Type 'tag_snp' OR if all 'final_weight' are NA
          is_only_tag_after_filter <- all(result$Type == "tag_snp") || all(is.na(result$final_weight))
          
          if (is_only_tag_after_filter && nrow(result) > 0) {
            # Check if the original data had more than just the tag SNP
            original_had_candidates <- any(data_in$Type != "tag_snp", na.rm = TRUE)
            
            if(original_had_candidates) {
              # Only show notification if candidates were present initially but filtered out
              showNotification(
                "Applied filters resulted in only the Tag SNP remaining.",
                type = "warning",
                duration = 7
              )
            } else {
              # If original data also had only tag SNP, a different message might be appropriate
              # or the message from the analysis run might suffice.
              # We can also check if the original data triggering this reactive had only tag SNP
              if(all(data_in$Type == "tag_snp") || all(is.na(data_in$final_weight))) {
                showNotification(
                  "Displaying Tag SNP only (no candidates found in analysis).",
                  type = "info",
                  duration = 7
                )
              }
            }
            
          }
        } else if (is.null(result) || nrow(result) == 0) {
          showNotification("No data matches the current filter criteria.", type = "warning")
        }
        # --- End Notification ---
        
        return(result)
      })
      
      
      # Render data table with error handling and specific messages
      output$dataTable <- renderDT({
        data_to_display <- filtered_data() # Get the filtered data
        
        # Case 1: No data at all (either initially or after filtering)
        if (is.null(data_to_display) || nrow(data_to_display) == 0) {
          return(datatable(data.frame(Message = "No data available based on current filters."), options = list(dom = 't'), rownames = FALSE)) # Simple message table
        }
        
        # Case 2: Data exists, render the table
        datatable(data_to_display,
                  options = list(
                    pageLength = 10, # Shorter default page length
                    scrollX = TRUE,
                    scrollY = "350px", # Adjust height as needed
                    #scroller = TRUE, # Scroller can sometimes conflict with Buttons
                    deferRender = TRUE,
                    scrollCollapse = TRUE,
                    dom = 'Bfrtip', # Ensure Buttons extension is active
                    buttons = c('copy', 'csv', 'excel') # Standard export buttons
                  ),
                  filter = 'top', # Enable column filters
                  rownames = FALSE,
                  extensions = c('Buttons'), # Explicitly list extensions
                  selection = 'single' # Allow single row selection if needed later
        )
        
      })
      
      # Render plotly with error handling and specific messages
      output$data_plotly <- renderPlotly({
        data_to_plot <- filtered_data() # Get the filtered data
        
        # Case 1: No data to plot
        if (is.null(data_to_plot) || nrow(data_to_plot) == 0) {
          return(
            plot_ly() %>%
              layout(title = "No Data to Display", xaxis = list(visible = FALSE), yaxis = list(visible = FALSE)) %>%
              add_annotations(
                text = "No data available based on current filters.",
                x = 0.5, y = 0.5, xref = "paper", yref = "paper", showarrow = FALSE, font = list(size=14)
              )
          )
        }
        
        # Case 2: Only Tag SNP(s) remain
        # Check if all rows are tag_snp OR if all final_weights are NA
        is_only_tag <- all(data_to_plot$Type == "tag_snp") || all(is.na(data_to_plot$final_weight))
        
        if (is_only_tag) {
          return(
            plot_ly() %>%
              layout(title = "Tag SNP Only", xaxis = list(visible = FALSE), yaxis = list(visible = FALSE)) %>%
              add_annotations(
                text = "Only the Tag SNP remains after filtering.\nCannot generate comparative plot.",
                x = 0.5, y = 0.5, xref = "paper", yref = "paper", showarrow = FALSE, font = list(size=14)
              )
          )
        }
        
        # Case 3: Data available for plotting
        tryCatch({
          # Pass Bonferroni parameters from input if needed dynamically
          panvar_plotly_function(
            panvar_results_table = data_to_plot,
            nrows_in_gwas = input$total_rows,         # Pass from input
            pvalue_threshold = input$pvalue_threshold # Pass from input
          )
        }, error = function(e) {
          # Generic error plot
          plot_ly() %>%
            layout(title = "Plotting Error") %>%
            add_annotations(
              text = paste("Error creating plot:", e$message),
              x = 0.5, y = 0.5, xref = "paper", yref = "paper", showarrow = FALSE
            )
        })
        
      }) # End renderPlotly
      
    } # End function(input, output, session)
  ) # End moduleServer
} # End output_dashboard_Server
# --- END SERVER ---