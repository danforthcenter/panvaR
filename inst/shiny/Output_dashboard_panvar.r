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
  
  # If no tag_snp found, return NULL
  if (is.null(tag_snp)) {
    warning("No tag_snp found in data")
    return(NULL)
  }
  
  # Safely apply filters with error handling
  filtered_data <- tryCatch({
    result <- my_data
    
    # Apply numeric filters if input exists
    if (!is.null(input$LD_min) && !is.null(input$LD_max)) {
      result <- result %>%
        filter(LD >= as.numeric(input$LD_min) & LD <= as.numeric(input$LD_max))
    }
    
    if (!is.null(input$BP_LHS) && !is.null(input$BP_RHS) && !is.null(tag_snp)) {
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
    return(NULL)
  })
  
  return(filtered_data)
}


# Rijan: Ploty object generator
panvar_plotly_function <- function(panvar_results_table, nrows_in_gwas = NULL, pvalue_threshold = 0.05, point_size = 3, alpha_base = 0.7){
  
  # Rijan: What is the tag_snp
  tag_df <- panvar_results_table %>%
    filter(Type == 'tag_snp')
  
  tag_snp <- tag_df %>%
    pull(BP) %>%
    unique()
  
  # Rijan: Did the user supply a default value for the bonforoni correction?
  if(is.null(nrows_in_gwas)){
    print("You did not supply a value for the number of tests that were in the GWAS - a place holder value will be used. This is not ideal.")
    nrows_in_gwas <- 5e6
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
      size = 20 # Adjust size here (in pixels)
    ),
    showlegend = TRUE,
    name = ~IMPACT
  )
  
  panvar_plotly <- layout(panvar_plotly, title = "Interactive PanvaR plot.",
                          xaxis = list(title = "Position", titlefont = list(size = 18)),
                          yaxis = list(title = "-log[10](p-value)", titlefont = list(size = 18)))
  
  # Rijan: Function for horizontal line
  hline <- function(y = 0, color = "blue") {
    list(
      type = "line",
      x0 = 0,
      x1 = 1,
      xref = "paper", # Stretches across the plot
      y0 = y,
      y1 = y,
      line = list(color = color),
      name = "Bonferroni threshold"
    )
  }
  
  # Rijan: Function for vertical line
  vline <- function(x = 0, color = "green") {
    list(
      type = "line",
      y0 = 0,
      y1 = 1,
      yref = "paper", # Stretches across the plot
      x0 = x,
      x1 = x,
      line = list(color = color, dash = "dot"),
      name = "Tag_SNP"
    )
  }
  
  # Rijan: Add lines to the layout
  panvar_plotly <- layout(panvar_plotly,
                          title = "Interactive PanvaR plot.",
                          xaxis = list(title = "Position", titlefont = list(size = 18)),
                          yaxis = list(title = "-log[10](p-value)", titlefont = list(size = 18)),
                          shapes = list(
                            vline(tag_snp),
                            hline(bonf.cor) # Add horizontal line at y = 0.3
                          ),
                          legend = list(title = list(text = "Legend"))
                          
  )
  
  return(panvar_plotly)
  
}
# ---

# ---
# The UI logic

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
          condition = sprintf("input['%s'] === 'dynamic' || (input['%s'] === 'file' && input['%s'] !== null)", 
                              ns("data_source"), ns("data_source"), ns("Pre_existing_panvaR_results_path")),
          
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
              id = "selected_effect_types_tooltip"
            )
          ),
          bsTooltip(
            id = "selected_effect_types_tooltip",
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
              id = "selected_amino_acid_tooltip"
            )
          ),
          bsTooltip(
            id = "selected_amino_acid_tooltip",
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
              id = "selected_ALT_types_tooltip"
            )
          ),
          bsTooltip(
            id = "selected_ALT_types_tooltip",
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
              id = "selected_REF_types_tooltip"
            )
          ),
          bsTooltip(
            id = "selected_REF_types_tooltip",
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
              id = "LD_range_tooltip"
            )
          ),
          bsTooltip(
            id = "LD_range_tooltip",
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
                                    "BP on the LHS of tag_SNP:",
                                    value = 0,
                                    min = 0,
                                    max = 0,
                                    step = 1)
                ),
                column(6,
                       numericInput(ns("BP_RHS"),
                                    "BP on the RHS of tag_snp:",
                                    value = 0,
                                    min = 0,
                                    max = 0,
                                    step = 1)
                )
              )
            ),
            span(
              style = "flex-shrink: 0;",
              icon("question-circle", style = "color: green;"),
              id = "BP_range_tooltip"
            )
          ),
          bsTooltip(
            id = "BP_range_tooltip",
            title = "What BP range should the results be filtered for? The LHS value will be added to the left of the tag_snp BP and the RHS value will be added to the right of the tag_snp BP.",
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
                                    "P-value Min:",
                                    value = 0)),
                column(6,
                       numericInput(ns("pvalue_max"),
                                    "P-value Max:",
                                    value = 30))
              )
            ),
            span(
              style = "flex-shrink: 0;",
              icon("question-circle", style = "color: green;"),
              id = "Pvalue_range_tooltip"
            )
          ),
          bsTooltip(
            id = "Pvalue_range_tooltip",
            title = "What Pvalue range should the results be filtered for?",
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
                                    "P-value threshold:",
                                    value = 0.05)),
                column(6,
                       numericInput(ns("total_rows"),
                                    "Total tests:",
                                    value = 450000000))
              )
            ),
            span(
              style = "flex-shrink: 0;",
              icon("question-circle", style = "color: green;"),
              id = "Bonferroni_correction_tooltip"
            )
          ),
          bsTooltip(
            id = "Bonferroni_correction_tooltip",
            title = "What are the parameters for your bonferroni hline?",
            placement = "right",
            trigger = "hover"
          )
        )
      ),
      
      mainPanel(
        plotlyOutput(ns("data_plotly")),
        DTOutput(ns("dataTable"))
      )
    )
  )
}

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
      
      # Read data reactively with error handling
      my_data <- reactive({
        # Check data source
        if (input$data_source == "dynamic") {
          if (!is.null(shared$analysis_results)) {
            return(shared$analysis_results$table)
          } else {
            showNotification(
              "No dynamic analysis results available",
              type = "warning"
            )
            return(NULL)
          }
        } else {
          # Handle file input
          req(input$Pre_existing_panvaR_results_path)
          file_path <- shinyFiles::parseFilePaths(
            c(Home = fs::path_home()),
            input$Pre_existing_panvaR_results_path
          )$datapath
          
          if (length(file_path) > 0) {
            tryCatch({
              data <- fread(file_path)
              if (nrow(data) == 0) {
                showNotification("Loaded file contains no data", type = "warning")
                return(NULL)
              }
              return(data)
            }, error = function(e) {
              showNotification(
                paste("Error loading file:", e$message),
                type = "error"
              )
              return(NULL)
            })
          }
          return(NULL)
        }
      })
      
      # Update UI elements when data changes
      observeEvent(my_data(), {
        req(my_data())
        
        # Safely update selection inputs
        tryCatch({
          updateSelectizeInput(
            session,
            "selected_genes",
            choices = sort(unique(my_data()$GENE)),
            selected = character(0)
          )
          
          updateSelectizeInput(
            session,
            "selected_effect_types",
            choices = sort(unique(my_data()$EFFECT)),
            selected = character(0)
          )
          
          updateSelectizeInput(
            session,
            "selected_amino_acid",
            choices = sort(unique(my_data()$AA)),
            selected = character(0)
          )
          
          updateSelectizeInput(
            session,
            "selected_REF_types",
            choices = sort(unique(my_data()$REF)),
            selected = character(0)
          )
          
          updateSelectizeInput(
            session,
            "selected_ALT_types",
            choices = sort(unique(my_data()$ALT)),
            selected = character(0)
          )
          
          # Calculate and update range inputs
          tag_snp <- my_data() %>% 
            filter(Type == "tag_snp") %>%
            pull(BP) %>% 
            unique() %>% 
            as.numeric()
          
          if (!is.null(tag_snp)) {
            min_BP <- min(my_data()$BP)
            max_BP <- max(my_data()$BP)
            
            bp_LHS_max <- tag_snp - min_BP
            bp_RHS_max <- max_BP - tag_snp
            
            updateNumericInput(session, "BP_LHS", 
                               max = bp_LHS_max,
                               value = bp_LHS_max)
            
            updateNumericInput(session, "BP_RHS", 
                               max = bp_RHS_max,
                               value = bp_RHS_max)
          }
          
          # Update p-value range
          updateNumericInput(session, "pvalue_min",
                             min = min(my_data()$Pvalues),
                             value = min(my_data()$Pvalues))
          
          updateNumericInput(session, "pvalue_max",
                             max = max(my_data()$Pvalues),
                             value = max(my_data()$Pvalues))
          
        }, error = function(e) {
          showNotification(
            paste("Error updating UI elements:", e$message),
            type = "error"
          )
        })
      })
      
      # Filter data with error handling
      filtered_data <- reactive({
        req(my_data())
        result <- load_and_filter_module(my_data(), input)
        if (is.null(result)) {
          showNotification(
            "No data available after applying filters",
            type = "warning"
          )
        }
        result
      })
      
      # Render data table with error handling
      output$dataTable <- renderDT({
        req(filtered_data())
        if (is.null(filtered_data()) || nrow(filtered_data()) == 0) {
          return(datatable(data.frame(Message = "No data available")))
        }
        
        datatable(filtered_data(),
                  options = list(
                    pageLength = 25,
                    scrollX = TRUE,
                    scrollY = "400px",
                    scroller = TRUE,
                    deferRender = TRUE,
                    scrollCollapse = TRUE,
                    dom = 'Bfrtip',
                    buttons = c('copy', 'csv', 'excel')
                  ),
                  filter = 'top',
                  rownames = FALSE,
                  extensions = c('Scroller', 'Buttons'))
      })
      
      # Render plotly with error handling
      output$data_plotly <- renderPlotly({
        req(filtered_data())
        if (is.null(filtered_data()) || nrow(filtered_data()) == 0) {
          plot_ly() %>%
            add_annotations(
              text = "No data available",
              x = 0.5,
              y = 0.5,
              xref = "paper",
              yref = "paper",
              showarrow = FALSE
            )
        } else {
          tryCatch({
            panvar_plotly_function(panvar_results_table = filtered_data())
          }, error = function(e) {
            plot_ly() %>%
              add_annotations(
                text = paste("Error creating plot:", e$message),
                x = 0.5,
                y = 0.5,
                xref = "paper",
                yref = "paper",
                showarrow = FALSE
              )
          })
        }
      })
    }
  )
}