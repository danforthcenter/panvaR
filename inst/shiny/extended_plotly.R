library(ggplot2)
library(plotly)
library(data.table)
library(shiny)
library(dplyr)
library(DT)

# Enhanced panvar plot function with gene name tooltips and alpha control
panvar_plot <- function(reports_table, total_rows, pvalue_threshold, point_size = 3, alpha_base = 0.7) {
  hline_value = -log10(pvalue_threshold / total_rows)
  
  panvar_plots <- reports_table %>%
    ggplot(aes(x = BP, 
               y = Pvalues, 
               color = IMPACT,
               alpha = LD,
               text = paste("Gene:", GENE,
                            "<br>BP:", BP,
                            "<br>P-value:", Pvalues,
                            "<br>Impact:", IMPACT,
                            "<br>LD:", LD))) +
    theme_classic() +
    geom_point(size = point_size) +
    geom_hline(aes(yintercept = hline_value, 
                   linetype = "Bonferroni threshold"), 
               color = "purple") +
    scale_alpha(range = c(0.3 * alpha_base, alpha_base)) +
    labs(x = "Base Position",
         y = "-log10(P-value)") +
    theme(legend.position = "right",
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14))
  
  return(panvar_plots)
}

# Add the missing renderPlotFunction
renderPlotFunction <- function(plot_obj, height = 600) {
  if ("ggplot" %in% class(plot_obj)) {
    ggplotly(plot_obj, tooltip = "text", height = height) %>%
      layout(hoverlabel = list(bgcolor = "white"),
             legend = list(x = 1.05, y = 0.9))
  } else if ("plotly" %in% class(plot_obj)) {
    plot_obj
  } else {
    plot_ly() %>%
      add_annotations(text = "Plot format not recognized",
                      showarrow = FALSE)
  }
}

# Enhanced load and filter module with additional filtering options
load_and_filter_module <- function(data, input) {
  filtered_data <- data %>%
    filter(
      IMPACT %in% input$impact,
      LD >= as.numeric(input$ld_min) & LD <= as.numeric(input$ld_max),
      BP >= as.numeric(input$bp_min) & BP <= as.numeric(input$bp_max),
      Pvalues >= as.numeric(input$pvalue_min) & Pvalues <= as.numeric(input$pvalue_max)
    )
  
  # Additional filters
  if (!is.null(input$selected_genes) && length(input$selected_genes) > 0) {
    filtered_data <- filtered_data %>% filter(GENE %in% input$selected_genes)
  }
  
  if (!is.null(input$selected_effects) && length(input$selected_effects) > 0) {
    filtered_data <- filtered_data %>% filter(EFFECT %in% input$selected_effects)
  }
  
  if (!is.null(input$selected_aa) && length(input$selected_aa) > 0) {
    filtered_data <- filtered_data %>% filter(AA %in% input$selected_aa)
  }
  
  if (!is.null(input$selected_ref) && length(input$selected_ref) > 0) {
    filtered_data <- filtered_data %>% filter(REF %in% input$selected_ref)
  }
  
  if (!is.null(input$selected_alt) && length(input$selected_alt) > 0) {
    filtered_data <- filtered_data %>% filter(ALT %in% input$selected_alt)
  }
  
  return(filtered_data)
}

# Define UI with enhanced controls
ui <- fluidPage(
  titlePanel("Enhanced Interactive GWAS Plot"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      # Impact selection
      checkboxGroupInput("impact", 
                         "Select Impact Levels:",
                         choices = c("HIGH", "MODERATE", "MODIFIER"),
                         selected = c("HIGH", "MODERATE", "MODIFIER")),
      
      # Gene selection
      selectizeInput("selected_genes",
                     "Select Genes:",
                     choices = NULL,
                     multiple = TRUE),
      
      # Effect type selection
      selectizeInput("selected_effects",
                     "Select Effect Types:",
                     choices = NULL,
                     multiple = TRUE),
      
      # AA selection
      selectizeInput("selected_aa",
                     "Select AA Changes:",
                     choices = NULL,
                     multiple = TRUE),
      
      # REF selection
      selectizeInput("selected_ref",
                     "Select REF Alleles:",
                     choices = NULL,
                     multiple = TRUE),
      
      # ALT selection
      selectizeInput("selected_alt",
                     "Select ALT Alleles:",
                     choices = NULL,
                     multiple = TRUE),
      
      # LD range as numeric inputs
      fluidRow(
        column(6,
               numericInput("ld_min",
                            "LD Min:",
                            value = 0,
                            min = 0,
                            max = 1,
                            step = 0.01)),
        column(6,
               numericInput("ld_max",
                            "LD Max:",
                            value = 1,
                            min = 0,
                            max = 1,
                            step = 0.01))
      ),
      
      # Base position range
      fluidRow(
        column(6,
               numericInput("bp_min",
                            "BP Min:",
                            value = 6500000)),
        column(6,
               numericInput("bp_max",
                            "BP Max:",
                            value = 6800000))
      ),
      
      # P-value range as numeric inputs
      fluidRow(
        column(6,
               numericInput("pvalue_min",
                            "P-value Min:",
                            value = 0)),
        column(6,
               numericInput("pvalue_max",
                            "P-value Max:",
                            value = 30))
      ),
      
      # Bonferroni correction inputs
      fluidRow(
        column(6,
               numericInput("pvalue_threshold",
                            "P-value threshold:",
                            value = 0.05)),
        column(6,
               numericInput("total_rows",
                            "Total tests:",
                            value = 450000000))
      ),
      
      # Visual controls
      numericInput("pointSize",
                   "Point Size:",
                   value = 3,
                   min = 1,
                   max = 10),
      
      numericInput("alphaValue",
                   "Transparency:",
                   value = 0.7,
                   min = 0.1,
                   max = 1,
                   step = 0.1),
      
      downloadButton("downloadData", "Download Filtered Data")
    ),
    
    mainPanel(
      width = 9,
      
      # Summary statistics panels
      fluidRow(
        column(12,
               wellPanel(
                 style = "background-color: #f5f5f5; padding: 10px; margin-bottom: 20px;",
                 h4("Summary Statistics:"),
                 fluidRow(
                   column(4, uiOutput("impactCounts")),
                   column(4, uiOutput("geneCounts")),
                   column(4, uiOutput("aaCounts"))
                 )
               ))
      ),
      
      plotlyOutput("filteredPlot", height = "400px"),
      
      tags$div(style = "margin-top: 150px;"),
      
      DTOutput("dataTable")
    )
  )
)

# Enhanced server function
server <- function(input, output, session, data) {
  # Update selection choices
  observe({
    updateSelectizeInput(session, "selected_genes",
                         choices = sort(unique(data$GENE)),
                         selected = character(0))
    
    updateSelectizeInput(session, "selected_effects",
                         choices = sort(unique(data$EFFECT)),
                         selected = character(0))
    
    updateSelectizeInput(session, "selected_aa",
                         choices = sort(unique(data$AA)),
                         selected = character(0))
    
    updateSelectizeInput(session, "selected_ref",
                         choices = sort(unique(data$REF)),
                         selected = character(0))
    
    updateSelectizeInput(session, "selected_alt",
                         choices = sort(unique(data$ALT)),
                         selected = character(0))
  })
  
  # Reactive filtered data
  filtered_data <- reactive({
    load_and_filter_module(data, input)
  })
  
  # Impact counts summary
  output$impactCounts <- renderUI({
    counts <- filtered_data() %>%
      group_by(IMPACT) %>%
      summarise(count = n(), .groups = 'drop')
    
    tags$div(
      style = "display: flex; flex-direction: column; gap: 10px;",
      h5("Impact Distribution:"),
      lapply(1:nrow(counts), function(i) {
        tags$div(
          style = paste0("text-align: center; padding: 5px; border-radius: 5px; ",
                         "background-color: ", switch(counts$IMPACT[i],
                                                      "HIGH" = "#ff9999",
                                                      "MODERATE" = "#99ff99",
                                                      "MODIFIER" = "#9999ff")),
          p(style = "margin: 0;", 
            paste0(counts$IMPACT[i], ": ", counts$count[i]))
        )
      })
    )
  })
  
  # Gene counts summary
  output$geneCounts <- renderUI({
    counts <- filtered_data() %>%
      group_by(GENE) %>%
      summarise(count = n(), .groups = 'drop') %>%
      arrange(desc(count)) %>%
      head(5)
    
    tags$div(
      style = "display: flex; flex-direction: column; gap: 10px;",
      h5("Top 5 Genes:"),
      lapply(1:nrow(counts), function(i) {
        tags$div(
          style = "text-align: center; padding: 5px; border-radius: 5px; background-color: #e6e6e6;",
          p(style = "margin: 0;", 
            paste0(counts$GENE[i], ": ", counts$count[i]))
        )
      })
    )
  })
  
  # AA counts summary
  output$aaCounts <- renderUI({
    counts <- filtered_data() %>%
      group_by(AA) %>%
      summarise(count = n(), .groups = 'drop') %>%
      arrange(desc(count)) %>%
      head(5)
    
    tags$div(
      style = "display: flex; flex-direction: column; gap: 10px;",
      h5("Top 5 AA Changes:"),
      lapply(1:nrow(counts), function(i) {
        tags$div(
          style = "text-align: center; padding: 5px; border-radius: 5px; background-color: #e6e6e6;",
          p(style = "margin: 0;", 
            paste0(counts$AA[i], ": ", counts$count[i]))
        )
      })
    )
  })
  
  # Render plot with alpha control
  output$filteredPlot <- renderPlotly({
    req(filtered_data())
    plot <- panvar_plot(filtered_data(), 
                        input$total_rows,
                        input$pvalue_threshold,
                        point_size = input$pointSize,
                        alpha_base = input$alphaValue)
    renderPlotFunction(plot)
  })
  
  output$dataTable <- renderDT({
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
  
  # Download handler
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("filtered_gwas_data_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(filtered_data(), file, row.names = FALSE)
    }
  )
  
  # Update BP range based on data
  observe({
    updateNumericInput(session, "bp_min",
                       value = min(data$BP))
    updateNumericInput(session, "bp_max",
                       value = max(data$BP))
  })
}

# Run the application
runApp <- function(data) {
  required_cols <- c("BP", "Pvalues", "IMPACT", "LD", "GENE", "EFFECT", "AA", "REF", "ALT")
  missing_cols <- setdiff(required_cols, colnames(data))
  
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", 
               paste(missing_cols, collapse = ", ")))
  }
  
  shinyApp(ui = ui,
           server = function(input, output, session) {
             server(input, output, session, data)
           })
}

sample_table <- fread("~/work_in_progress/panvar_local_run/setaria_sample_table.tsv")
runApp(sample_table)
