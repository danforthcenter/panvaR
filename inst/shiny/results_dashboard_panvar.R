

# Define the User Interface (UI)
plot_ui <- function(id){
  ns <- NS(id)
  tagList(
    
    # Application title
    titlePanel("Plot Results"),
    # Sidebar layout with input widgets
    sidebarLayout(
      sidebarPanel(
        width = 3,
        
        # general options
        div(class = "input-section",
            h4("General options"),
        # Window 
        div(
          style = "display: flex; align-items: center; gap: 10px;",
          numericInput(ns("window"), "Window around tag SNP in kilobases:", value = 50, min = 0, step = 5)
        ),
        hr()
        ),
        
        # Manhattan options
        div(class = "input-section",
          h4("Manhattan options"),
          
          # plot.r2.thresh 
          div(
            style = "display: flex; align-items: center; gap: 10px;",
            numericInput(ns("plot.r2.thresh"), "LD threshold:", value = .2, min = 0, step = .1)
          ),
          # unplotted.alpha 
          div(
            style = "display: flex; align-items: center; gap: 10px;",
            checkboxInput(ns("remove.low.ld.points"), "Display low LD points?", TRUE)
          ),
          
          # sig.line
          div(
            style = "display: flex; align-items: center; gap: 10px;",
            numericInput(ns("unplotted.alpha"), "Significance line:", value = 6)
          ),
          
          # qualitative.annotation
          div(
            style = "display: flex; align-items: center; gap: 10px;",
            textInput(ns("qualitative.annotation"), "Qualitative variable (shape)", value = NULL)
          ),
          # quantitative.annotation
          div(
            style = "display: flex; align-items: center; gap: 10px;",
            textInput(ns("quantitative.annotation"), "Quantitative variable (fill)", value = NULL)
          ),
          
          hr(),
          ),
        
        # Annotation options
        div(class = "input-section",
          h4("Annotation options"),
          
          # annotation.point.variable
          div(
            style = "display: flex; align-items: center; gap: 10px;",
            textInput(ns("annotation.point.variable"), "Annotation point variable (fill)", value = NULL)
          )
          
          )
      ),
      
      # Main panel for displaying the plot
      mainPanel(
        tabsetPanel(
          id = "tabs",
          tabPanel("Plot", plotlyOutput(ns("panvar_plot"))),
          tabPanel("Genes Table", dataTableOutput(ns("table_genes"))),
          tabPanel("SNPs Table", DTOutput(ns("table_snps")))
        )
      )
    )
  )
}

# Define the Server logic
plot_server <- function(id, shared) {
  moduleServer(
    id,
    function(input, output, session) {

      # set to null if empty for text inputs
      # input$qualitative.annotation <- ifelse(input$qualitative.annotation == "", NULL, input$qualitative.annotation)
      # input$quantitative.annotation <- ifelse(input$quantitative.annotation == "", NULL, input$quantitative.annotation)
      # input$annotation.point.variable <- ifelse(input$annotation.point.variable == "", NULL, input$annotation.point.variable)
      
      
      # panvar plot
      output$panvar_plot <- renderPlotly({
        if (input$qualitative.annotation == ""){
          # qc_qualitative.annotation <- ifelse(input$qualitative.annotation == "", NULL, input$qualitative.annotation)
          qc_qualitative.annotation <- NULL
        } else {
          qc_qualitative.annotation <- input$qualitative.annotation
        }
        if (input$quantitative.annotation == ""){
          # qc_quantitative.annotation <- ifelse(input$quantitative.annotation == "", NULL, input$quantitative.annotation)
          qc_quantitative.annotation <- NULL
        } else {
          qc_quantitative.annotation <- input$quantitative.annotation
        }
        if (input$annotation.point.variable == ""){
          # qc_annotation.point.variable <- ifelse(input$annotation.point.variable == "", NULL, input$annotation.point.variable)
          qc_annotation.point.variable <- NULL
        } else {
          qc_annotation.point.variable <- input$annotation.point.variable
        }
        # # if(input$quantitative.annotation == "")
        #   qc_quantitative.annotation <- ifelse(input$quantitative.annotation == "", NULL, input$quantitative.annotation)
        # # if(input$annotation.point.variable == "")
        #   qc_annotation.point.variable <- ifelse(input$annotation.point.variable == "", NULL, input$annotation.point.variable)
        
        plotly_panvar(shared$analysis_results,
                    window = input$window,
                    plot.r2.thresh = input$plot.r2.thresh,
                    # unplotted.alpha = input$unplotted.alpha,
                    remove.low.ld.points = !input$remove.low.ld.points,
                    sig.line = input$sig.line,
                    qualitative.annotation = qc_qualitative.annotation,
                    quantitative.annotation = qc_quantitative.annotation,
                    annotation.point.variable = qc_annotation.point.variable,
                    pvals.in.log = F)
        
      })
      
      # genes table
      output$table_genes <- renderDataTable({
        datatable(shared$analysis_results$anno)
      })
      
      # snps table
      output$table_snps <- renderDataTable({
        datatable(shared$analysis_results$gwas)
      })
    }
  )
}