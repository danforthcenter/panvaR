library(shiny)
library(bslib)
library(data.table)

ui <- page_sidebar(
  sidebar = sidebar(
    fileInput("file", "Choose a file", accept = c(".csv", ".tsv", ".txt"))
  ),
  mainPanel(
    column(
      width = 8,
      tableOutput("data_table")
    ),
    column(
      width = 4,
      verbatimTextOutput("id_count"),
      verbatimTextOutput("significance_count")
    )
  )
)

server <- function(input, output, session) {
  df <- reactive({
    req(input$file)
    fread(input$file$datapath)
  })

  output$data_table <- renderTable({
    req(df())
    df()
  })

  output$id_count <- renderPrint({
    req(df())
    id_counts <- df()[, .N, by = ID]
    paste("Unique ID counts:")
    print(id_counts)
  })

  output$significance_count <- renderPrint({
    req(df())
    significance_true_count <- sum(df()$significance)
    paste("Number of TRUE in the 'significance' field:", significance_true_count)
  })

  output$plot
}

shinyApp(ui, server)
