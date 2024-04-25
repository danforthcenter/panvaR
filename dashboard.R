library(shiny)
library(bslib)
library(data.table)
library(ggplot2)
library(ggExtra)

ui <- page_sidebar(
  sidebar = sidebar(
    fileInput("file", "Choose a file", accept = c(".csv", ".tsv", ".txt")),
    selectInput("xvar", "X variable", choices = NULL),
    selectInput("yvar", "Y variable", choices = NULL),
    checkboxGroupInput(
      "species", "Filter by species",
      choices = NULL,
      selected = NULL
    ),
    hr(), # Add a horizontal rule
    checkboxInput("by_species", "Show species", TRUE),
    checkboxInput("show_margins", "Show marginal plots", TRUE),
    checkboxInput("smooth", "Add smoother"),
  ),
  plotOutput("scatter")
)

server <- function(input, output, session) {
  df <- reactive({
    req(input$file)
    fread(input$file$datapath)
  })

  observe({
    req(df())
    updateSelectInput(session, "xvar", choices = names(df()[, sapply(.SD, is.numeric), .SDcols = -"Species"]))
    updateSelectInput(session, "yvar", choices = names(df()[, sapply(.SD, is.numeric), .SDcols = -"Species"]))
    updateCheckboxGroupInput(session, "species", choices = unique(df()$Species))
  })

  subsetted <- reactive({
    req(input$species)
    df()[Species %in% input$species]
  })

  output$scatter <- renderPlot({
    req(input$xvar, input$yvar)
    p <- ggplot(subsetted(), aes(!!sym(input$xvar), !!sym(input$yvar))) + list(
      theme(legend.position = "bottom"),
      if (input$by_species) aes(color = Species),
      geom_point(),
      if (input$smooth) geom_smooth()
    )

    if (input$show_margins) {
      margin_type <- if (input$by_species) "density" else "histogram"
      p <- ggExtra::ggMarginal(p, type = margin_type, margins = "both",
        size = 8, groupColour = input$by_species, groupFill = input$by_species)
    }

    p
  }, res = 100)
}

shinyApp(ui, server)
