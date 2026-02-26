

# Define the User Interface (UI)
plot_ui <- function(id){
  ns <- NS(id)
  tagList(
    # Application title
    titlePanel("Interactive Plot with User Inputs"),
    
    # Sidebar layout with input widgets
    sidebarLayout(
      sidebarPanel(
        # Input 1: Slider for number of observations
        sliderInput(ns("n_obs"),
                    "Number of observations:",
                    min = 10,
                    max = 100,
                    value = 50),
        
        # Input 2: Select box for variable to plot on x-axis
        selectInput(ns("x_var"),
                    "X Variable:",
                    choices = names(iris),
                    selected = "Sepal.Length"),
        
        # Input 3: Select box for variable to plot on y-axis
        selectInput(ns("y_var"),
                    "Y Variable:",
                    choices = names(iris),
                    selected = "Sepal.Width")
      ),
      
      # Main panel for displaying the plot
      mainPanel(
        # plotOutput(ns("distPlot"))
        plotOutput(ns("panvar_plot"))
        
      )
    )
  )
}

# Define the Server logic
plot_server <- function(id, shared) {
  moduleServer(
    id,
    function(input, output, session) {
      
      # panvar plot
      output$panvar_plot <- renderPlot({
        plot_panvar(shared$analysis_results,
                    window = 100,
                    sig.line = 6,
                    pvals.in.log = F)
        
      })
      
      # Create the plot output
      output$distPlot <- renderPlot({
        # Ensure inputs are available before plotting (optional but good practice)
        req(input$n_obs, input$x_var, input$y_var) 
        
        # Generate the plot based on user inputs
        ggplot(data = head(iris, input$n_obs), 
               aes_string(x = input$x_var, y = input$y_var, color = "Species")) +
          geom_point() +
          labs(title = paste("Plot of", input$y_var, "vs", input$x_var))
      })
    }
  )
}