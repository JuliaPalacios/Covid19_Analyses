fluidPage(
  # Application title
  # titlePanel("Palacios COVID"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      selectInput(
        inputId = "selectedCountry",
        label = strong("Country"),
        choices = sort(names(trees)),
        selected = "China"
      )
    ),

    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("tree_plot"),
      plotOutput("eps_plot"),
      plotOutput("eps_plot_ps"),
      plotOutput("case_plot")
    )
  )
)
