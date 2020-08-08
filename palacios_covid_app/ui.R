fluidPage(
  fluidRow(style="display:flex; align-items:center; margin-bottom:20px",
    column(8, style="padding-left: 20px", withTags({
      div(
        h3("Phylodynamics of SARS-CoV-2"),
        p("Maintained by the Palacios Lab at Stanford. Check out our",
          a(href="https://github.com/juliapalacios/Covid19_Analyses", 
            "GitHub repository"),
          "for more information on these analyses. Pre-processing relies on",
          a(href="https://nextstrain.org/", "Nextstrain's"),
          a(href="https://github.com/nextstrain/ncov", "ncov project."),
          "We are thankful for their effort and ",
          "for making their work publicly available."),
        p("Enabled by data from",
          img(src="https://www.gisaid.org/fileadmin/gisaid/img/schild.png",
              alt="GISAID", style="height:1.5em"))
      )
    })),
    column(4, style="padding-top: 10px; padding-right: 40px", 
      wellPanel(style="margin-bottom:0", selectInput(
        inputId = "selectedCountry",
        label = strong("Country"),
        choices = sort(names(countries)),
        selected = "China"
      ))
    )
  ),
  fluidRow(
    column(12, withSpinner(plotOutput("tree_plots")))
  ),
  fluidRow(
    column(12, withSpinner(plotOutput("case_plot")))
  )
)
