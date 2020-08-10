


function(input, output) {
  # tree_cache <- reactiveValues()
  # for (c in names(countries)) {
  #   tree_cache[[]] = NULL
  # }
  output$case_plot <- renderPlot({
    is_state <- is_state(input$selectedCountry)
    country <- reformat_country(input$selectedCountry)
    
    if (is_state) {
      cases <- cases_usa[cases_usa$state == country, ]
    } else {
      cases <- cases_all[cases_all$location == country, ]  
    }
    
    # use as.Date() instead of date() to avoid tz warning
    dates <- decimal_date(as.Date(cases$date))
    date.labs <- seq(min(dates), max(dates), by = 4 / 365)

    ymax <- as.integer(max(max(cases$new_cases), max(cases$total_deaths)) * 1.1)
    
    # total cases will exceed ylim in the plot, you can add
    # ylim=c(0, max(cases$total_cases)) but it will dwarf everything else
    plot(dates, cases$new_cases,
      type = "l", xaxt = "n", main = paste("Case Data -", country),
      ylab = "Count", xlab = "", lwd = 2, ylim = c(0,ymax)
    ) # ---------------------------- CASES PLOT
    axlabs2 <- list(
      x = date.labs, labs = format(date_decimal(date.labs), "%b-%d"),
      cexlab = .1
    )
    axis(1, at = axlabs2$x, labels = axlabs2$labs, cex.axis = 1, las = 1)
    # abline(v = trees[[input$selectedCountry]]$lastdate, lty = 2)
    #points(dates, cases$total_cases, type = "l",
    #    xaxt = "n", col = "blue", lwd = 2)
    points(dates, cases$total_deaths, type = "l",
        xaxt = "n", col = "red", lwd = 2)

    # TODO need to plot legend outside of plot because we don't know where
    # we can put it to stay out of the way if it's inside
    legend(2020.0,
      y = ymax, bty = "n", c("New cases", "Total deaths"),
      col = c("black", "red"), lty = 1, lwd = 2
    )
  })
  
  output$tree_plots <- renderPlot({
    country <- input$selectedCountry
    tree_meta <- compute_tree(country, as.numeric(input$selectedMu))
    
    # 2-row matrix, full-width plot on row 1, row 2 split.
    # layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
    layout(matrix(c(1,2,1,3), 2, 2, byrow = TRUE))
    par(mar=c(1,1,1,1)) # set margins
    
    # tree plot
    par(mai=c(0,0,0.5,0)) # bottom, left, top, right margin of plot
    plot(tree_meta$tree, show.tip.label = FALSE, cex = .3, main = "UPGMA Tree")
    
    # eps plot
    bnp <- BNPR(tree_meta$tree)
    axlabs <- axis_label(bnp, tree_meta$lastdate, byy = 4 / 365)
    par(mai=c(0.75,0.5,0.5,0.5))
    plot_BNPR2(bnp, axlabs = axlabs, log = "", xlab = NULL,
               main = "Effective Population Size (EPS)")
    
    # eps ps plot
    bnp_ps <- BNPR_PS(tree_meta$tree)
    axlabs <- axis_label(bnp_ps, tree_meta$lastdate, byy = 4 / 365)
    par(mai=c(0.5,0.5,0.75,0.5))
    plot_BNPR2(bnp_ps, axlabs = axlabs, log = "", xlab = NULL,
               main = "EPS - Preferential Sampling")
  })
}
