


function(input, output) {
  output$case_plot <- renderPlot({
    country <- reformat_country(input$selectedCountry)
    cases <- cases_all[cases_all$location == country, ]

    # use as.Date() instead of date() to avoid tz warning
    dates <- decimal_date(as.Date(cases$date))
    date.labs <- seq(min(dates), max(dates), by = 4 / 365)

    # total cases will exceed ylim in the plot, you can add
    # ylim=c(0, max(cases$total_cases)) but it will dwarf everything else
    plot(dates, cases$new_cases,
      type = "l", xaxt = "n", main = "Case Data",
      ylab = "Count", xlab = "", lwd = 2
    ) # ---------------------------- CASES PLOT
    axlabs2 <- list(
      x = date.labs, labs = format(date_decimal(date.labs), "%b-%d"),
      cexlab = .1
    )
    axis(1, at = axlabs2$x, labels = axlabs2$labs, cex.axis = 1, las = 1)
    abline(v = trees[[input$selectedCountry]]$lastdate, lty = 2)
    points(dates, cases$total_cases, type = "l",
        xaxt = "n", col = "blue", lwd = 2)
    points(dates, cases$total_deaths, type = "l",
        xaxt = "n", col = "red", lwd = 2)

    # TODO don't just do y = 4000 for this because it won't work in every plot
    # TODO need to plot legend outside of plot because we don't know where
    # we can put it for x value anyway.
    legend_y <- as.integer(max(cases$new_cases) * 0.9)
    legend(2020.15,
      y = legend_y, bty = "n", c("New cases", "Total cases", "Deaths"),
      col = c("black", "blue", "red"), lty = 1, lwd = 2
    )
  })

  output$tree_plot <- renderPlot({
    plot(trees[[input$selectedCountry]]$tree, show.tip.label = FALSE, cex = .3,
         main = "UPGMA Tree")
  })

  output$eps_plot <- renderPlot({
    country <- input$selectedCountry
    # Cannot abstract this out into a common function because R annoyingly
    # copies list vars rather than just passing refs
    if (is.element(country, names(bnp_cache))) {
      bnp <- bnp_cache[[country]]
    } else {
      bnp <- BNPR(trees[[country]]$tree)
      bnp_cache[[country]] <<- bnp # set global var
    }
    axlabs <- axis_label(bnp, trees[[country]]$lastdate, byy = 4 / 365)
    plot_BNPR2(bnp, axlabs = axlabs, log = "", 
               main = "Effective Population Size (EPS)")
  })

  output$eps_plot_ps <- renderPlot({
    country <- input$selectedCountry
    if (is.element(country, names(bnp_ps_cache))) {
      bnp_ps <- bnp_ps_cache[[country]]
    } else {
      bnp_ps <- BNPR_PS(trees[[country]]$tree)
      bnp_ps_cache[[country]] <<- bnp_ps # set global var
    }
    axlabs <- axis_label(bnp_ps, trees[[country]]$lastdate, byy = 4 / 365)
    plot_BNPR2(bnp_ps, axlabs = axlabs, log = "",
        main = "EPS - Preferential Sampling")
  })
}
