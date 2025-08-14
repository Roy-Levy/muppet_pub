# Turns off the option to open the HTML file automatically
   mcmcplot.mod <- function(...) {
    # Save the current browser options
    old_browser <- getOption("browser")

    # Set a dummy browser function that does nothing
    options(browser = function(...) {})

    # Call the original mcmcplot function
    result <- mcmcplots::mcmcplot(...)

    # Restore the original browser option
    options(browser = old_browser)

    # Return the result
    return(result)
  }
