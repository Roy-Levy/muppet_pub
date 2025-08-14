# Turns off the option to open the HTML file automatically
#' Plotting function mcmcplot() from mcmcplots() package modified to not automatically open the html
#'
#' @param See mcmcplot() from mcmcplots() package.
#'
#' @return See mcmcplot() from mcmcplots() package.
#' @export
#'
#' @examples See mcmcplot() from mcmcplots() package.

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
