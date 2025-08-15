# Turns off the option to open the HTML file automatically
#' Plot MCMC draws
#'
#' @description A modified version of the mcmcplot() function from the mcmcplots. Modified to not automatically open the html file that is created. This version creates and saves various plots and an HTML file collecting them (but does not automatically open the HTML file.)

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
