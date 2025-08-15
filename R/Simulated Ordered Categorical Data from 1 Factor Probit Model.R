# Function to simulate ordered categorical data from a factor model
# Probit model
# 1 Latent variable
# with underlying latent response formulation

library(lavaan)
library(MASS)  # For mvrnorm function

# Function to simulate categorical data from a 1-factor probit model
#' Simulate categorical data from a 1-factor probit model
#'
#' @description
#' Simulates categorical data from a 1-factor probit model. The model assumes the latent variable is normally distributed with a mean of 0. The model is akin to a latent response variable formulation where the underlying latent response variable for a person's item response is formed by the product of the loading and latent variable, and then a N(0,1) residual is added. The underlying latent variable is then discretized via the item thresholds, yielding a categorical response variable.
#'
#' @param n_obs The number of cases, subjects, participants in the dataset to be generated. Default is 500.
#' @param n_items The number of items in the dataset to be generated. Default is 5.
#' @param n_categories The number of response categories for the item responses. One value applies to all items. Default is 2.
#' @param lowest_category The value of the lowest response category (usually 0 or 1). Default is 0.
#' @param factor_loadings A vector of length n_items with factor loadings (1 per item). If no vector is provided, defaults to simulating from a Uniform distribution(.5, .8) distribution.
#' @param thresholds A matrix of thresholds. Each row in the matrix contains the thresholds for that item. There should be one fewer columns than the number of categories. If no matrix is provided, defaults to each item having thresholds equally spaced out between -2 and 2.
#' @param latent_var_sd The value of the standard deviation of the latent variable. Default is 1.
#'
#' @return A list, including: data, the simulated item response data; model, an expression of the model used for generating the data; true_loadings, the loadings used for simulating the data; true_thresholds, the thresholds used for simulating the data; underlying_latent_variable_responses, the matrix of simulated latent response variables (one per case per item).
#' @export
#'
simulate_ordinal_factor_data <- function(n_obs = 500,
                                         n_items = 5,
                                         n_categories = 2,
                                         lowest_category = 0,
                                         factor_loadings = NULL,
                                         thresholds = NULL,
                                         latent_var_sd = 1) {

  # Set default factor loadings if not provided (between 0.5 and 0.8)
  if (is.null(factor_loadings)) {
    factor_loadings <- runif(n_items, 0.5, 0.8)
  }

  # Create the factor model definition
  model_def <- "
  # Define the latent factor
  F1 =~ "

  # Add items with their loadings
  for (i in 1:n_items) {
    model_def <- paste0(model_def, factor_loadings[i], "*Y", i)
    if (i < n_items) model_def <- paste0(model_def, " + ")
  }

  # Close the model definition
  model_def <- paste0(model_def, "\n")

  # Generate continuous latent factor scores
  latent_variable_values <- rnorm(n_obs, 0, latent_var_sd)

  # Generate continuous responses before categorization
  underlying_latent_variable_responses <- matrix(0, nrow = n_obs, ncol = n_items)
  for (j in 1:n_items) {
    # Each item is determined by the latent factor plus noise
    # underlying_latent_variable_responses[, j] <- factor_loadings[j] * latent_variable_values + rnorm(n_obs, 0, sqrt(1 - factor_loadings[j]^2))
    underlying_latent_variable_responses[, j] <- factor_loadings[j] * latent_variable_values + rnorm(n_obs, 0, 1)
  }

  # Convert continuous responses to categorical
  categorical_data <- matrix(0, nrow = n_obs, ncol = n_items)

  # Set default thresholds if not provided
  if (is.null(thresholds)) {
    # Create (n_categories-1) thresholds for each item
    thresholds <- matrix(0, nrow = n_items, ncol = n_categories - 1)
    for (j in 1:n_items) {
      # Evenly spaced thresholds between -2 and 2 by default
      thresholds[j, ] <- seq(-2, 2, length.out = n_categories - 1)
    }
  }

  # Apply thresholds to create categorical data
  for (j in 1:n_items) {
    for (i in 1:n_obs) {
      categorical_data[i, j] <- lowest_category  # Default lowest category
      for (k in 1:(n_categories - 1)) {
        if (underlying_latent_variable_responses[i,j] > thresholds[j, k]) {
          categorical_data[i, j] <- k + lowest_category
        }
      }
    }
  }

  # Create a dataframe with the categorical responses
  df <- as.data.frame(categorical_data)
  colnames(df) <- paste0("Y", 1:n_items)

  # Return the data and model
  return(list(
    data = df,
    model = model_def,
    true_loadings = factor_loadings,
    true_thresholds = thresholds,
    underlying_latent_variable_responses = underlying_latent_variable_responses
  ))
}
