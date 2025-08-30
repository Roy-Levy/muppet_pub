#' @importFrom coda mcmc.list
#'
#' @title Converts a dataframe to an mcmc.list object
#'
#' @description This is an internal function, and is not intended to be called by users directly.

#' @param df The dataframe to be converted to an mcmc.list object
#' @param chain_col Do not alter
#' @param iter_col Do not alter
#'
#' @return An mcmc.list object
#'
#'
df_to_mcmc_list <- function(df, chain_col = "Chain.number", iter_col = "Iteration.number") {
  # # Check if coda package is installed and load it
  # if (!requireNamespace("coda", quietly = TRUE)) {
  #   install.packages("coda")
  #   library(coda)
  # } else {
  #   library(coda)
  # }

  # Ensure the required columns exist
  if (!(chain_col %in% colnames(df))) {
    stop(paste("Chain column", chain_col, "not found in the CSV file"))
  }
  if (!(iter_col %in% colnames(df))) {
    stop(paste("Iteration column", iter_col, "not found in the CSV file"))
  }

  # Identify the unique chains
  chains <- unique(df[[chain_col]])
  n_chains <- length(chains)

  # Create a list to store each chain's MCMC object
  chain_list <- list()

  for (i in 1:n_chains) {
    # Subset df for current chain
    chain_df <- df[df[[chain_col]] == chains[i], ]

    # Sort by iteration
    chain_df <- chain_df[order(chain_df[[iter_col]]), ]

    # Extract parameter columns (excluding chain and iteration columns)
    param_cols <- setdiff(colnames(chain_df), c(chain_col, iter_col))

    # Create an mcmc object with just the parameter columns
    chain_list[[i]] <- coda::mcmc(chain_df[, param_cols])
  }

  # Convert to mcmc.list
  mcmc_list_object <- coda::mcmc.list(chain_list)

  return(mcmc_list_object)
}
