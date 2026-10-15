#' @title Process Pairwise Correlations
#'
#' @description Internal function for calculating pairwise correlations. Optionally, may
#' return a mapping between correlation statistic names (e.g., "cor_var1_var2")
#' and the corresponding variable names (e.g., "var1" and "var2").
#'
#'
#' @param data Data frame used to calculate pairwise correlations.
#' @param return_map Logical, indicating whether to return a mapping
#' of correlation statistic names to the names of the corresponding variable pairs. Default is FALSE.
#'
#' @return If \code{return_map = FALSE}, a named numeric vector of pairwise
#'   correlations.
#'
#'   If \code{return_map = TRUE}, a list containing:
#'   \itemize{
#'     \item \code{cor_data}: A named numeric vector of pairwise correlations.
#'     \item \code{cor_map}: A data frame mapping correlation statistic names
#'       to the names of the corresponding variable pairs.
#'   }
pairwise.correlation <- function(
    data,
    return_map = FALSE
    )
{
  # Correlation matrix
  cor_mat <- cor(
    data,
    use = "pairwise.complete.obs"
  )

  # Find the row and column indices for upper triangle elements
  upper_indices <- which(
    upper.tri(cor_mat),
    arr.ind = TRUE
  )

  # Sort the order so it becomes correlations between x1 and other variables, then x2 with other variables
  upper_indices <- upper_indices[
    order(upper_indices[, "row"],
          upper_indices[, "col"]),
    ,
   drop = FALSE
   ]

  # List of row (first variable) names in pairwise correlations
  var1 <- rownames(cor_mat)[upper_indices[, "row"]]
  # List of column (second variable) names in pairwise correlations
  var2 <- colnames(cor_mat)[upper_indices[, "col"]]
  # Define variable name for correlations
  cor_names <- paste0(
    "cor_",
    var1,
    "_",
    var2
  )
  # Match the correlation values with correlation variable names
  cor_data <- setNames(
    cor_mat[upper_indices],
    cor_names
  )

  if (!return_map) {
    return(cor_data)
  }

  # May also return a map, matching correlation variable name with raw variable names
  list(
    cor_data = cor_data,
    cor_map = tibble::tibble(
      cor_name = cor_names,
      var1 = var1,
      var2 = var2
    )
  )
}
