#' @title Process Residual Correlations
#'
#' @description Internal function for reordering pairwise residual correlations
#' according to the variable order in the observed data.
#'
#' @param res_cor_out Matrix of residual correlations extracted from an Mplus output file.
#' @param variable_list Character string containing the variable names in the observed data.
#' @param prefix Character string used as prefix of residual correlation names.
#'
#' @return A one-row data frame of pairwise residual correlations.
#'
residual.cor.reorder <- function(
    res_cor_out,
    variable_list,
    prefix
)
{
  # Turn this into a full correlation matrix (lower-triangle only in .out)
  res_cor_out_full <- as.matrix(res_cor_out)

  res_cor_out_full[upper.tri(res_cor_out_full)] <-
    t(res_cor_out_full)[upper.tri(res_cor_out_full)]

  # Reorder the columns and rows so that the variables are in the same order as in the input data
  var_list <- unlist(strsplit(variable_list, "\\s+"))
  var_idx <- match(
    tolower(var_list),
    tolower(rownames(res_cor_out_full))
  )

  res_cor_out_full <- res_cor_out_full[var_idx, var_idx]

  # Find the row and column indices for upper triangle elements
  upper_indices <- which(
    upper.tri(res_cor_out_full),
    arr.ind = TRUE
  )

  # Sort the order so it becomes correlations between x1 and other variables, then x2 with other variables
  upper_indices <- upper_indices[
    order(upper_indices[, "row"],
          upper_indices[, "col"]),
    ,
    drop = FALSE
  ]

  # Extract correlations
  res_cor <- res_cor_out_full[upper_indices]

  # Define variable name for correlations
  names(res_cor) <- paste0(
    prefix,
    rownames(res_cor_out_full)[upper_indices[, "row"]],
    "_",
    colnames(res_cor_out_full)[upper_indices[, "col"]]
  )

  # Convert to a one-row data frame
  res_cor_df <- as.data.frame(t(res_cor))

  return (res_cor_df)
}
