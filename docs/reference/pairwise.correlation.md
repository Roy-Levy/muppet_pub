# Process Pairwise Correlations

Internal function for calculating pairwise correlations. Optionally, may
return a mapping between correlation statistic names (e.g.,
"cor_var1_var2") and the corresponding variable names (e.g., "var1" and
"var2").

## Usage

``` r
pairwise.correlation(data, return_map = FALSE)
```

## Arguments

- data:

  Data frame used to calculate pairwise correlations.

- return_map:

  Logical, indicating whether to return a mapping of correlation
  statistic names to the names of the corresponding variable pairs.
  Default is FALSE.

## Value

If `return_map = FALSE`, a named numeric vector of pairwise
correlations.

If `return_map = TRUE`, a list containing:

- `cor_data`: A named numeric vector of pairwise correlations.

- `cor_map`: A data frame mapping correlation statistic names to the
  names of the corresponding variable pairs.
