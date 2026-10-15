# Process Residual Correlations

Internal function for reordering pairwise residual correlations
according to the variable order in the observed data.

## Usage

``` r
residual.cor.reorder(res_cor_out, variable_list, prefix)
```

## Arguments

- res_cor_out:

  Matrix of residual correlations extracted from an Mplus output file.

- variable_list:

  Character string containing the variable names in the observed data.

- prefix:

  Character string used as prefix of residual correlation names.

## Value

A one-row data frame of pairwise residual correlations.
