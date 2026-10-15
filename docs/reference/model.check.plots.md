# MUPPET Model Checking Plotting Functions

Internal functions for generating posterior predictive model checking
plots.

## Usage

``` r
model.checking.plot(
  var.name,
  summary.simulate.data,
  summary.input.data,
  bins = 30
)

model.checking.cor.plot(cor.simulate.long, cor.input.long, bins = 30)

model.checking.res.cor.plot(res.cor.long)

model.checking.fit.plot(fit.long, fit.name)

ppp.plot(ppp_data, cor_type, cor_map, var_order)
```

## Arguments

- var.name:

  Name of the statistic to be plotted.

- summary.simulate.data:

  Data frame containing univariate statistics calculated from posterior
  predicted data.

- summary.input.data:

  Data frame containing univariate statistics calculated from the
  observed data.

- bins:

  Number of bins in the histogram.

- cor.simulate.long:

  Data frame containing pairwise correlations calculated from posterior
  predicted data.

- cor.input.long:

  Data frame containing pairwise correlations calculated from the
  observed data.

- res.cor.long:

  Data frame containing pairwise residual correlations calculated from
  posterior predicted data and the observed data.

- fit.long:

  Data frame containing fit indices calculated from posterior predicted
  data and the observed data.

- fit.name:

  Name of the fit index to be plotted.

- ppp_data:

  Data frame containing posterior predictive p-values.

- cor_type:

  Type of correlations to be plotted (pairwise correlations or pairwise
  residual correlations).

- cor_map:

  Data frame mapping correlation statistic names to the names of the
  corresponding variable pairs.

- var_order:

  Order of variables in the observed data.

## Value

A ggplot object.

## Functions

- `model.checking.plot()`: Generate a histogram of posterior predictive
  values for a univariate statistic, with the realized value indicated
  by a dashed vertical line.

- `model.checking.cor.plot()`: Generate histograms of posterior
  predictive values for all pairwise correlations in a single faceted
  plot, with the realized values indicated by dashed vertical lines.

- `model.checking.res.cor.plot()`: Generate scatter plots of posterior
  predictive residual correlations against realized residual
  correlations for all variable pairs in a single faceted plot.

- `model.checking.fit.plot()`: Generate a scatter plot of posterior
  predictive values against realized values for a fit index.

- `ppp.plot()`: Generate a heat map of posterior predictive p-values for
  pairwise correlations or pairwise residual correlations.
