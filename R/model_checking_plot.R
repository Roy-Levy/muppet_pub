#' @title MUPPET Model Checking Plotting Functions
#'
#' @description Internal functions for generating posterior predictive model checking plots.
#'
#' @name model.check.plots
NULL

#' @param var.name Name of the statistic to be plotted.
#' @param summary.simulate.data Data frame containing univariate statistics calculated from posterior predicted data.
#' @param summary.input.data Data frame containing univariate statistics calculated from the observed data.
#' @param bins Number of bins in the histogram.
#'
#' @return A ggplot object.
#' @describeIn model.check.plots Generate a histogram of posterior predictive values for a univariate statistic,
#' with the realized value indicated by a dashed vertical line.
#'
model.checking.plot <- function(
    var.name,
    summary.simulate.data,
    summary.input.data,
    bins = 30
) {
    ggplot(
    summary.simulate.data,
    aes(x = .data[[var.name]])
    ) +
    geom_histogram(
      bins = bins,
      color = "black",
      fill = "grey80"
    ) +
    geom_vline(
      xintercept = summary.input.data[[var.name]][1],
      linetype = "dashed",
      linewidth = 1
    ) +
    labs(
      title = paste("Posterior predictive model checking of ", var.name),
      x = var.name,
      y = "Frequency"
    ) +
    theme_classic()
}


#' @param cor.simulate.long Data frame containing pairwise correlations calculated from posterior predicted data.
#' @param cor.input.long Data frame containing pairwise correlations calculated from the observed data.
#' @param bins Number of bins in the histogram.
#'
#' @describeIn model.check.plots Generate histograms of posterior predictive values for all pairwise correlations
#' in a single faceted plot, with the realized values indicated by dashed vertical lines.
#'
model.checking.cor.plot <- function(
    cor.simulate.long,
    cor.input.long,
    bins = 30
) {
  ggplot(
    data = cor.simulate.long,
    aes(x = value)
  ) +
  geom_histogram(
    bins = bins,
    color = "black",
    fill = "white"
  ) +
  geom_vline(
    data = cor.input.long,
    aes(xintercept = realized_value),
    linetype = "dashed"
  ) +
  facet_grid(
    var2 ~ var1,
    drop = FALSE
  ) +
  labs(
    title = "Posterior predictive model checking of pairwise correlations",
    x = "Correlation",
    y = "Frequency"
  ) +
  guides(
    x = guide_axis(check.overlap = TRUE) # prevent overlapping of ticks labels
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.1, 0.1))
  ) +
  theme_bw()
}


#' @param res.cor.long Data frame containing pairwise residual correlations calculated from posterior predicted data
#' and the observed data.
#'
#' @describeIn model.check.plots Generate scatter plots of posterior predictive residual correlations
#'  against realized residual correlations for all variable pairs in a single faceted plot.
#'
model.checking.res.cor.plot <- function(
    res.cor.long
) {
  lims <- range(
    c(res.cor.long$realized, res.cor.long$post_pred),
    na.rm = TRUE
  )

  ggplot(res.cor.long,
         aes(x = realized, y = post_pred)
         ) +
  geom_point(alpha = 0.5) +
  geom_abline(
    data = res.cor.long %>% distinct(var1, var2), # do not plot reference line in empty cell
    aes(intercept = 0, slope = 1),
    linetype = "dashed"
  ) +
  facet_grid(
    var2 ~ var1,
    drop = FALSE
  ) +
  coord_equal(
    xlim = lims,
    ylim = lims
  ) +
  guides(
    x = guide_axis(check.overlap = TRUE) # prevent overlapping of ticks labels
  ) +
  labs(
    title = "Posterior predictive model checking of residual correlations",
    x = "Realized residual correlation",
    y = "Posterior predictive residual correlation"
  ) +
  theme_bw()
}


#' @param fit.long Data frame containing fit indices calculated from posterior predicted data
#' and the observed data.
#' @param fit.name Name of the fit index to be plotted.
#'
#' @describeIn model.check.plots Generate a scatter plot of posterior predictive values against realized values
#' for a fit index.
#'
model.checking.fit.plot <- function(
    fit.long,
    fit.name
) {

  fit_labels <- c(
    Chisq = "chi-square"
  )

  fit_label <- if (fit.name %in% names(fit_labels)) {
    fit_labels[[fit.name]]
  } else {
    fit.name
  }

  plot_data <- fit.long %>%
  filter(fit == fit.name)

  lims <- range(
    c(plot_data$realized, plot_data$post_pred),
    na.rm = TRUE
  )

  ggplot(plot_data,
    aes(x = realized, y = post_pred)
  ) +
  geom_point(alpha = 0.5) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed"
  ) +
  coord_equal(
    xlim = lims,
    ylim = lims
  ) +
  labs(
    title = paste(
      "Posterior predictive model checking of",
      fit_label
    ),
    x = paste("Realized", fit_label),
    y = paste("Posterior predictive", fit_label)
  ) +
  theme_classic()
}


#' @param ppp_data Data frame containing posterior predictive p-values.
#' @param cor_type Type of correlations to be plotted (pairwise correlations or pairwise residual correlations).
#' @param cor_map Data frame mapping correlation statistic names to the names of the corresponding variable pairs.
#' @param var_order Order of variables in the observed data.
#'
#' @describeIn model.check.plots Generate a heat map of posterior predictive p-values for pairwise correlations
#' or pairwise residual correlations.
#'
ppp.plot <- function(
    ppp_data,
    cor_type,
    cor_map,
    var_order
){

  prefix <- if (cor_type == "res_cor") {
    "res_"
  } else if (cor_type == "cor") {
    ""
  }

  plot_title <- if (cor_type == "res_cor") {
    "Posterior predictive p-values of residual correlations"
  } else if (cor_type == "cor") {
    "Posterior predictive p-values of pairwise correlations"
  }

  ppp_cor <- ppp_data %>%
    select(statistic, ppp) %>%
    mutate(
      statistic = sub(paste0("^", prefix), "", statistic)
    ) %>%
    mutate(cor_join = tolower(statistic)) %>%
    inner_join(
      cor_map %>%
        mutate(cor_join = tolower(cor_name)),
      by = "cor_join"
    ) %>%
    select(-cor_join) %>%
    mutate(
      var1 = factor(var1, levels = var_order),
      var2 = factor(var2, levels = var_order)
    )

  grid <- expand.grid(
    var1 = var_order,
    var2 = var_order
  )

  # heat map
  p <- ggplot() +
    geom_tile(
      data = grid,
      aes(x = var1, y = var2),
      fill = NA,
      color = "grey80",
      linewidth = 0.4
    ) +
    geom_tile(
      data = ppp_cor,
      aes(x = var2, y = var1, fill = ppp),
      color = "grey80",
      linewidth = 0.4
    ) +
    scale_fill_gradient2(
      low = "#3B7EA1",
      mid = "#D9D9D9",
      high = "#D95F45",
      midpoint = 0.5,
      limits = c(0, 1),
      name = "ppp"
    ) +
    scale_x_discrete(
      limits = var_order,
      drop = FALSE, # show all boxes even if no data
      expand = c(0, 0) # no gap between matrix and axis
    ) +
    scale_y_discrete(
      limits = var_order,
      drop = FALSE,
      expand = c(0, 0)
    ) +
    coord_fixed() +
    labs(x = NULL,
         y = NULL,
         title = plot_title) +
    theme_classic()

  return (p)
}
