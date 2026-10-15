# MUPPET Analysis

Conduct modular MUPPET modeling.

## Usage

``` r
MUPPET.modular(
  fragments = NULL,
  n.chains = 2,
  n.warmup = 0,
  n.burnin = 500,
  n.iters.per.chain.after.warmup.and.burnin = 2500,
  n.thin = 1,
  n.estimation.batches = NULL,
  convergence.assessment = "unstandardized",
  save.summary.stats.from.MUPPET = TRUE,
  save.summary.stats.from.MUPPET.in.Word = FALSE,
  digits.to.round.for.Word = 2,
  save.summary.plots.from.MUPPET = "unstandardized",
  save.draws.from.MUPPET = TRUE,
  retain.iteration.files = FALSE,
  model.check = FALSE,
  n.model.check.iterations = NULL,
  save.post.pred.data = FALSE,
  save.summary.stats.from.model.check = TRUE,
  save.summary.plots.from.model.check = TRUE,
  retain.model.check.iteration.files = FALSE,
  return.R.object = FALSE
)
```

## Arguments

- fragments:

  A list with elements defining the model fragments. Each element in
  this list contains a set of specifications for the corresponding
  fragment. The 1st element in this list pertains to the 1st fragment,
  the 2nd element in this list pertains to the 2nd fragment, and so on.
  The order of the list also communicates the order each model fragment
  will be fit. Each element of this list is itself a list with
  specifications for the fragment. Each set of specifications can
  include:

  - name: Character, giving a name to the fragment. No default. If
    nothing is supplied, the fragment name will be the fragment number.

  - fragment.folder: Character, giving the name of the folder to be
    created storing the results for the fragment. If no value is
    supplied, will default to a name made up of the fragment number and
    fragment name.

  - model.syntax: Character, containing syntax for the Mplus MODEL
    statement.

  - variable.syntax: Character, containing syntax for the Mplus VARIABLE
    statement.

  - priors.syntax: Character, expression containing Mplus syntax for the
    prior specifications for the model in this fragment, corresponding
    to the Mplus MODELPRIORS statement. If nothing is specified,
    defaults to using Mplus's default priors.

  - conditioning: Number, or vector of numbers, indicating which
    fragments to condition on when fitting this fragment. Examples: a
    value of 0 indicates the fragment is unconditional; a value of 1
    indicates the fragment is conditional on the first fragment in the
    fragments list; a value of c(1, 2) indicates the fragment is
    conditional on fragment 1 and fragment 2 in the fragments list.

  - parameters.to.exclude.in.conditioning: A list, of length equal to
    the number of fragments that are conditioned on (i.e., the length of
    the collection supplied as the argument for `conditioning`). Each
    element of the list is a character, or character vector, naming the
    kinds of parameters from the associated antecedent fragment that
    should be ignored in fitting this fragment. Possible values are:

    - "none"

    - "loadings"

    - "(co)variances of observables"

    - "means and intercepts of observables"

    - "thresholds"

    - "structural coefficients"

    - "(co)variances of latent variables"

    - "means and intercepts of latent variables"

    A value of "none" indicates that all fitted parameters from the
    antecedent fragment should be conditioned on. For example, if
    `conditioning = c(1,2)`, and
    `parameters.to.exclude.in.conditioning = list("none", c("(co)variances of latent variables", "means and intercepts of latent variables"))`,
    then none of the parameters from fragment 1 will be ignored, while
    there will be parameters from fragment 2 that are ignored. These
    parameters include the (co)variances of latent variables and the
    means and intercepts of latent variables. If no value for
    `parameters.to.exclude.in.conditioning` is supplied, the default is
    "none" for all antecedent fragments.

  - to.fit: Logical, indicating if this fragment should be fit. If no
    value is supplied, defaults to TRUE. If value is FALSE, the fragment
    will not be fit. This can be useful when a fragment has previously
    been fit, and is now being used as an antecdent to another fragment.

  - estimating.lvs: Logical, indicating if this fragment pertains to
    estimating latent variable values. Default is FALSE. But using
    default may be finicky, so recommend specifying explicitly.
    `estimating.lvs = TRUE` should not be specified for unconditional
    fragment.

  - data: The dataset to use in fitting this fragment.

  - model.check.this.fragment: Logical, indicating if model checking
    should be conducted for this fragment. If no value is supplied,
    defaults to the value of `model.check`. Model checking is only
    available when `estimating.lvs = FALSE`.

- n.chains:

  The number of Markov chains to run in each unconditional fragment.
  Default is 2.

- n.warmup:

  The number of warmup iterations used by the Markov chains. This is
  applicable for stan and related software. Not needed, and need not be
  specified, when using Mplus. Default is 0.

- n.burnin:

  The number of iterations to discard as burn-in iterations for the
  Markov chains. See also the argument `n.thin`. Default is 500.

- n.iters.per.chain.after.warmup.and.burnin:

  The number of iterations desired, for each chain, after any warmup or
  burn-in period. See also the argument `n.thin`. Default is 2500.

- n.thin:

  The thinning parameter, for each Markov chain. Importantly, when using
  Mplus, this applies to all iterations. For example, suppose
  `n.burnin = 500`, `n.iters.per.chain.after.warmup.and.burnin = 2500`,
  and `n.thin = 2`. This effectively instructs Mplus to run 1000
  iterations as burn-in (i.e., 1000 iterations thinned by 2 is 500) and
  5000 iterations after burn-in (i.e., 5000 iterations thinned by 2 is
  2500). Default is 1.

- n.estimation.batches:

  The number of estimation batches when fitting fragments that are
  conditional on antecedent fragments. Default is NULL, in which case
  the function will divide the total number of post burn-in iterations
  by 1000, effectively doing batches of 1000 iterations at a time. Using
  more estimation batches leads to fewer iterations per estimation
  batch. This can be helpful for reducing computational burden and run
  time, and also offers future possibilities for checkpointing.

  If `n.estimation.batches` is specified, the total number of post
  burn-in iterations must be a multiple of `n.estimation.batches`. If
  left as NULL, the total number of post burn-in iterations must be a
  multiple of 1000 when exceeding 1000 iterations, otherwise, a single
  batch is used.

- convergence.assessment:

  Character, indicating if convergence assessment is desired for
  fragments that are not conditional on antecedents. Options are:
  "unstandardized" to only plot unstandardized parameters, "all" to plot
  standardized and unstandardized parameters, "none" for turning off
  convergence assessment. This argument is irrelevant for fragments that
  are conditional on antecedent fragments. Default is "unstandardized".

- save.summary.stats.from.MUPPET:

  Logical, indicating if numerical summaries of posterior distributions
  be saved. Default is TRUE.

- save.summary.stats.from.MUPPET.in.Word:

  Logical, indicating if slimmed down numerical summary statistics of
  the posterior distributions should be saved as Word documents. Default
  is FALSE.

- digits.to.round.for.Word:

  Numerical. Digits to round if saving results in Word documents. Only
  relevant if `save.summary.stats.from.MUPPET.in.Word = TRUE`. Default =
  2.

- save.summary.plots.from.MUPPET:

  Character, indicating if summary plots of the posterior should be
  saved. Options are: "unstandardized" to only plot unstandardized
  parameters, "all" to plot standardized and unstandardized parameters,
  "none" for turning off these plots. Default is "unstandardized".

- save.draws.from.MUPPET:

  Logical, indicating whether the draws from fitting the MUPPET model
  should be written out as an output file. Default is TRUE.

- retain.iteration.files:

  Logical. For conditional fragments, should the iteration-specific
  files be saved? Saving all these files may take up a lot of space. If
  set to FALSE, will still save the files for iteration 1, which may be
  useful for troubleshooting if needed. Default is FALSE.

- model.check:

  Logical, indicating if posterior predictive model checking should be
  conducted. Default is FALSE.

  In the current version, model checking is available for univariate
  statistics (mean, SD, minimum, and maximum), pairwise correlations,
  pairwise residual correlations, and global fit indices (Chi-square,
  CFI, TLI, RMSEA, and SRMR). Note that pairwise residual correlations
  and global fit indices are not available for models with categorical
  outcomes.

- n.model.check.iterations:

  Number of iterations (number of posterior predicted datasets) used in
  model checking. Default to use all iterations.

- save.post.pred.data:

  Logical, indicating if posterior predicted datasets from the model
  checking procedure should be saved. Default is FALSE.

- save.summary.stats.from.model.check:

  Logical, indicating if numerical summaries of model checking should be
  saved. Default is TRUE.

- save.summary.plots.from.model.check:

  Logical, indicating if plots of model checking should be saved.
  Default is TRUE.

- retain.model.check.iteration.files:

  Logical. Should the iteration-specific files be saved when conducting
  model checking? If set to FALSE, will still save the
  simulation-related files for iteration 1, which may be useful for
  troubleshooting if needed. Default is FALSE.

- return.R.object:

  Logical, indicating whether to return results as an R object. Default
  is FALSE. This was supported in earlier versions for a limited class
  of models. This is not current, and should not be set to TRUE. For
  future development.

## Value

A statement of how long it took to run. The requested output files will
be written out to the hard drive. Future development will include
returning objects to R.

## See also

For examples and tips, see [MUPPET Examples and
Tips](https://roy-levy.github.io/muppet_pub/articles/index.html)
