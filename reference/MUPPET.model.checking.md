# MUPPET Model Checking

Conduct posterior predictive model checking. This is an internal
function called by `MUPPET.modular`, and is not intended to be called
directly by users.

## Usage

``` r
MUPPET.model.checking(
  fragments,
  which.fragment,
  n.estimation.batches,
  iters.per.batch,
  home.folder,
  n.model.check.iterations,
  save.post.pred.data,
  save.summary.stats.from.model.check,
  save.summary.plots.from.model.check,
  retain.model.check.iteration.files
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

- which.fragment:

  The fragment number for which model checking is conducted.

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

- iters.per.batch:

  Number of iterations per batch, calculated in `MUPPET.modular`.

- home.folder:

  Working directory obtained from
  [`getwd()`](https://rdrr.io/r/base/getwd.html) in `MUPPET.modular`.

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

## Value

Model checking progress and completion messages. Model checking output
files will be written out to the hard drive.
