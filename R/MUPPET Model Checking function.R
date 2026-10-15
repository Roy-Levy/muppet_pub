
#' @title MUPPET Model Checking
#'
#' @description  Conduct posterior predictive model checking. This is an internal function called by \code{MUPPET.modular}, and is not intended to be called directly by users.
#'
#' @inheritParams MUPPET.modular

#' @param which.fragment The fragment number for which model checking is conducted.
#' @param iters.per.batch Number of iterations per batch, calculated in \code{MUPPET.modular}.
#' @param home.folder Working directory obtained from \code{getwd()} in \code{MUPPET.modular}.
#'
#' @returns Model checking progress and completion messages. Model checking output files
#' will be written out to the hard drive.
#'
MUPPET.model.checking <- function(
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
){

  old.wd <- getwd()
  # Restore the working directory when the model checking function exits
  on.exit(setwd(old.wd), add = TRUE)

  setwd(home.folder)

  ##
  current.time <- Sys.time()
  print(paste0("Model checking for Fragment ", which.fragment, " initiated ", round(current.time, units="secs")))

  sink("R Session Info.out", append = TRUE)
  print(paste0("Model checking for Fragment ", which.fragment, " initiated ", round(current.time, units="secs")))
  sink()

  # Set up model checking folder
  model.checking.folder <- paste0(
    fragments[[which.fragment]]$fragment.folder,
    "/Model Checking/")
  if(!dir.exists(model.checking.folder)) dir.create(model.checking.folder)

  ######################################
  # Check if simulated data already exist
  post.pred.data.folder <- file.path(
    model.checking.folder,
    "Posterior Predicted Data"
  )

  # Check if simulated data are saved in Posterior Predicted Data folder
  simulation.files.post.pred.data.folder <- character(0)
  if (dir.exists(post.pred.data.folder)) {
    simulation.files.post.pred.data.folder <- list.files(
      post.pred.data.folder,
      pattern = "\\.iter\\.[0-9]+\\.simulation\\.dat$"
    )
  }

  # Check if simulated data are saved in Model Checking folder
  # This may happen if the user retains all iteration-specific files
  simulation.files.model.checking.folder <- list.files(
    model.checking.folder,
    pattern = "\\.iter\\.[0-9]+\\.simulation\\.dat$"
  )

  # Use existing simulated data when available; otherwise, regenerate them
  if (length(simulation.files.post.pred.data.folder) >0){
    need.generate.post.pred.data <- FALSE
    simulation.files.folder <- post.pred.data.folder
  }
  # Because first iteration files are always retained, so need to check if there are more than 1 retained
  else if (length(simulation.files.model.checking.folder) >1) {
    need.generate.post.pred.data <- FALSE
    simulation.files.folder <- model.checking.folder
  } else {
    need.generate.post.pred.data <- TRUE
    simulation.files.folder <- model.checking.folder
  }


  ###############
  # Variable list is needed even if simulated datasets already exist

  # Read in the Mplus output object from the current fragment

  # Define the filename if the fragment is NOT conditional on other fragments
  if(sum(fragments[[which.fragment]]$conditioning) == 0){
    file.name <- paste0(
      fragments[[which.fragment]]$name, ".Mplus output object.rds"
    )
  }

  # Define the filename if the fragment is conditional on other fragments
  if(sum(fragments[[which.fragment]]$conditioning) != 0){
    file.name <- paste0(
      fragments[[which.fragment]]$name, ".shell Mplus output object.rds"
    )
  }

  current.fragment.Mplus.output <- readRDS(file.name)

  # Extract variable list
  variable_list <- current.fragment.Mplus.output$input$variable$names
  # Extract categorical variable list
  categorical_variable_syntax <- current.fragment.Mplus.output$input$variable$categorical

  ################################
  ################################
  # A list of univariate statistics to compute
  stats <- list(
    mean = ~ mean(.x, na.rm = TRUE),
    sd = ~ sd(.x, na.rm = TRUE),
    min = ~ min(.x, na.rm = TRUE),
    max = ~ max(.x, na.rm = TRUE)
  )

  # * * * * Read in the input data and calculate statistics ------------
  input.data.file <- current.fragment.Mplus.output$input$data$file

  # First look in the home folder
  input.data.file.path <- file.path(home.folder, input.data.file)

  # If not found, look in the fragment folder
  if (!file.exists(input.data.file.path)) {
    input.data.file.path <- file.path(
      fragments[[which.fragment]]$fragment.folder,
      input.data.file
    )
  }

  input.data <- read.table(
    input.data.file.path,
    header = FALSE,
    na.strings = "."
  )

  colnames(input.data) <- unlist(strsplit(variable_list, "\\s+"))

  # Summary statistics for input data
  summary.input.data.univariate <- input.data %>%
    summarise(
      across(everything(), stats)
    )

  # Pairwise correlations
  input_cor_result <- pairwise.correlation(
    input.data,
    return_map = TRUE
  )

  input_cor_df <- as.data.frame(as.list(input_cor_result$cor_data))
  input_cor_df <- input_cor_df %>%
    # This is to drop the statistic if it is not applicable to the input data
    # In IRT example, correlations between some variables may be NA if the items are not administered on the same persons
    select(where(~ !all(is.na(.x))))
  # Map the correlation name to the names of two variables
  cor_map <- input_cor_result$cor_map

  # Univariate statistics and pairwise correlations
  summary.input.data <-
    bind_cols(
      summary.input.data.univariate,
      input_cor_df
    )

  #######################################
  # If need to generate posterior predicted data
  if (need.generate.post.pred.data){

    # Extract number of observations
    nobs <- current.fragment.Mplus.output$summaries$Observations

    # For categorical variables
    if (!is.null(categorical_variable_syntax)){
      # Fix the residual variance of the underlying latent response variables to 1 (theta parameterization)
      categorical_variables <- strsplit(trimws(categorical_variable_syntax), "\\s+")[[1]] # Remove leading/trailing whitespace, and split by space
      categorical_syntax <- paste0(categorical_variables, "@1;", collapse = "\n")
    } else categorical_syntax <- NULL

    ####
    # If the fragment is NOT conditional on other fragments
    # Extract the model syntax where the parameters are fixed
    # Keep only lines containing "@"
    if(sum(fragments[[which.fragment]]$conditioning) == 0){
      fixed_lines <- grep("@", current.fragment.Mplus.output$input$model, value = TRUE)
      current.fragment.fixed.syntax <- paste(fixed_lines, collapse = "\n")
    }

    ####
    # Read in draws of the current fragment
    current.fragment.draws <- read_csv(
      file = paste0(
        fragments[[which.fragment]]$name, " draws.csv")
    )
    # Only keep unstandardized values by discarding columns contain STD
    # Also discard columns contain R-SQ
    current.fragment.draws.unstd <- current.fragment.draws %>%
      select(-contains("STD"), -contains("R-SQ"))


    ##################
    # Set working directory to the model checking folder
    setwd(model.checking.folder)

    # * * * Clean up existing parallel environment ------
    unregister_dopar <- function() {
      env <- foreach:::.foreachGlobals
      rm(list = ls(name = env), pos = env)
    }

    # * * * Setup parallel environment ------

    # * * * Detect and set up the number of cores
    available_cores <- future::availableCores() # This can also be used in HPC (available cores to use)

    # If using HPC then use all cores
    if (nzchar(Sys.getenv("SLURM_CPUS_PER_TASK"))) {
      numCores <- available_cores
    } else {
      numCores <- max(1L, available_cores - 1L)
    }

    registerDoParallel(numCores)

    # * * * Loop over estimation batches ------
    for(which.estimation.batch in 1:n.estimation.batches){

      # Define the start and end iteration for each batch
      start.iter = (which.estimation.batch-1)*iters.per.batch + 1
      end.iter = (which.estimation.batch-1)*iters.per.batch + iters.per.batch

      # Loop through each iteration
      foreach(
        which.iter = start.iter:end.iter,
        .packages = c("dplyr", "MplusAutomation"),
        .combine = rbind
      ) %dopar%
        {
          which.iter = as.integer(which.iter)

          # * * *  * ******************************************************************* ----
          # * * * * Define the Mplus MODEL portion of the input file -----

          # * * * * * Initiate an empty object. Will append text ------
          temp <- NULL

          # * * * * * A list of parameters -------
          # 1) loadings
          # 2) (co)variances of observables
          # 3) means and intercepts of observables
          # 4) thresholds
          # 5) structural coefficients
          # 6) (co)variances of latent variables
          # 7) means and intercepts of latent variables


          # * * * * * Add code for loadings ('lambda' in Mplus)-------

          # Extract the parameter specification from Mplus as a data frame
          parameter.specification <- as.data.frame(
            current.fragment.Mplus.output$tech1$parameterSpecification$lambda
          )

          # Proceed if there are any fitted parameters
          if(any(parameter.specification > 0, na.rm = TRUE)){

            # Go through each variable (row in the data frame)
            for(which.variable in 1:nrow(parameter.specification)){
              # Go through each other variable (column in the data frame)
              for(which.other.variable in 1:ncol(parameter.specification)){
                # Proceed if it's a a fitted parameter
                if(!is.na(parameter.specification[which.variable, which.other.variable]) && parameter.specification[which.variable, which.other.variable] >0){
                  # Define the parameter name for the Mplus input code
                  parameter.name <- paste0(
                    colnames(parameter.specification)[which.other.variable],
                    " by ",
                    rownames(parameter.specification)[which.variable]
                  )
                  # Define the parameter name in the matrix of draws
                  parameter.name.in.draws <- paste0(
                    colnames(parameter.specification)[which.other.variable],
                    ".by.",
                    rownames(parameter.specification)[which.variable]
                  )

                  # Define the value of the entity
                  matched_cols <- select(
                    current.fragment.draws.unstd,
                    matches(parameter_string_pattern(parameter.name.in.draws))
                  )

                  if (ncol(matched_cols) != 1) {
                    stop(
                      "Expect exactly one column matching '",
                      parameter.name,
                      "', but found ",
                      ncol(matched_cols),
                      "."
                    )
                  }

                  value.to.paste <- matched_cols[which.iter, 1]

                  # Define the line to append to syntax
                  this.line <- paste0(
                    parameter.name,
                    " @",
                    value.to.paste,
                    ";"
                  )

                  # Append the line
                  temp <- paste0(
                    temp, "\n",
                    this.line
                  )

                } # closes if this is a fitted parameter
              } # closes loop over "other" variable (columns)
            } # closes loop over each variable with fitted parameters (rows)
          } # closes if there are fitted loadings


          # * * * * * Add code for (co)variances of observables ('theta' in Mplus)  -------

          # Extract the parameter specification from Mplus as a data frame
          parameter.specification <- as.data.frame(
            current.fragment.Mplus.output$tech1$parameterSpecification$theta
          )

          # Proceed if there are any fitted parameters
          if(any(parameter.specification > 0, na.rm = TRUE)){
            # Go through each variable (row in the data frame)
            for(which.variable in 1:nrow(parameter.specification)){
              # Go through each other variable (column in the data frame)
              for(which.other.variable in 1:ncol(parameter.specification)){

                # Proceed if it's a a fitted parameter
                if(!is.na(parameter.specification[which.variable, which.other.variable]) && parameter.specification[which.variable, which.other.variable] >0){

                  # If it's a variance
                  if(rownames(parameter.specification)[which.variable] ==
                     colnames(parameter.specification)[which.other.variable]){

                    # Define the parameter name for the Mplus input code
                    parameter.name <- paste0(
                      rownames(parameter.specification)[which.variable]
                    )

                    # Define the parameter name in the matrix of draws
                    parameter.name.in.draws <- parameter.name

                    # Define the value of the entity
                    matched_cols <- select(
                      current.fragment.draws.unstd,
                      matches(parameter_string_pattern(parameter.name.in.draws))
                    )

                    if (ncol(matched_cols) != 1) {
                      stop(
                        "Expect exactly one column matching '",
                        parameter.name,
                        "', but found ",
                        ncol(matched_cols),
                        "."
                      )
                    }

                    value.to.paste <- matched_cols[which.iter, 1]

                    # Define the line to append to syntax
                    this.line <- paste0(
                      parameter.name,
                      " @",
                      value.to.paste,
                      ";"
                    )

                    # Append the line
                    temp <- paste0(
                      temp, "\n",
                      this.line
                    )

                  } # closes if it's a variance

                  # If it's a covariance
                  if(rownames(parameter.specification)[which.variable] !=
                     colnames(parameter.specification)[which.other.variable]){

                    # Define the parameter name for the Mplus input code
                    parameter.name <- paste0(
                      colnames(parameter.specification)[which.variable],
                      " with ",
                      rownames(parameter.specification)[which.other.variable]
                    )

                    # Define the parameter name in the matrix of draws
                    parameter.name.in.draws <- paste0(
                      rownames(parameter.specification)[which.variable],
                      ".with.",
                      colnames(parameter.specification)[which.other.variable]
                    )

                    # check if that parameter name is in the matrix of draws
                    # If it ISN'T, reverse the ordering of the variables in the name
                    if(ncol(
                      select(
                        current.fragment.draws.unstd,
                        matches(parameter_string_pattern(parameter.name.in.draws))
                      )
                    ) == 0){
                      parameter.name.in.draws <- paste0(
                        colnames(parameter.specification)[which.other.variable],
                        ".with.",
                        rownames(parameter.specification)[which.variable]
                      )
                    } # closes if the parameter name is not in the matrix of draws

                    # Define the value of the entity
                    matched_cols <- select(
                      current.fragment.draws.unstd,
                      matches(parameter_string_pattern(parameter.name.in.draws))
                    )

                    if (ncol(matched_cols) != 1) {
                      stop(
                        "Expect exactly one column matching '",
                        parameter.name,
                        "', but found ",
                        ncol(matched_cols),
                        "."
                      )
                    }

                    value.to.paste <- matched_cols[which.iter, 1]

                    # Define the line to append to syntax
                    this.line <- paste0(
                      parameter.name,
                      " @",
                      value.to.paste,
                      ";"
                    )

                    # Append the line
                    temp <- paste0(
                      temp, "\n",
                      this.line
                    )

                  } # closes if it's a covariance
                } # closes if this is a fitted parameter
              } # closes loop over "other" variable (columns)
            } # closes loop over each variable with fitted parameters (rows)
          } # closes if there are fitted (co)variances of observables


          # * * * * * Add code for observable means and intercepts ('nu' in Mplus)-------

          # Extract the parameter specification from Mplus as a data frame
          parameter.specification <- as.data.frame(
            current.fragment.Mplus.output$tech1$parameterSpecification$nu
          )

          # Proceed if there are any fitted parameters
          if(any(parameter.specification > 0, na.rm = TRUE)){

            # Loop over observables
            for(which.variable in 1:ncol(parameter.specification)){

              # proceed if this parameter is fitted (i.e., the value of the parameter specification table is > 0)
              if(!is.na(parameter.specification[1, which.variable]) && parameter.specification[1, which.variable] >0){

                # Define the parameter name for the Mplus input code
                parameter.name <- colnames(parameter.specification)[which.variable]

                # Define the parameter name in the matrix of draws
                parameter.name.in.draws <- paste0(
                  "MEAN.",
                  colnames(parameter.specification)[which.variable]
                )

                # Define the value of the entity
                matched_cols <- select(
                  current.fragment.draws.unstd,
                  matches(parameter_string_pattern(parameter.name.in.draws))
                )

                if (ncol(matched_cols) != 1) {
                  stop(
                    "Expect exactly one column matching '",
                    parameter.name,
                    "', but found ",
                    ncol(matched_cols),
                    "."
                  )
                }

                value.to.paste <- matched_cols[which.iter, 1]

                this.line <- paste0(
                  "[", parameter.name,"@",
                  value.to.paste,
                  "];"
                )

                # Append the line
                temp <- paste0(
                  temp, "\n",
                  this.line
                )

              } # closes if this parameter is fitted
            } # closes loop over variable
          } # closes if there are fitted observable means and intercepts


          # * * * * * Add code for thresholds ('tau' in Mplus) -------

          # Extract the parameter specification from Mplus as a data frame
          parameter.specification <- as.data.frame(
            current.fragment.Mplus.output$tech1$parameterSpecification$tau
          )

          # Proceed if there are any fitted parameters
          if(any(parameter.specification > 0, na.rm = TRUE)){

            # Loop over observables
            # Note the name of the observables comes from the loading matrix
            # This is because the names in the parameter.specification table can be too long, and ends get cut off
            for(which.observable in 1:nrow(current.fragment.Mplus.output$tech1$parameterSpecification$lambda)){

              # proceed if this parameter is fitted (i.e., the value of the parameter specification table is > 0)
              if(!is.na(parameter.specification[1, which.observable]) && parameter.specification[1, which.observable] >0){

                # Define the parameter name for the Mplus input code
                parameter.name <- rownames(current.fragment.Mplus.output$tech1$parameterSpecification$lambda)[which.observable]

                # Define the parameter name in the matrix of draws
                parameter.name.in.draws <- paste0(
                  parameter.name, "$"
                )

                # select the thresholds for this observable
                draws.for.thresholds <- select(
                  current.fragment.draws.unstd,
                  contains(paste0(".", parameter.name.in.draws))
                )

                n.thresholds.this.observable <- ncol(draws.for.thresholds)

                # Loop over thresholds
                for(which.threshold in 1:n.thresholds.this.observable){

                  value.to.paste <- draws.for.thresholds[which.iter, which.threshold]

                  this.line <- paste0(
                    "[", parameter.name, "$", which.threshold, "@",
                    value.to.paste,
                    "];"
                  )

                  # Append the line
                  temp <- paste0(
                    temp, "\n",
                    this.line
                  )

                } # closes loop over thresholds
              } # closes if this parameter is fitted
            } # closes loop over observables
          } # closes if there are fitted thresholds


          # * * * * * Add code for structural coefficients ('beta' in Mplus)-------

          # Extract the parameter specification from Mplus as a data frame
          parameter.specification <- as.data.frame(
            current.fragment.Mplus.output$tech1$parameterSpecification$beta
          )

          # Proceed if there are any fitted parameters
          if(any(parameter.specification > 0, na.rm = TRUE)){

            # Go through each variable (row in the data frame)
            for(which.variable in 1:nrow(parameter.specification)){
              # Go through each other variable (column in the data frame)
              for(which.other.variable in 1:ncol(parameter.specification)){

                # Proceed if it's a fitted parameter
                if(!is.na(parameter.specification[which.variable, which.other.variable]) && parameter.specification[which.variable, which.other.variable] >0){

                  # Define the parameter name for the Mplus input code
                  parameter.name <- paste0(
                    rownames(parameter.specification)[which.variable],
                    " on ",
                    colnames(parameter.specification)[which.other.variable]
                  )

                  # Define the parameter name in the matrix of draws
                  parameter.name.in.draws <- paste0(
                    rownames(parameter.specification)[which.variable],
                    ".on.",
                    colnames(parameter.specification)[which.other.variable]
                  )

                  # Define the value of the entity
                  matched_cols <- select(
                    current.fragment.draws.unstd,
                    matches(parameter_string_pattern(parameter.name.in.draws))
                  )

                  if (ncol(matched_cols) != 1) {
                    stop(
                      "Expect exactly one column matching '",
                      parameter.name,
                      "', but found ",
                      ncol(matched_cols),
                      "."
                    )
                  }

                  value.to.paste <- matched_cols[which.iter, 1]

                  # Define the line to append to syntax
                  this.line <- paste0(
                    parameter.name,
                    " @",
                    value.to.paste,
                    ";"
                  )

                  # Append the line
                  temp <- paste0(
                    temp, "\n",
                    this.line
                  )

                } # closes if this is a fitted parameter
              } # closes loop over "other" variable (columns)
            } # closes loop over each variable with fitted parameters (rows)
          } # closes if there are fitted structural coefficients



          # * * * * * Add code for (co)variances of latent variables ('psi' in Mplus)  -------

          # Extract the parameter specification from Mplus as a data frame
          parameter.specification <- as.data.frame(
            current.fragment.Mplus.output$tech1$parameterSpecification$psi
          )

          # Proceed if there are any fitted parameters
          if(any(parameter.specification > 0, na.rm = TRUE)){

            # Go through each variable (row in the data frame)
            for(which.variable in 1:nrow(parameter.specification)){
              # Go through each other variable (column in the data frame)
              for(which.other.variable in 1:ncol(parameter.specification)){

                # Proceed if it's a a fitted parameter
                if(!is.na(parameter.specification[which.variable, which.other.variable]) && parameter.specification[which.variable, which.other.variable] >0){

                  # If it's a variance
                  if(rownames(parameter.specification)[which.variable] ==
                     colnames(parameter.specification)[which.other.variable]){

                    # Define the parameter name for the Mplus input code
                    parameter.name <- paste0(
                      rownames(parameter.specification)[which.variable]
                    )

                    # Define the parameter name in the matrix of draws
                    parameter.name.in.draws <- parameter.name

                    # Define the value of the entity
                    matched_cols <- select(
                      current.fragment.draws.unstd,
                      matches(parameter_string_pattern(parameter.name.in.draws))
                    )

                    if (ncol(matched_cols) != 1) {
                      stop(
                        "Expect exactly one column matching '",
                        parameter.name,
                        "', but found ",
                        ncol(matched_cols),
                        "."
                      )
                    }

                    value.to.paste <- matched_cols[which.iter, 1]

                    # Define the line to append to syntax
                    this.line <- paste0(
                      parameter.name,
                      " @",
                      value.to.paste,
                      ";"
                    )

                    # Append the line
                    temp <- paste0(
                      temp, "\n",
                      this.line
                    )

                  } # closes if it's a variance

                  # If it's a covariance
                  if(rownames(parameter.specification)[which.variable] !=
                     colnames(parameter.specification)[which.other.variable]){

                    # Define the parameter name for the Mplus input code
                    parameter.name <- paste0(
                      colnames(parameter.specification)[which.variable],
                      " with ",
                      rownames(parameter.specification)[which.other.variable]
                    )

                    # Define the parameter name in the matrix of draws
                    parameter.name.in.draws <- paste0(
                      rownames(parameter.specification)[which.variable],
                      ".with.",
                      colnames(parameter.specification)[which.other.variable]
                    )

                    # check if that parameter name is in the matrix of draws
                    # If it ISN'T, reverse the ordering of the variables in the name
                    if(ncol(
                      select(
                        current.fragment.draws.unstd,
                        matches(parameter_string_pattern(parameter.name.in.draws))
                      )
                    ) == 0){
                      parameter.name.in.draws <- paste0(
                        colnames(parameter.specification)[which.other.variable],
                        ".with.",
                        rownames(parameter.specification)[which.variable]
                      )
                    }

                    # Define the value of the entity
                    matched_cols <- select(
                      current.fragment.draws.unstd,
                      matches(parameter_string_pattern(parameter.name.in.draws))
                    )

                    if (ncol(matched_cols) != 1) {
                      stop(
                        "Expect exactly one column matching '",
                        parameter.name,
                        "', but found ",
                        ncol(matched_cols),
                        "."
                      )
                    } # closes if the parameter name is not in the matrix of draws

                    value.to.paste <- matched_cols[which.iter, 1]

                    # Define the line to append to syntax
                    this.line <- paste0(
                      parameter.name,
                      " @",
                      value.to.paste,
                      ";"
                    )

                    # Append the line
                    temp <- paste0(
                      temp, "\n",
                      this.line
                    )

                  } # closes if it's a covariance
                } # closes if this is a fitted parameter
              } # closes loop over "other" variable (columns)
            } # closes loop over each variable with fitted parameters (rows)
          } # closes if there are fitted (co)variances of latent variables


          # * * * * * Add code for latent means and intercepts ('alpha' in Mplus) -------

          # Extract the parameter specfication from Mplus as a data frame
          parameter.specification <- as.data.frame(
            current.fragment.Mplus.output$tech1$parameterSpecification$alpha
          )

          # Proceed if there are any fitted parameters
          if(any(parameter.specification > 0, na.rm = TRUE)){

            # Loop over variables
            for(which.variable in 1:ncol(parameter.specification)){

              # Proceed if this parameter is fitted (i.e., the value of the parameter specification table is > 0)
              if(!is.na(parameter.specification[1, which.variable]) && parameter.specification[1, which.variable] >0){

                # Define the parameter name for the Mplus input code
                parameter.name <- colnames(parameter.specification)[which.variable]

                # Define the parameter name in the matrix of draws
                parameter.name.in.draws <- paste0(
                  "MEAN.", parameter.name
                )

                # Define the value of the entity
                matched_cols <- select(
                  current.fragment.draws.unstd,
                  matches(parameter_string_pattern(parameter.name.in.draws))
                )

                if (ncol(matched_cols) != 1) {
                  stop(
                    "Expect exactly one column matching '",
                    parameter.name,
                    "', but found ",
                    ncol(matched_cols),
                    "."
                  )
                }

                value.to.paste <- matched_cols[which.iter, 1]

                this.line <- paste0(
                  "[", parameter.name,"@",
                  value.to.paste,
                  "];"
                )

                # Append the line
                temp <- paste0(
                  temp, "\n",
                  this.line
                )

              } # closes if this parameter is fitted (i.e., the value of the parameter specification table is > 0)
            } # closes loop over variable
          } # closes if there are fitted latent means and intercepts


          # * * *  * ******************************************************************* ----
          # # * * * * Define the syntax for the MODEL POPULATION portion of Mplus sytnax -------

          # If the fragment is conditional on other fragments
          # The fixed portion is different by iterations, so need to extract the output object for each iteration
          # Extract the model syntax where the parameters are fixed
          # Keep only lines containing "@"
          if(sum(fragments[[which.fragment]]$conditioning) != 0){
            file.name <- paste0(
              fragments[[which.fragment]]$name, ".iter.", which.iter, ".Mplus output object.rds")

            current.fragment.Mplus.output <- readRDS(file.path(home.folder, file.name))

            fixed_lines <- grep("@", current.fragment.Mplus.output$input$model, value = TRUE)
            current.fragment.fixed.syntax <- paste(fixed_lines, collapse = "\n")
          }

          # first portion: parameters fixed at the values from the posterior draw
          # second portion: parameters fixed as in the user input (or fixed based on draws from antecedent fragments)
          # third portion: add categorical syntax for categorical dependent variables
          Mplus.MODELPOPULATION.syntax <- paste0(
            temp,
            "\n",
            "\n",
            current.fragment.fixed.syntax,
            "\n",
            categorical_syntax
          )

          # Construct data frame for categorical variables and the number of thresholds
          matches_pattern <- gregexpr(
            "\\b[A-Za-z][A-Za-z0-9_]*\\$[0-9]+",
            Mplus.MODELPOPULATION.syntax,
            perl = TRUE
          )
          # Extract syntax for thresholds
          threshold_terms <- regmatches(Mplus.MODELPOPULATION.syntax, matches_pattern)[[1]]

          item <- sub("\\$[0-9]+$", "", threshold_terms)
          threshold_order <- as.integer(
            sub("^.*\\$", "", threshold_terms)
          )

          # Find the number of thresholds for each item
          threshold_counts <- tapply(
            threshold_order,
            factor(item, levels = unique(item)),
            max
          )

          # Put in a data frame
          number_thresholds <- data.frame(
            item = names(threshold_counts),
            num_thresholds = as.integer(threshold_counts),
            row.names = NULL
          )

          # * * * * Define the syntax for the MONTECARLO portion of Mplus sytnax -----

          # Syntax for categorical variables
          # Note that the estimates are in probit scale, so need to specify the categorical variables be generated using probit link
          if (nrow(number_thresholds)!= 0) {
            generate_terms <- paste0(
              number_thresholds$item,
              "(",
              number_thresholds$num_threshold,
              " p)"
            )

            generate_syntax <- paste0(
              "GENERATE = ",
              paste(generate_terms, collapse = " "),
              ";\n"
            )
          } else {
            generate_syntax <- NULL
          }

          Mplus.MONTECARLO.syntax <-
            paste0("NAMES ARE ", variable_list, ";\n",
                   generate_syntax,
                   paste0("NOBS = ", nobs, ";\n"),
                   "NREP = 1;\n",
                   paste0("SEED = ", which.iter,";\n"), # use a different seed for each iteration
                   "REPSAVE = ALL;\n",
                   paste0("SAVE = ", fragments[[which.fragment]]$name, ".iter.", which.iter, ".simulation.dat;\n")
            )

          # * * * * Define the title for Mplus model -----
          Mplus.TITLE <- paste0(fragments[[which.fragment]]$name, " simulation;")

          # Define the Mplus object
          fragment.simulation.mplus.object <- mplusObject(
            TITLE = Mplus.TITLE,
            MONTECARLO = Mplus.MONTECARLO.syntax,
            MODELPOPULATION = Mplus.MODELPOPULATION.syntax
          )

          # * * * * Fit the model in Mplus ------------
          fitted.model.mplus <- mplusModeler(
            object = fragment.simulation.mplus.object,
            modelout = paste0(fragments[[which.fragment]]$name, ".iter.", which.iter, ".simulation.inp"),
            run = 1L
          )


          ##########################
          # Model checking for fit statistics
          # Need to do a separate Mplus run, fixing the parameter values at each draw
          # Fit on input data, as well as iteration-specific simulated data

          # ONLY for the case when all variables are continuous
          if (is.null(categorical_variable_syntax)){

            # Define simulated dataset for each iteration
            simulate.data <- read.table(
              paste0(fragments[[which.fragment]]$name, ".iter.", which.iter, ".simulation.dat"),
              header = FALSE
            )
            names(simulate.data) <- unlist(strsplit(variable_list, "\\s+"))

            # The purpose is not to do model estimation so call it dummy estimation step
            # * * * * Define the title for Mplus model -----
            Mplus.TITLE.obs <- paste0(fragments[[which.fragment]]$name, " dummy estimation on input data;")
            Mplus.TITLE.sim <- paste0(fragments[[which.fragment]]$name, " dummy estimation on simulated data;")

            # The model syntax is the same as MONTE CARLO simulation (except that there is no syntax for categorical variables)
            Mplus.MODEL.syntax <- paste0(
              temp,
              "\n",
              "\n",
              current.fragment.fixed.syntax,
              "\n"
            )

            # Fitting on input data
            # * * * Define the Mplus object ------------
            fragment.est.obs.mplus.object <- mplusObject(
              TITLE = Mplus.TITLE.obs,
              VARIABLE = fragments[[which.fragment]]$variable.syntax,
              MODEL = Mplus.MODEL.syntax,
              OUTPUT = "RESIDUAL;",
              rdata = input.data
            )

            # * * * * Fit the model in Mplus ------------
            fitted.model.mplus <- mplusModeler(
              object = fragment.est.obs.mplus.object,
              hashfilename = FALSE, # do not add a hash to the data file name
              modelout = paste0("fragment", which.fragment, ".iter.", which.iter, ".dummy.est.obs.inp"),
              run = 1L
            )

            # Fitting on simulated data
            # * * * Define the Mplus object ------------
            fragment.est.sim.mplus.object <- mplusObject(
              TITLE = Mplus.TITLE.sim,
              VARIABLE = fragments[[which.fragment]]$variable.syntax,
              MODEL = Mplus.MODEL.syntax,
              OUTPUT = "RESIDUAL;",
              rdata = simulate.data
            )

            # * * * * Fit the model in Mplus ------------
            fitted.model.mplus <- mplusModeler(
              object = fragment.est.sim.mplus.object,
              hashfilename = FALSE, # do not add a hash to the data file name
              modelout = paste0("fragment", which.fragment, ".iter.", which.iter, ".dummy.est.sim.inp"),
              run = 1L
            )

          } # For continuous variables only
        } # Closes loop over each iteration in a batch
    } # Closes loop over estimation batches
  } # Ends the process of simulating data


  # Set working directory to the model checking folder
  # This is added so that if the simulation are done in a previous run, the following process can continue on reading in the output files
  setwd(model.checking.folder)

  #################################
  # Calculate summary statistics on simulated data

  # Set up number of iterations (number of posterior predicted datasets) used in model checking
  n.total.iterations <- n.estimation.batches * iters.per.batch

  if (is.null(n.model.check.iterations)){
    n.model.check.iterations <- n.total.iterations
  }

  if (n.model.check.iterations > n.total.iterations){
    warning("Requested number of model checking iterations exceeds the total number of iterations. Using all available iterations.")
    n.model.check.iterations <- n.total.iterations
  }

  # If not all iterations are used, space out the iterations
  # This does not require n.total.iterations to be divisible by n.model.check.iterations
  model.check.iteration.indices <- floor(
    seq(
      0,
      n.total.iterations,
      length.out = n.model.check.iterations + 1
    )
  )[-1]


  # Set up a list to store the summary data
  summary.simulate.data.list <- vector(
    "list",
    length = n.model.check.iterations
  )

  if (is.null(categorical_variable_syntax)){
    # Set up a list to store realized values of fit statistics
    real.fit.list <- vector(
      "list",
      length = n.model.check.iterations
    )
    # Set up a list to store posterior predictive values of fit statistics
    pred.fit.list <- vector(
      "list",
      length = n.model.check.iterations
    )
  } # For continuous variables only



  for(which.iter in model.check.iteration.indices){

    # * * * * Read in the simulated data and calculate statistics ------------
    simulate.data <- read.table(
      file.path(simulation.files.folder,
      paste0(fragments[[which.fragment]]$name, ".iter.", which.iter, ".simulation.dat")), header = FALSE)

    colnames(simulate.data) <- unlist(strsplit(variable_list, "\\s+"))

    summary.simulate.data.univariate <- simulate.data %>%
      summarise(
      across(everything(), stats),
      iteration = which.iter
      )

    # Pairwise correlations
    cor_result <- pairwise.correlation(simulate.data)
    cor_df <- as.data.frame(as.list(cor_result))

    # Univariate statistics and pairwise correlations
    summary.simulate.data.list[[which.iter]] <-
      bind_cols(
        summary.simulate.data.univariate,
        cor_df
      )

    ##
    # Process fit statistics
    if (is.null(categorical_variable_syntax)){

      # Realized values
      mplus_real_out <- readModels(
        paste0("fragment", which.fragment, ".iter.", which.iter, ".dummy.est.obs.out"),
        what = c("summaries", "residuals"))
      # Posterior predictive values
      mplus_pred_out <- readModels(
        paste0("fragment", which.fragment, ".iter.", which.iter, ".dummy.est.sim.out"),
        what = c("summaries", "residuals"))

      # Extract residual correlations
      real_res_cor_out <- mplus_real_out$residuals$correlationResid
      pred_res_cor_out <- mplus_pred_out$residuals$correlationResid

      # To make sure the variable order is consistent with the input data
      real_res_cor_df <- residual.cor.reorder (
        real_res_cor_out,
        variable_list,
        prefix = "realized_res_cor_"
      )
      pred_res_cor_df <- residual.cor.reorder (
        pred_res_cor_out,
        variable_list,
        prefix = "post_pred_res_cor_"
      )

      # Only keep the pairs of correlations that are observed in input data
      # For realized residual correlations
      cor_names <- sub(
        "^cor_",
        "realized_res_cor_",
        names(input_cor_df)
      )
      real_res_cor_df <- real_res_cor_df %>%
        select(
          all_of(
            names(real_res_cor_df)[
              match(tolower(cor_names), tolower(names(real_res_cor_df)))
            ]
          )
        ) %>%
        mutate(iteration = which.iter, .before = 1)


      # for posterior predictive residual correlations
      cor_names <- sub(
        "^cor_",
        "post_pred_res_cor_",
        names(input_cor_df)
      )
      pred_res_cor_df <- pred_res_cor_df %>%
        select(
          all_of(
            names(pred_res_cor_df)[
              match(tolower(cor_names), tolower(names(pred_res_cor_df)))
            ]
          )
        ) %>%
        mutate(iteration = which.iter, .before = 1)


      # Extract fit indices
      fit_indices <- c("ChiSqM_Value", "CFI", "TLI", "RMSEA_Estimate", "SRMR")

      # realized values
      real_fit_df <- as.data.frame(
        mplus_real_out$summaries[fit_indices]
      ) %>%
        rename_with(
          ~ paste0(
          "realized_",  # add prefix
          sub(
            "^ChiSqM_Value$", # clean ChiSq name
            "Chisq",
            sub("_Estimate$", "", .x) # clean RMSEA name
        )))

      # posterior predictive values
      pred_fit_df <- as.data.frame(
        mplus_pred_out$summaries[fit_indices]
      ) %>%
        rename_with(
          ~ paste0(
            "post_pred_",  # add prefix
            sub(
              "^ChiSqM_Value$", # clean ChiSq name
              "Chisq",
              sub("_Estimate$", "", .x) # clean RMSEA name
        )))

      # Realized values of fit statistics
      real.fit.list[[which.iter]] <- bind_cols(
        real_res_cor_df,
        real_fit_df
      )

      # Posterior predictive values of fit statistics
      pred.fit.list[[which.iter]] <- bind_cols(
        pred_res_cor_df,
        pred_fit_df
      )

    } # For continuous variables only
  } # Closes loop over each iteration


  # Summary statistics for simulated data, of all iterations
  summary.simulate.data <- bind_rows(
    summary.simulate.data.list
  )

  # Ensure that we only report those statistics applicable to the input data
  summary.simulate.data <-
    summary.simulate.data %>%
    select(
      iteration,
      all_of(names(summary.input.data))
    )

  if (is.null(categorical_variable_syntax)){

    # Realized values of fit statistics, of all iterations
    real.fit <- bind_rows(
      real.fit.list
    )
    # Posterior predictive values of fit statistics, of all iterations
    pred.fit <- bind_rows(
      pred.fit.list
    )

    # Combine the above two in one data frame
    real.pred.fit <- bind_cols(
      real.fit,
      pred.fit %>% select(-iteration)
    )

  } # For continuous variables only

  # Calculate ppp value for univariate statistics and pairwise correlations
  summary.simulate.data.clean <- select(summary.simulate.data, -iteration)

  summary.stats.with.ppp <- data.frame(
    statistic = names(summary.simulate.data.clean),
    realized_value = unlist(summary.input.data[1, ]),
    post_pred_mean = colMeans(summary.simulate.data.clean, na.rm = TRUE),
    post_pred_sd = apply(summary.simulate.data.clean, 2, sd, na.rm = TRUE),
    post_pred_min = apply(summary.simulate.data.clean, 2, min, na.rm = TRUE),
    post_pred_max = apply(summary.simulate.data.clean, 2, max, na.rm = TRUE),
    ppp = sapply(names(summary.simulate.data.clean), function(x) {
      mean(summary.simulate.data.clean[[x]] >= summary.input.data[[x]], na.rm = TRUE)
    }),
    row.names = NULL
  )

  if (is.null(categorical_variable_syntax)){
    # Calculate ppp value for fit statistics
    # Check if the variable pairs match
    stopifnot(
       sub("^realized_", "", names(real.fit)) ==
       sub("^post_pred_", "", names(pred.fit))
    )

    fit.ppp <- colMeans(
      select(pred.fit, -iteration) >
        select(real.fit, -iteration)
    )

    real.fit.summary <- data.frame(
      realized_mean = sapply(real.fit[, names(real.fit) != "iteration"], mean),
      realized_sd   = sapply(real.fit[, names(real.fit) != "iteration"], sd),
      realized_min  = sapply(real.fit[, names(real.fit) != "iteration"], min),
      realized_max  = sapply(real.fit[, names(real.fit) != "iteration"], max)
    )
    pred.fit.summary <- data.frame(
      post_pred_mean = sapply(pred.fit[, names(pred.fit) != "iteration"], mean),
      post_pred_sd   = sapply(pred.fit[, names(pred.fit) != "iteration"], sd),
      post_pred_min  = sapply(pred.fit[, names(pred.fit) != "iteration"], min),
      post_pred_max  = sapply(pred.fit[, names(pred.fit) != "iteration"], max)
    )

    fit.summary <- data.frame(
      statistic = sub(
        "^realized_",
        "",
        names(real.fit)[names(real.fit) != "iteration"]
      ),
      real.fit.summary,
      pred.fit.summary,
      ppp = fit.ppp,
      row.names = NULL
    )

  } # For continuous variables only


  ##############
  # Plotting
  if (save.summary.plots.from.model.check){

    post.pred.plot.folder <- paste0(
      model.checking.folder,
      "/Plots/")
    if(!dir.exists(post.pred.plot.folder)) dir.create(post.pred.plot.folder)

    setwd(post.pred.plot.folder)

    plot_variable_list <- names(summary.input.data)

    # A list for univariate statistics
    other_variable_list <- plot_variable_list[
      !grepl("^cor_", plot_variable_list)
    ]

    # Save plots for univariate statistics individually
    for (var.name in other_variable_list) {

      p <- model.checking.plot(
        var.name,
        summary.simulate.data,
        summary.input.data
      )

      ggsave(
        filename = paste0(var.name, ".png"),
        plot = p,
        width = 7,
        height = 5,
        units = "in",
        dpi = 300
      )
    }

    #############
    # Correlations from simulated data
    # Convert to long format, and match with raw variable names
    # Also make sure in the plot, the variables appear in the order as in the dataset
    var_order <- colnames(input.data)

    cor.simulate.long <- summary.simulate.data %>%
      select(starts_with("cor_")) %>%
      tidyr::pivot_longer(
        cols = everything(),
        names_to = "cor_name",
        values_to = "value"
      ) %>%
      inner_join(
        cor_map,
        by = "cor_name"
      ) %>%
      mutate(
        var1 = factor(var1, levels = var_order[var_order %in% var1]),
        var2 = factor(var2, levels = var_order[var_order %in% var2])
      )

    # Correlations from input data
    # Convert to long format, and match with raw variable names
    # Also make sure in the plot, the variables appear in the order as in the dataset
    cor.input.long <- summary.input.data %>%
      select(starts_with("cor_")) %>%
      tidyr::pivot_longer(
        cols = everything(),
        names_to = "cor_name",
        values_to = "realized_value"
      ) %>%
      inner_join(
        cor_map,
        by = "cor_name"
      ) %>%
      mutate(
        var1 = factor(var1, levels = var_order[var_order %in% var1]),
        var2 = factor(var2, levels = var_order[var_order %in% var2])
      )

    p <- model.checking.cor.plot(
      cor.simulate.long,
      cor.input.long
    )

    # Set plot size according to number of variables
    plot.size <- max(8, length(var_order))

    ggsave(
      filename = paste0("Correlations.png"),
      plot = p,
      width = plot.size,
      height = plot.size,
      units = "in",
      dpi = 300
    )

    ##################
    if (is.null(categorical_variable_syntax)){

      # Plot for residual correlations
      res.cor.long <- real.pred.fit %>%
        select(
          iteration,
          matches("^(realized|post_pred)_res_cor_")
        ) %>%
        tidyr::pivot_longer(
          cols = matches(
            "^(realized|post_pred)_res_cor_"
          ),
          names_to = c("type", "correlation"),
          names_pattern = "(realized|post_pred)_res_(.*)", # correlation names in the format of "cor_var1_var2"
          values_to = "value"
        ) %>%
        tidyr::pivot_wider(
          names_from = type, # type: realized or post_pred
          values_from = value
        ) %>%
        mutate(cor_join = tolower(correlation)) %>%
        inner_join(
          cor_map %>%
            mutate(cor_join = tolower(cor_name)),
          by = "cor_join"
        ) %>%
        select(-cor_join) %>%
        mutate(
          var1 = factor(var1, levels = var_order[var_order %in% var1]),
          var2 = factor(var2, levels = var_order[var_order %in% var2])
        )

      p <- model.checking.res.cor.plot(res.cor.long)

      # Set plot size according to number of variables
      plot.size <- max(8, length(var_order))

      ggsave(
        filename = paste0("Residual correlations.png"),
        plot = p,
        width = plot.size,
        height = plot.size,
        units = "in",
        dpi = 300
      )

      ###
      # Plot for fit indices
      fit.long <- real.pred.fit %>%
        select(
          iteration,
          matches("^(realized|post_pred)_(Chisq|CFI|TLI|RMSEA|SRMR)$")
        ) %>%
        tidyr::pivot_longer(
          cols = matches(
            "^(realized|post_pred)_(Chisq|CFI|TLI|RMSEA|SRMR)$"
          ),
          names_to = c("type", "fit"),
          names_pattern = "(realized|post_pred)_(.*)",
          values_to = "value"
        ) %>%
        tidyr::pivot_wider(
          names_from = type,
          values_from = value
        )

      fit_indices_list <- c("Chisq", "CFI", "TLI", "RMSEA", "SRMR")

      # Save plots for fit indices individually
      for (fit.name in fit_indices_list) {

        p <- model.checking.fit.plot(
          fit.long,
          fit.name
        )

        ggsave(
          filename = paste0(fit.name, ".png"),
          plot = p,
          width = 7,
          height = 5,
          units = "in",
          dpi = 300
        )
      }

    } # For continuous variables only


    # Plot heat maps for ppp values
    # For correlations
    cor_ppp_plot <- ppp.plot (
      ppp_data = summary.stats.with.ppp,
      cor_type = "cor",
      cor_map = cor_map,
      var_order = var_order
    )
    ggsave(
      filename = "Correlation ppp values.png",
      plot = cor_ppp_plot,
      width = 7,
      height = 5,
      units = "in",
      dpi = 300
    )

    if (is.null(categorical_variable_syntax)){
      # For residual correlations
      res_cor_ppp_plot <- ppp.plot (
        ppp_data = fit.summary,
        cor_type = "res_cor",
        cor_map = cor_map,
        var_order = var_order
      )
      ggsave(
        filename = "Residual correlation ppp values.png",
        plot = res_cor_ppp_plot,
        width = 7,
        height = 5,
        units = "in",
        dpi = 300
      )
    } # For continuous variables only

  } # Closes plotting

  setwd(model.checking.folder)

  #########################
  # Save to csv

  if (save.summary.stats.from.model.check){

    # Posterior predicted data statistics
    file.name.simulate <- paste0(
      fragments[[which.fragment]]$name,
      " model checking statistics for posterior predicted data",
      ".csv"
    )
    write.csv(
      x = summary.simulate.data,
      file = file.name.simulate,
      row.names = FALSE
    )

    # Input data statistics
    file.name.input <- paste0(
      fragments[[which.fragment]]$name,
      " model checking statistics for input data",
      ".csv"
    )
    write.csv(
      x = summary.input.data,
      file = file.name.input,
      row.names = FALSE
    )


    # Realized values, summary statistics, and ppp values
    file.name.summary <- paste0(
      fragments[[which.fragment]]$name,
      " model checking summary statistics",
      ".csv"
    )
    write.csv(
      x = summary.stats.with.ppp,
      file = file.name.summary,
      row.names = FALSE
    )

    if (is.null(categorical_variable_syntax)){
      # Summary and ppp values for fit statistics
      file.name.fit.summary <- paste0(
        fragments[[which.fragment]]$name,
        " model checking fit summary",
        ".csv"
      )
      write.csv(
        fit.summary,
        file.name.fit.summary,
        row.names = FALSE
      )

      # Fit statistics for each iteration
      file.name.fit <- paste0(
        fragments[[which.fragment]]$name,
        " model checking fit",
        ".csv"
      )
      write.csv(
        real.pred.fit,
        file.name.fit,
        row.names = FALSE
      )

    } # For continuous variables only
  } # Closes save summary stats


  # * * * Copy the simulation datasets into separate folder if desired -----
  if(need.generate.post.pred.data && save.post.pred.data){

    post.pred.data.folder <- paste0(
      model.checking.folder,
      "/Posterior Predicted Data/")
    if(!dir.exists(post.pred.data.folder)) dir.create(post.pred.data.folder)

    files_to_copy <- list.files(pattern = "simulation\\.dat$")

    file.copy(files_to_copy, file.path(post.pred.data.folder, files_to_copy))
  }


  # * * * Delete iteration-specific files if desired -----
  if(need.generate.post.pred.data && !retain.model.check.iteration.files){

    # make a copy of the first iteration, to be saved
    file.name <- list.files(getwd(), pattern = "\\.iter\\.1\\.")

    temp.file.name <- paste0("temp_", seq_along(file.name))

    # copy the file
    file.copy(file.name, temp.file.name)

    # delete files
    to_be_deleted <- list.files(getwd(), pattern = ".iter")
    file.remove(to_be_deleted)

    # rename copy of first iteration
    file.rename(temp.file.name, file.name)

    ##
    if (is.null(categorical_variable_syntax)){
      # * * * Delete dummy estimation files -----
      dummy_files <- list.files(getwd(), pattern = ".dummy")
      file.remove(dummy_files)
    }
  }


  ##
  setwd(home.folder)

  current.time <- Sys.time()
  print(paste0("Model checking for Fragment ", which.fragment, " completed ", round(current.time, units="secs")))

  sink("R Session Info.out", append = TRUE)
  print(paste0("Model checking for Fragment ", which.fragment, " completed ", round(current.time, units="secs")))
  sink()

} # closes MUPPET.model.checking function
