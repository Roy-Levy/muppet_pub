# Function to conduct MUPPET analyses
# Via a fully modular approach


#' @import dplyr
#' @import MCMCvis
#' @import coda
#' @import parallel
#' @import foreach
#' @import doParallel
#' @import stringr
#' @import MplusAutomation
#' @import tinytable
#' @import readr
#' @import tibble
#' @import rlang

.onLoad <- function(libname, pkgname) {
# Ensure rlang is loaded when your package loads, b/c otherwise tinytable won't access it properly to write out tables
library(rlang)
}

#' @import tidyverse
#'
#' @importFrom Hmisc capitalize
#'
#' @title MUPPET Analysis
#'
#' @description  Conducts modular MUPPET modeling.
#'
#' @param n.chains The number of Markov chains to run in each unconditional fragment. Default is 2.

#' @param n.warmup The number of warmup iterations used by the Markov chains. This is applicable for stan and related software. Not needed, and need not be specified, when using Mplus. Default is 0..

#' @param n.burnin The number of iterations to discard as burn-in iterations for the Markov chains. See also the argument n.thin. Default is 500.

#' @param n.iters.per.chain.after.warmup.and.burnin The number of iterations desired, for each chain, after any warmup or burn-in period. See also the argument n.thin. Default is 2500.

#' @param n.thin The thinning parameter, for each Markov chain. Importantly, when using Mplus, this applies to all iterations. For example, suppose n.burnin = 500,  n.iters.per.chain.after.warmup.and.burnin = 2500, and n.thin = 2. This effectively instructs Mplus to run 1000 iterations as burn-in (i.e., 1000 iterations thinned by 2 is 500) and 5000 iterations after burn-in (i.e., 5000 iterations thinned by 2 is 2500). Default is 1.

#' @param n.estimation.batches The number of estimation batches when fitting
#'   fragments that are conditional on antecedent fragments. Default is NULL, in
#'   which case the function will divide the total number post burn-in
#'   iterations by 1000, effectively doing batches of 1000 iterations at a time.
#'   Using more estimation batches leads to fewer iterations per estimation
#'   batch. This can be helpful for reducing computational burden and run time,
#'   and also offers future possibilities for checkpointing.

#' @param convergence.assessment Character, indicating if convergence assessment
#'   is desired for fragments that are not conditional on antecedents. Options
#'   are: "unstandardized" to only plot
#'   unstandardized parameters, "all" to plot standardized and unstandardized
#'   parameters, "none" for turning off convergence assessment. This argument is
#'   irrelevant for fragments that are conditional on antecedent fragments. Default is "unstandardized".

#' @param save.summary.stats.from.MUPPET Logical, indicating if numerical
#'   summaries of posterior distributions be saved. Default is TRUE.

#' @param save.summary.stats.from.MUPPET.in.Word Logical, indicating if slimmed
#'   down numerical summary statistics of the posterior distributions should be
#'   saved as Word documents. Default is FALSE.

#' @param digits.to.round.for.Word Numerical. Digits to round if saving results in Word documents. Only relevant if save.summary.stats.from.MUPPET.in.Word = TRUE. Default = 2.

#' @param save.summary.plots.from.MUPPET Character, indicating if summary plots of the posterior should be saved. Options
#'   are: "unstandardized" to only plot
#'   unstandardized parameters, "all" to plot standardized and unstandardized
#'   parameters, "none" for turning off these plots. Default is "unstandardized".

#' @param save.draws.from.MUPPET Logical, indicating whether the draws from fitting the MUPPET model should be written out as an output file. Default is TRUE.

#' @param model.check Logical, indicating if model-data fit should be conducted. Default is FALSE. This feature was supported in earlier versions for limited classes of models. It is not currently, but is planned for future development.

#' @param save.post.pred.data Logical, indicating if posterior predicted datasets from the model-data fit procedure should be
#'   saved Default is FALSE. This feature was supported in earlier versions
#'   for limited classes of models. It is not currently, but is planned for future
#'   development.

#' @param fragments A list with elements defining the model fragments. Each
#'   element in this list contains a set of specifications for the corresponding
#'   fragment. The 1st element in this list pertains to the 1st fragment, the
#'   2nd element in this list pertains to the 2nd fragment, and so on. The order
#'   of the list also communicates the order each model fragment will be fit. Each element of this list is itself a list with specifications for the fragment. Each set of specifications can include:
#'  \itemize{
  #'  \item{name: Character, giving a name to the fragment. No default. If nothing is supplied, the fragment name will be the fragment number.}
  #'  \item{fragment.folder: Character, giving the name of the folder to be created storing the results for the fragment. If no value is supplied, will default to a name made up of the fragment number and fragment name.}
  #'  \item{model.syntax: Character, containing syntax for the Mplus MODEL statement.}
  #'  \item{variable.syntax: Character, containing syntax for the Mplus VARIABLE statement.}
  #'  \item{priors.syntax: Character, expression containing Mplus syntax for the prior specifications for the model in this fragment, corresponding to the Mplus MODELPRIORS statement. If nothing is specified, defaults to using Mplus's default priors.}
  #'  \item{conditioning: Number, or vector of numbers, indicating which fragments
  #'  to condition on when fitting this fragment. Examples: a value of 0 indicates
  #'  the fragment is unconditional; a value of 1 indicates the fragment is
  #'  conditional on the first fragment in the fragments list; a value of c(1, 2)
  #'  indicates the fragment is conditional on fragment 1 and fragment 2 in the
  #'  fragments list.}
#'  \item{parameters.to.exclude.in.conditioning: A list, of length equal to the number of fragments that are conditioned on (i.e., the length of the collection supplied as the argument for pkg{conditioning}).
#'  Each element of the list is a character, or character vector, naming the kinds of parameters from the associated antecedent fragment that should be ignored in fitting this fragment.
#'  Possible values are:
#'
#'  \itemize{
#'    \item "none"
#'    \item "loadings"
#'    \item "(co)variances of observables"
#'    \item "means and intercepts of observables"
#'    \item "thresholds"
#'    \item "structural coefficients"
#'    \item "(co)variances of latent variables"
#'    \item "means and intercepts of latent variables"
#'  }
#'  A value of "none" indicates that all fitted parameters from the antecedent
#'  fragment should be conditioned on. For example, if pkg{conditioning} = c(1,2),
#'  and pkg{parameters.to.exclude.in.conditioning} = list("none", c("(co)variances of
#'  latent variables", "means and intercepts of latent variables")),then none of
#'  the parameters from fragment 1 will be ignored, while there will be
#'  parameters from fragment 2 that are ignored. These parameters include the
#'  (co)variances of latent variables and the means and intercepts of latent
#'  variables.

#'
#'  If no value for pkg{parameters.to.exclude.in.conditioning}is not supplied, the default is "none" for all antecedent fragments.
#'  }
#'
#'  \item{to.fit: Logical, indicating if this fragment should be fit. If no value is supplied, defaults to TRUE. If value is FALSE, the fragment will not be fit. This can be useful when a fragment has previously been fit, and is now being used as an antecdent to another fragment.}
#'  \item{estimating.lvs: Logical, indicating if this fragment pertains to estimating latent variable values. Default is FALSE.}
#'  \item{data: The dataset to use in fitting this fragment.}


#'}
#'






#' @param retain.iteration.files Logical. For conditional fragments, should the iteration-specific files be saved? Saving all these files may take up a lot of space. If set to FALSE, will still save the files for iteration 1, which may be useful for troubleshooting if needed. Default is FALSE.

#' @param return.R.object Logical, indicating whether to return results as an R object. Default is FALSE. This was supported in earlier versions for a limited class of models. This is not current, and should not be set to TRUE. For future development.


#'
#' @return A statement of how long it took. Future development will include returning objects to R.
#' @export
#'
MUPPET.modular.function <- function(
    # software.environment = "Mplus", # Currently supporting Mplus, others in development

    n.chains = 2,
    n.warmup = 0,
    n.burnin = 500, # For Mplus, thin applies to all iterations. This burn-in is how many of the *thinned* iterations to burn in. So number actually burned in is (n.warmup*n.thin)
    n.iters.per.chain.after.warmup.and.burnin = 2500, # For Mplus, thin applies to all iterations. This is how many of the *thinned* iterations to save as post-burnin.
    n.thin = 1, # For Mplus, thin applies to all iterations. So number actually burned in is (n.burnin*n.thin)

    # number of batches for splitting up estimation of conditional fragments
    n.estimation.batches = NULL,

    # for fragments that are not conditional on antecedents, do we want convergence assessment
    # options are:
    #   "unstandardized" (default if not specified) to only plot unstandardized parameters
    #   "all" to plot standardized and unstandardized parameters
    #   "none" for turning off convergence assessment
    # ignored if the fragment is conditional on antecedents
    convergence.assessment = "unstandardized",


    save.summary.stats.from.MUPPET=TRUE,
    save.summary.stats.from.MUPPET.in.Word=FALSE,
    digits.to.round.for.Word = 2,

    # Summary plots from the posterior
    # options are:
    #   "unstandardized" (default if not specified) to only plot unstandardized parameters
    #   "all" to plot standardized and unstandardized parameters
    #   "none" for turning off these plots
    save.summary.plots.from.MUPPET = "unstandardized",

    save.draws.from.MUPPET=TRUE,
    model.check=FALSE, # Supported in earlier versions, not current, for future development
    save.post.pred.data=FALSE, # Supported in earlier versions, not current, for future development

    fragments = NULL,

    retain.iteration.files = FALSE,     # should we retain all the iteration specific files from Mplus

    return.R.object=FALSE # Supported in earlier versions, not current, for future development
){


  # To delete

  # # Define the operator from rlang that tinytable seems to have trouble importing -----
  # `%||%` <- function(x, y) if (is.null(x)) y else x




  # Sink out file with some session info -----
      sink("R Session Info.out")

  cat(dQuote(paste0("Session Info:"), FALSE), '\n')

  s<-sessionInfo()
  cat(dQuote(paste0(s[["R.version"]][["version.string"]]), FALSE), '\n')
  cat(dQuote(paste0(s[["R.version"]][["platform"]]), FALSE), '\n')
  cat(dQuote(paste0(s[["running"]]), FALSE), '\n')

  # library(parallel)
  cat(dQuote(paste0("Cores: ",parallel::detectCores(),""), FALSE), '\n')
  # cat(dQuote(paste0("Cores: ",strtoi(system("nproc",intern=TRUE))-1,""), FALSE), '\n')

  print(warnings())

  sink()

  # Define the "home" folder -----
  home.folder <- getwd()




  # ***************************** ----
  # IF MODULAR = TRUE ----
  # ***************************** ----

  # if(modular==TRUE){

    # obtain the number of fragments
    n.fragments = length(fragments)


    # For each fragment, define the types of parameters for convergence assessment
    # and
    # If not supplied, define
    #   the fragment to fit as TRUE
    #   the names of the fragments
    #   the names of the fragment specific folders, which are created
    # which.fragment = 1
    for(which.fragment in 1:n.fragments){

      fragments[[which.fragment]]$convergence.assessment = convergence.assessment

      if(is.null(fragments[[which.fragment]]$convergence.assessment)){
        fragments[[which.fragment]]$convergence.assessment = "unstandardized"
      } # closes setting default if convergence assessment argument not supplied



      if(is.null(fragments[[which.fragment]]$to.fit)){
        fragments[[which.fragment]]$to.fit = TRUE
      } #closes if to.fit is not supplied

      if(is.null(fragments[[which.fragment]]$name)){
        fragments[[which.fragment]]$name = paste0("Fragment ", which.fragment)
      } #closes if name is not supplied

      if(is.null(fragments[[which.fragment]]$fragment.folder)){
        fragments[[which.fragment]]$fragment.folder = paste0(home.folder, "/Fragment ", which.fragment, " ", fragments[[which.fragment]]$name, "/")
      } #closes if fragment folder name is not supplied



    } # closes loop over fragments from the supplied list


    # Calculate the number of iterations that will be stored
    # These may need to be split up into batches
    n.iters.total.to.estimate.posterior = n.chains*n.iters.per.chain.after.warmup.and.burnin

    # Set the number of batches for splitting up estimation of conditional fragments
    if(is.null(n.estimation.batches)){
      # n.estimation.batches = 10

      # Based on the number of iterations, decide on the batches
        # if <= 1000, 1 batch
        if(n.iters.total.to.estimate.posterior <= 1000){ n.estimation.batches = 1 }

        # if > 1000, 1 batch per thousand
        if(n.iters.total.to.estimate.posterior > 1000){ n.estimation.batches = n.iters.total.to.estimate.posterior/1000 }

    }

    # Define the iterations per batch
    iters.per.batch = n.iters.total.to.estimate.posterior/n.estimation.batches


    # Declare the software environment -----

    # If using Mplus ------
    # if(software.environment == "Mplus"){


      # ***************************** ----

      # BEGIN FRAGMENTS 1,2,... -----
      # ***************************** ----
      # Loop over fragments ------
      which.fragment = 0
      which.fragment = which.fragment + 1
      for(which.fragment in 1:n.fragments){

        # * Start the timer ----
        start.time <- Sys.time()
        # Save the start time of the first fragment separately
        if(which.fragment == 1){
          start.time.fragment.1 <- start.time
        }



        # * Proceed if this fragment is to be fit ------
        if(fragments[[which.fragment]]$to.fit == TRUE){

          # * * Create the fragment folder ------
          fragment.folder <- fragments[[which.fragment]]$fragment.folder
          if(!dir.exists(fragment.folder)) dir.create(fragment.folder)

          # * * ******************************************************************* ----
          # * * If fragment is NOT conditional on other fragments   ------
          if(sum(fragments[[which.fragment]]$conditioning) == 0){


            # * * * Define the syntax for the ANALYSIS portion of Mplus sytnax -----
            Mplus.ANALYSIS.syntax <-
              paste0("ESTIMATOR = BAYES;",  "\n",
                     paste0("CHAINS = ", n.chains, ";"), "\n",
                     paste0("POINT = MEAN;"), "\n",
                     paste0("FBITERATIONS = ", n.warmup + n.burnin + n.iters.per.chain.after.warmup.and.burnin, ";"),  "\n",
                     paste0("THIN = ", n.thin, ";"),  "\n",
                     paste0("PROCESSORS = ", n.chains*2, ";"),  "\n",
                     paste0("BSEED = ", sample.int(50000,1), ";")

              )


            # * * * Define the title for Mplus model -----
            Mplus.TITLE <- paste0(fragments[[which.fragment]]$name, ";")

            # * * * Define the name of the BPARAMETER file saved from Mplus, which has the draws -----
            BPARAMETER.line <- paste0("BPARAMETER = ", fragments[[which.fragment]]$name, " draws.out;")

            # * * * Define the mplus object for the fragment ------------
            fragment.mplus.object <- mplusObject(
              TITLE = Mplus.TITLE,
              VARIABLE = fragments[[which.fragment]]$variable.syntax,
              ANALYSIS = Mplus.ANALYSIS.syntax,
              MODEL = fragments[[which.fragment]]$model.syntax,
              MODELPRIORS = fragments[[which.fragment]]$priors.syntax,
              SAVEDATA = BPARAMETER.line,


              OUTPUT = "
            STANDARDIZED TECH1 TECH8;
          ",
              PLOT = "
            TYPE = PLOT3;
          ",

              rdata = fragments[[which.fragment]]$data
            )


            # * * * Fit the model in Mplus and read in (some) results ------------
            current.time <- Sys.time()
            print(paste0("Sampling for Fragment ", which.fragment, " initiated ", round(current.time, units="secs")))

            sink("R Session Info.out", append = TRUE)
            print(paste0("Sampling for Fragment ", which.fragment, " initiated ", round(current.time, units="secs")))
            sink()


            fitted.model.mplus <- mplusModeler(
              object = fragment.mplus.object,

              dataout = paste0(fragments[[which.fragment]]$name, " data.dat"),

              writeData = "always", # write out the data if it's not there
              hashfilename = FALSE, # do not add a hash to the data file name

              modelout = paste0(fragments[[which.fragment]]$name, ".inp"),

              run = 1L
            )

            current.time <- Sys.time()
            print(paste0("Sampling for Fragment ", which.fragment, " completed ", round(current.time, units="secs")))

            sink("R Session Info.out", append = TRUE)
            print(paste0("Sampling for Fragment ", which.fragment, " completed ", round(current.time, units="secs")))
            sink()




            # * * * Read in the Mplus output -----
            # Note that this assumes the Mplus .out file has been produced with all lowercase letters
            # Removing here because Linux is case sensitive for this, Windows is not
            # So *don't* lowercase it for Linux
            # Windows should be fine either way

            # Mplus.output <- readModels(
            #   paste0(tolower(fragments.names[which.fragment]),
            #          ".out"))

            Mplus.output <- readModels(
              paste0((fragments[[which.fragment]]$name),
                     ".out")
              )

            # Save the Mplus.output from the fragment as a separate object
            fragment.Mplus.output <- Mplus.output

            # Save the Mplus output object
            file.name <- paste0(fragments[[which.fragment]]$name, ".Mplus output object.rds")
            saveRDS(
              object = fragment.Mplus.output,
              file = file.name
            )


            # Get the names of the file with MCMC draws
            Mplus.MCMC.file.name <- Mplus.output$savedata_info$bayesFile

            # Read in the MCMC draws
            Mplus.bparameters <- read_table(Mplus.MCMC.file.name, col_names = FALSE)
            colnames(Mplus.bparameters) <- Mplus.output$savedata_info$bayesVarNames

            # * * * Construct convergence assessment plots if desired -----
            if(fragments[[which.fragment]]$convergence.assessment != "none"){

            current.time <- Sys.time()
            print(paste0("Convergence assessment for Fragment ", which.fragment, " initiated ", round(current.time, units="secs")))
            sink("R Session Info.out", append = TRUE)
            print(paste0("Convergence assessment for Fragment ", which.fragment, " initiated ", round(current.time, units="secs")))
            sink()

              # Define a folder for convergence assessment results
              convergence.folder <- paste0(fragment.folder, "/Convergence Assessment")
              if(!dir.exists(convergence.folder)) dir.create(convergence.folder)

              # store the draws temporarily
              draws <- Mplus.bparameters

              # Drop the "R-SQUARE" columns
              draws <- draws %>%
                select(
                  !contains("R-SQUARE")
                )

              # Select only unstandardized values if desired
              if(fragments[[which.fragment]]$convergence.assessment == "unstandardized"){
                draws <- draws %>%
                select(
                  !contains("STD")
                )
              } # closes selecting only unstandardized values if desired

              # Need to convert to an MCMC.list
              draws.to.analyze <- df_to_mcmc_list(draws, chain_col = "Chain.number", iter_col = "Iteration.number")

              # Tidy up
              rm(draws)

              # Plot the results
              # trace, densities, and autocorrelations
              # mcmcplots::mcmcplot(
              mcmcplot.mod(
                mcmcout = draws.to.analyze,
                #parms = parameters.for.convergence,
                dir = convergence.folder,
                style="plain",
                filename = "MCMC Plots"
              )

              current.time <- Sys.time()
              print(paste0("Convergence assessment for Fragment ", which.fragment, " completed ", round(current.time, units="secs")))

              sink("R Session Info.out", append = TRUE)
              print(paste0("Convergence assessment for Fragment ", which.fragment, " completed ", round(current.time, units="secs")))
              sink()


              } # closes if doing convergence assessment





            # * * * Select the draws to be those after warmup and burnin -----
            Mplus.bparameters.after.warmup.and.burnin <- dplyr::filter(Mplus.bparameters, Iteration.number > n.warmup + n.burnin)

            # * * * Store and save those draws separately -----
            draws.all.parameters.Mplus <- Mplus.bparameters.after.warmup.and.burnin
            if(save.draws.from.MUPPET) write_csv(
              draws.all.parameters.Mplus,
              paste0(fragments[[which.fragment]]$name, " draws.csv")
            )

            # * * * Define a folder for storing output
            current.time <- Sys.time()
            print(paste0("Posterior summary for Fragment ", which.fragment, " initiated ", round(current.time, units="secs")))

            sink("R Session Info.out", append = TRUE)
            print(paste0("Posterior summary for Fragment ", which.fragment, " initiated ", round(current.time, units="secs")))
            sink()


            posterior.summary.folder <- paste0(fragment.folder, "/Summary of Posterior/")
            if(!dir.exists(posterior.summary.folder)) dir.create(posterior.summary.folder)
            setwd(posterior.summary.folder)


            # Compute the summary statistics for the draws
            draws.to.analyze <- as.mcmc(draws.all.parameters.Mplus)
            summary.statistics.MUPPET <- MCMCsummary(
              draws.to.analyze,
              HPD=TRUE,
              n.eff=FALSE,
              Rhat=FALSE,
              round=8,
              func=median,
              func_name = "median"
            )

            # * * * Write out the summary statistics for the MUPPET model ------
            if(save.summary.stats.from.MUPPET){
              file.name=paste0(
                fragments[[which.fragment]]$name,
                " summary statistics",
                ".csv"
              )
              write.csv(
                x=summary.statistics.MUPPET,
                file=file.name
              )
            } # closes if saving summary statistics from MUPPET model

            # * * * Write out the summary statistics for the MUPPET model in Word ------
            if(save.summary.stats.from.MUPPET.in.Word){

              # Function to remove text up until and including a character
              remove_prefix_function <- function(
              text_vector,
              delimiter # will replace up to (and including) first instance of this
                        ) {
                          sub(paste0("^.*?", delimiter), "", text_vector)
              }

              # Prepare the summary statistics table for publishing ------
              # digits.to.round.for.Word = 2

              # which.fragment = 3
              # fragments[[which.fragment]]$name

              # * Read in the summary statistics table ----
              # file.name <- paste0(fragments[[which.fragment]]$name, " summary statistics.csv")
              # raw.table <-  read_csv(
              #   file = file.name
              # )
              raw.table <- summary.statistics.MUPPET

              # Set the working table
              working.table <- raw.table
              working.table <- rownames_to_column(working.table, var = "parameter")
              # colnames(working.table)[1] = "parameter"

              # * Discarding chain, iterations, and R-squared rows -----
              working.table <- working.table %>%
                filter(
                  !grepl("chain", parameter, ignore.case = TRUE)
                )  %>%
                filter(
                  !grepl("iteration", parameter, ignore.case = TRUE)
                )   %>%
                filter(
                  !grepl("R-SQUARE", parameter, ignore.case = TRUE)
                )



              # * Round to selected digits -----
              working.table <- working.table %>%
                mutate_if(is.numeric, round, digits.to.round.for.Word)

              # * Create interval in single column -----
              working.table$`95% HPDI` <-
                paste0(
                  "(",
                  str_trim(format(working.table$`95%_HPDL`, nmall=digits.to.round.for.Word)), # remove whitespace
                  ", ",
                  str_trim(format(working.table$`95%_HPDU`, nmall=digits.to.round.for.Word)), # remove whitespace
                  ")"
                )


              # * Drop columns not needed -----
              working.table <- select(working.table, -median)
              working.table <- select(working.table, -`95%_HPDL`)
              working.table <- select(working.table, -`95%_HPDU`)

              # * Capitalize columns -----
              colnames(working.table) <- capitalize(colnames(working.table))


              # * Create table ----
              # working.table <- tt(
              #   working.table,
              #   notes = list(
              #     a = list(i = 0, j=ncol(working.table), text = "Highest posterior density interval")
              #   )
              # )


              # * Store the table with all parameters -----
              working.table.all.parameters <- working.table

              # Write out the table
              table.to.write <- working.table.all.parameters

              if(nrow(table.to.write) >0){
                # file.name <- paste0(fragments[[which.fragment]]$name, " summary statistics.docx")
                # tt(table.to.write) |> save_tt(file.name, overwrite=TRUE)

                file.name <- paste0(fragments[[which.fragment]]$name, " summary statistics reduced.csv")
                write.csv(table.to.write, paste0(file.name), row.names = FALSE, quote = TRUE)

              }

              # * * unstandardized solution -----
              if(1==1){

                  working.table <- working.table.all.parameters %>%
                  filter(
                    !grepl("STD", Parameter, ignore.case = TRUE)
                  )


                # * * Remove text before a character from the parameter names -----
                working.table$Parameter <- remove_prefix_function(
                  text_vector = working.table$Parameter,
                  delimiter = "_" # will remove text before (and including) the first instance of this character
                )


                # * * * loadings -----
                loadings.table <- working.table %>%
                  filter(
                    grepl("\\.BY\\.", Parameter, ignore.case = TRUE)
                  )


                # Write out the table
                table.to.write <- loadings.table

                # sort it
                table.to.write <- dplyr::arrange(table.to.write, Parameter)

                sink("Checkpoint.out", append = TRUE)
                print("Checkpoint 7")
                sink()

                # Write out the table
                if(nrow(table.to.write) >0){
                  file.name <- paste0(fragments[[which.fragment]]$name, " unstandardized loadings summary statistics.docx")

                  tt(table.to.write) |> save_tt(paste0(file.name), overwrite=TRUE)

                  # temp.table.tt <- tinytable::tt(table.to.write)

                  # print(paste0("Attempting to save to: ", file.name))

                   # tinytable::save_tt(temp.table.tt, paste0(file.name), overwrite=TRUE)
                  # tinytable::save_tt(temp.table.tt, "test.docx", overwrite=TRUE)

                }

                sink("Checkpoint.out", append = TRUE)
                print("Checkpoint 8")
                sink()

                # * * * intercepts and means -----
                intercepts.and.means.table <- working.table %>%
                  filter(
                    grepl("MEAN\\.", Parameter, ignore.case = TRUE)
                  ) %>%
                  filter(
                    !grepl("\\$", Parameter, ignore.case = TRUE)

                  )

                intercepts.and.means.table$Parameter <- remove_prefix_function(
                  text_vector = intercepts.and.means.table$Parameter,
                  delimiter = "MEAN\\." # will remove text before (and including) the first instance of this character
                )


                # Write out the table
                table.to.write <- intercepts.and.means.table

                # sort it
                # table.to.write <- dplyr::arrange(table.to.write, Parameter)

                if(nrow(table.to.write) >0){
                  file.name <- paste0(fragments[[which.fragment]]$name, " unstandardized intercepts and means summary statistics.docx")
                  tt(table.to.write) |> save_tt(paste0(file.name), overwrite=TRUE)

                }

                # * * * thresholds -----
                thresholds.table <- working.table %>%
                  filter(
                    grepl("MEAN\\.", Parameter, ignore.case = TRUE)
                  ) %>%
                  filter(
                    grepl("\\$", Parameter, ignore.case = TRUE)
                  )

                # Write out the table
                table.to.write <- thresholds.table

                # sort it
                # table.to.write <- dplyr::arrange(table.to.write, Parameter)

                if(nrow(table.to.write) >0){
                  file.name <- paste0(fragments[[which.fragment]]$name, " unstandardized thresholds summary statistics.docx")
                  tt(table.to.write) |> save_tt(paste0(file.name), overwrite=TRUE)

                }





                # * * * covariances -----
                covariances.table <- working.table %>%
                  filter(
                    grepl("\\.WITH\\.", Parameter, ignore.case = TRUE)
                  )

                # Write out the table
                table.to.write <- covariances.table

                # sort it
                table.to.write <- dplyr::arrange(table.to.write, Parameter)

                if(nrow(table.to.write) >0){
                  file.name <- paste0(fragments[[which.fragment]]$name, " unstandardized covariances summary statistics.docx")
                  tt(table.to.write) |> save_tt(paste0(file.name), overwrite=TRUE)

                }




                # * * *  structural coefficents -----
                coefficients.table <- working.table %>%
                  filter(
                    grepl("\\.ON\\.", Parameter, ignore.case = TRUE)
                  )

                # Write out the table
                table.to.write <- coefficients.table

                # sort it
                table.to.write <- dplyr::arrange(table.to.write, Parameter)

                if(nrow(table.to.write) >0){
                  file.name <- paste0(fragments[[which.fragment]]$name, " unstandardized coefficients summary statistics.docx")
                  tt(table.to.write) |> save_tt(paste0(file.name), overwrite=TRUE)

                }



                # * * *  variances and other parameters  -----
                variances.and.others.table <- working.table %>%
                  filter(
                    !grepl("\\.BY\\.", Parameter, ignore.case = TRUE)
                  ) %>%
                  filter(
                    !grepl("MEAN\\.", Parameter, ignore.case = TRUE)
                  ) %>%
                  filter(
                    !grepl("\\.WITH\\.", Parameter, ignore.case = TRUE)
                  ) %>%
                  filter(
                    !grepl(".ON\\.", Parameter, ignore.case = TRUE)
                  )

                # Write out the table
                table.to.write <- variances.and.others.table

                # sort it
                # table.to.write <- dplyr::arrange(table.to.write, Parameter)

                if(nrow(table.to.write) >0){
                  file.name <- paste0(fragments[[which.fragment]]$name, " unstandardized variances and remaining parameters summary statistics.docx")
                  tt(table.to.write) |> save_tt(paste0(file.name), overwrite=TRUE)

                }

              } # closes switch for unstandardized solution



              # * * STDYX solution -----
              if(1==1){
                working.table <- working.table.all.parameters %>%
                  filter(
                    grepl("STDYX", Parameter, ignore.case = TRUE)
                  )

                # * * Remove text before a character from the parameter names -----
                working.table$Parameter <- remove_prefix_function(
                  text_vector = working.table$Parameter,
                  delimiter = "_" # will remove text before (and including) the first instance of this character
                )


                # * * * loadings -----
                loadings.table <- working.table %>%
                  filter(
                    grepl("\\.BY\\.", Parameter, ignore.case = TRUE)
                  )

                # Write out the table
                table.to.write <- loadings.table

                # sort it
                table.to.write <- dplyr::arrange(table.to.write, Parameter)

                # Write out the table
                if(nrow(table.to.write) >0){
                  file.name <- paste0(fragments[[which.fragment]]$name, " standardized loadings summary statistics.docx")
                  tt(table.to.write) |> save_tt(paste0(file.name), overwrite=TRUE)

                }


                # * * * intercepts and means -----
                intercepts.and.means.table <- working.table %>%
                  filter(
                    grepl("MEAN\\.", Parameter, ignore.case = TRUE)
                  ) %>%
                  filter(
                    !grepl("\\$", Parameter, ignore.case = TRUE)

                  )

                intercepts.and.means.table$Parameter <- remove_prefix_function(
                  text_vector = intercepts.and.means.table$Parameter,
                  delimiter = "MEAN\\." # will remove text before (and including) the first instance of this character
                )


                # Write out the table
                table.to.write <- intercepts.and.means.table

                # sort it
                # table.to.write <- dplyr::arrange(table.to.write, Parameter)

                if(nrow(table.to.write) >0){
                  file.name <- paste0(fragments[[which.fragment]]$name, " standardized intercepts and means summary statistics.docx")
                  tt(table.to.write) |> save_tt(paste0(file.name), overwrite=TRUE)

                }

                # * * * thresholds -----
                thresholds.table <- working.table %>%
                  filter(
                    grepl("MEAN\\.", Parameter, ignore.case = TRUE)
                  ) %>%
                  filter(
                    grepl("\\$", Parameter, ignore.case = TRUE)
                  )

                # Write out the table
                table.to.write <- thresholds.table

                # sort it
                # table.to.write <- dplyr::arrange(table.to.write, Parameter)

                if(nrow(table.to.write) >0){
                  file.name <- paste0(fragments[[which.fragment]]$name, " standardized thresholds summary statistics.docx")
                  tt(table.to.write) |> save_tt(paste0(file.name), overwrite=TRUE)

                }





                # * * * covariances -----
                covariances.table <- working.table %>%
                  filter(
                    grepl("\\.WITH\\.", Parameter, ignore.case = TRUE)
                  )

                # Write out the table
                table.to.write <- covariances.table

                # sort it
                table.to.write <- dplyr::arrange(table.to.write, Parameter)

                if(nrow(table.to.write) >0){
                  file.name <- paste0(fragments[[which.fragment]]$name, " standardized covariances summary statistics.docx")
                  tt(table.to.write) |> save_tt(paste0(file.name), overwrite=TRUE)

                }




                # * * *  structural coefficents -----
                coefficients.table <- working.table %>%
                  filter(
                    grepl("\\.ON\\.", Parameter, ignore.case = TRUE)
                  )

                # Write out the table
                table.to.write <- coefficients.table

                # sort it
                table.to.write <- dplyr::arrange(table.to.write, Parameter)

                if(nrow(table.to.write) >0){
                  file.name <- paste0(fragments[[which.fragment]]$name, " standardized coefficients summary statistics.docx")
                  tt(table.to.write) |> save_tt(paste0(file.name), overwrite=TRUE)

                }



                # * * *  variances and other parameters  -----
                variances.and.others.table <- working.table %>%
                  filter(
                    !grepl("\\.BY\\.", Parameter, ignore.case = TRUE)
                  ) %>%
                  filter(
                    !grepl("MEAN\\.", Parameter, ignore.case = TRUE)
                  ) %>%
                  filter(
                    !grepl("\\.WITH\\.", Parameter, ignore.case = TRUE)
                  ) %>%
                  filter(
                    !grepl(".ON\\.", Parameter, ignore.case = TRUE)
                  )

                # Write out the table
                table.to.write <- variances.and.others.table

                # sort it
                # table.to.write <- dplyr::arrange(table.to.write, Parameter)

                if(nrow(table.to.write) >0){
                  file.name <- paste0(fragments[[which.fragment]]$name, " standardized variances and remaining parameters summary statistics.docx")
                  tt(table.to.write) |> save_tt(paste0(file.name), overwrite=TRUE)

                }

              } # closes switch for STDYX solution
            } # closes if saving summary statistics in Word




            # * * * Write out the summary plots for the MUPPET model ------
            if(save.summary.plots.from.MUPPET != "none"){

              # Create a folder for plots
              plots.folder <- paste0(posterior.summary.folder, "/Plots/")
              if(!dir.exists(plots.folder)) dir.create(plots.folder)
              setwd(plots.folder)

              # store the draws temporarily
              draws <- draws.all.parameters.Mplus

              # Select only unstandardized values if desired
              if(save.summary.plots.from.MUPPET == "unstandardized"){
                draws <- draws %>%
                  select(
                    !contains("STD")
                  )
              } # closes selecting only unstandardized values if desired

              # Need to convert to an MCMC.list
              draws.to.analyze <- df_to_mcmc_list(draws, chain_col = "Chain.number", iter_col = "Iteration.number")

              # Tidy up
              rm(draws)

              # Plot the results
              # trace, densities, and autocorrelations
              # mcmcplots::mcmcplot(
              mcmcplot.mod(
                mcmcout = draws.to.analyze,
                #parms = parameters.for.convergence,
                dir = plots.folder,
                style="plain",
                filename = "MCMC Plots"
              )

            } # closes if saving summary plots from MUPPET model

            current.time <- Sys.time()
            print(paste0("Posterior summary for Fragment ", which.fragment, " completed ", round(current.time, units="secs")))

            sink("R Session Info.out", append = TRUE)
            print(paste0("Posterior summary for Fragment ", which.fragment, " completed ", round(current.time, units="secs")))
            sink()


            # Move back to home folder
            setwd(home.folder)

          } # closes if the fragment is NOT conditional on other fragments


          # * * ******************************************************************* ----
          # * * If fragment is conditional on other fragments   ------
          if(sum(fragments[[which.fragment]]$conditioning) != 0){

            current.time <- Sys.time()
            print(paste0("Extracting draws from antecedent fragments for Fragment ", which.fragment, " ", round(current.time, units="secs")))

            sink("R Session Info.out", append = TRUE)
            print(paste0("Extracting draws from antecedent fragments for Fragment ", which.fragment, " ", round(current.time, units="secs")))
            sink()

          # * * * Obtain the draws from the antecedent fragments as a data frame ------
          if(length(fragments[[which.fragment]]$conditioning) > 0){

            which.antecedent.fragment = 0
            which.antecedent.fragment = which.antecedent.fragment + 1
          # for(which.antecedent.fragment in 2:length(fragments[[which.fragment]]$conditioning)){
            for(which.antecedent.fragment in 1:length(fragments[[which.fragment]]$conditioning)){

            # Define the antecedent fragment name and number
            antecedent.fragment.number = fragments[[which.fragment]]$conditioning[which.antecedent.fragment]
            antecedent.fragment.name = fragments[[antecedent.fragment.number]]$name

            # Read in the Mplus output object from that fragment

            # Define the filename if the antecedent fragment is conditional on other fragments
            if(sum(fragments[[antecedent.fragment.number]]$conditioning) == 0){
              file.name <- paste0(
                fragments[[antecedent.fragment.number]]$name, ".Mplus output object.rds"
              )
            } # closes if the antecedent fragment has no antecedents


            if(sum(fragments[[antecedent.fragment.number]]$conditioning) != 0){
              file.name <- paste0(
                fragments[[antecedent.fragment.number]]$name, ".shell Mplus output object.rds"
              )
            } # closes if the antecedent fragment has antecedents


            # file.exists(file.name)
            antecedent.fragment.Mplus.output <- readRDS(file.name)


            # * * * * * Define the parameters to condition on -------
            parameters.to.condition.on <- c(
              "loadings",
              "(co)variances of observables",
              "means and intercepts of observables",
              "thresholds",
              "structural coefficients",
              "(co)variances of latent variables",
              "means and intercepts of latent variables"
            )

            # Now exclude parameters from conditioning if requested
            # Check if there are parameters to exclude
            if(!is.null(fragments[[which.fragment]]$parameters.to.exclude.in.conditioning)){

              # If there are, extract the names of the parameter types to exclude
              temp.names <- unlist(fragments[[which.fragment]]$parameters.to.exclude.in.conditioning[which.antecedent.fragment])

              parameters.to.ignore <- temp.names
              # And now drop them
              parameters.to.condition.on <- parameters.to.condition.on[! parameters.to.condition.on %in% temp.names]

            } # closes if there parameters to condition on


            # Read in existing draws from first antecedent fragment
            temp.data.read.in <- read_csv(
              file = paste0(
                fragments[[antecedent.fragment.number]]$name, " draws.csv")
            )
            # colnames(temp.data.read.in)

            # Drop parameters from the results from the antecedent fragment that are excluded from conditioning in this fragment

            # loadings ('lambda' in Mplus)
            if(sum(parameters.to.condition.on == "loadings") == 0){

              temp.data.read.in <- dplyr::select(
                temp.data.read.in,
                !contains(".BY.")
              )
            } # closes if we are not conditioning on the loadings

            # (co)variances of observables ('theta' in Mplus)
            if(sum(parameters.to.condition.on == "(co)variances of observables") == 0){

              # Extract the parameter specfication from Mplus as a data frame
              parameter.specification <- as.data.frame(
                antecedent.fragment.Mplus.output$tech1$parameterSpecification$theta
              )

              # Proceed if there are any fitted parameters
              if(any(parameter.specification > 0, na.rm = TRUE)){

                # Go through each variable (row in the data frame)
                which.variable = 0
                which.variable = which.variable + 1
                for(which.variable in 1:nrow(parameter.specification)){

                  # Go through each other variable (column in the data frame)
                  which.other.variable = 0
                  which.other.variable = which.other.variable + 1
                  for(which.other.variable in 1:ncol(parameter.specification)){

                    # Proceed if it's a a fitted parameter

                    if(!is.na(parameter.specification[which.variable, which.other.variable]) && parameter.specification[which.variable, which.other.variable] >0){
                      # if(parameter.specification[which.variable, which.other.variable] >0){

                      # If it's a variance
                      if(rownames(parameter.specification)[which.variable] ==
                         colnames(parameter.specification)[which.other.variable]){

                        # Define the parameter name for the Mplus input code
                        antecedent.parameter.name <- paste0(
                          rownames(parameter.specification)[which.variable]
                        )

                          # first define an element of the parameter name in the matrix of antecedent draws
                          antecedent.parameter.name.in.antecedent.draws <- paste0(
                            rownames(parameter.specification)[which.variable]
                          )

                          # Now whittle down the antecendent draws to hopefully leave just one column
                          whittling <- select(
                            temp.data.read.in,
                            contains(antecedent.parameter.name.in.antecedent.draws)
                          )

                          whittling <- select(
                            whittling,
                            !contains(".ON.")
                          )

                          whittling <- select(
                            whittling,
                            !contains(".BY.")
                          )

                          whittling <- select(
                            whittling,
                            !contains("R-SQ")
                          )

                          whittling <- select(
                            whittling,
                            !contains("MEAN.")
                          )

                          whittling <- select(
                            whittling,
                            !contains(".WITH.")
                          )

                          # Now drop the columns that are left after whittling
                          temp.data.read.in <- dplyr::select(
                            temp.data.read.in,
                            !contains(colnames(whittling))
                          )

                      } # closes if it's a variance

                      # If it's a covariance
                      if(rownames(parameter.specification)[which.variable] !=
                         colnames(parameter.specification)[which.other.variable]){


                        # first define an element of the parameter name in the matrix of antecedent draws
                        antecedent.parameter.name.in.antecedent.draws <- paste0(
                          rownames(parameter.specification)[which.variable]
                        )

                        # Now whittle down the antecendent draws to hopefully leave just one column
                        whittling <- select(
                          temp.data.read.in,
                          contains(rownames(parameter.specification)[which.variable])
                        )

                        whittling <- select(
                          whittling,
                          contains(colnames(parameter.specification)[which.other.variable])
                        )

                        whittling <- select(
                          whittling,
                          contains(".WITH.")
                        )

                        # Now drop the columns that are left after whittling
                        temp.data.read.in <- dplyr::select(
                          temp.data.read.in,
                          !contains(colnames(whittling))
                        )

                      } # closes if it's a covariance



                    } # closes if this is a fitted parameter

                  } # closes loop over "other" variable (columns)

                } # closes loop over each variable with fitted parameters (rows)

              } # closes if there are fitted parameters

            } # closes if we are not conditioning on the (co)variances of observables

            # observable means and intercepts ('nu' in Mplus)
            if(sum(parameters.to.condition.on == "means and intercepts of observables") == 0){

              # Extract the parameter specfication from Mplus as a data frame
              parameter.specification <- as.data.frame(
                antecedent.fragment.Mplus.output$tech1$parameterSpecification$nu
              )

              # Proceed if there are any fitted parameters
              if(any(parameter.specification > 0, na.rm = TRUE)){

                # Loop over observables
                which.variable=1
                for(which.variable in 1:ncol(parameter.specification)){

                  name.of.variable <- colnames(parameter.specification)[which.variable]

                  # first define an element of the parameter name in the matrix of antecedent draws
                  antecedent.parameter.name.in.antecedent.draws <- paste0(
                    rownames(parameter.specification)[which.variable]
                  )

                  # Now whittle down the antecendent draws to hopefully leave just one column
                  whittling <- select(
                    temp.data.read.in,
                    contains(name.of.variable)
                  )

                  whittling <- select(
                    whittling,
                    contains("MEAN.")
                  )

                  # But not trying to identify thresholds
                  whittling <- select(
                    whittling,
                    !contains("$")
                  )

                  # Now drop the columns that are left after whittling
                  temp.data.read.in <- dplyr::select(
                    temp.data.read.in,
                    !contains(colnames(whittling))
                  )
                } # closes loop over variables
              } # closes if any parameters are fitted
            } # closes if we are not conditioning on the observable means and intercepts

            # observable thresholds ('tau' in Mplus)
            if(sum(parameters.to.condition.on == "thresholds") == 0){

              # Extract the parameter specfication from Mplus as a data frame
              parameter.specification <- as.data.frame(
                antecedent.fragment.Mplus.output$tech1$parameterSpecification$tau
              )

              # Proceed if there are any fitted parameters
              if(any(parameter.specification > 0, na.rm = TRUE)){

                names.of.observables <- rownames(antecedent.fragment.Mplus.output$tech1$parameterSpecification$lambda)


                # Now whittle down the antecendent draws to hopefully leave just one column
                whittling <- select(
                  temp.data.read.in,
                  contains(names.of.observables)
                )

                whittling <- select(
                  whittling,
                  contains("$")
                )

                  # Now drop the columns that are left after whittling
                  temp.data.read.in <- dplyr::select(
                    temp.data.read.in,
                    !contains(colnames(whittling))
                  )
              } # closes if any parameters are fitted



            } # closes if we are not conditioning on the observable thresholds

            # structural coefficients ('beta' in Mplus)
            if(sum(parameters.to.condition.on == "structural coefficients") == 0){
              temp.data.read.in <- dplyr::select(
                temp.data.read.in,
                !contains(".ON.")
              )
            } # closes if we are not conditioning on the structural coefficients

            # (co)variances of latents ('psi' in Mplus)
            if(sum(parameters.to.condition.on == "(co)variances of latent variables") == 0){

              # Extract the parameter specfication from Mplus as a data frame
              parameter.specification <- as.data.frame(
                antecedent.fragment.Mplus.output$tech1$parameterSpecification$psi
              )

              # Proceed if there are any fitted parameters

                # Go through each variable (row in the data frame)
                which.variable = 0
                which.variable = which.variable + 1
                for(which.variable in 1:nrow(parameter.specification)){

                  # Go through each other variable (column in the data frame)
                  which.other.variable = 0
                  which.other.variable = which.other.variable + 1
                  for(which.other.variable in 1:ncol(parameter.specification)){

                      # If it's a variance
                      if(rownames(parameter.specification)[which.variable] ==
                         colnames(parameter.specification)[which.other.variable]){

                        # Define the parameter name for the Mplus input code
                        # antecedent.parameter.name <- paste0(
                        #   rownames(parameter.specification)[which.variable]
                        # )

                        # first define an element of the parameter name in the matrix of antecedent draws
                        antecedent.parameter.name.in.antecedent.draws <- paste0(
                          rownames(parameter.specification)[which.variable]
                        )

                        # Now whittle down the antecendent draws to hopefully leave just one column
                        whittling <- select(
                          temp.data.read.in,
                          contains(antecedent.parameter.name.in.antecedent.draws)
                        )

                        whittling <- select(
                          whittling,
                          !contains(".ON.")
                        )

                        whittling <- select(
                          whittling,
                          !contains(".BY.")
                        )

                        whittling <- select(
                          whittling,
                          !contains("R-SQ")
                        )

                        whittling <- select(
                          whittling,
                          !contains("MEAN.")
                        )

                        whittling <- select(
                          whittling,
                          !contains(".WITH.")
                        )

                        # Now drop the columns that are left after whittling
                        temp.data.read.in <- dplyr::select(
                          temp.data.read.in,
                          !contains(colnames(whittling))
                        )

                      } # closes if it's a variance

                      # If it's a covariance
                      if(rownames(parameter.specification)[which.variable] !=
                         colnames(parameter.specification)[which.other.variable]){


                        # first define an element of the parameter name in the matrix of antecedent draws
                        antecedent.parameter.name.in.antecedent.draws <- paste0(
                          rownames(parameter.specification)[which.variable]
                        )

                        # Now whittle down the antecendent draws to hopefully leave just one column
                        whittling <- select(
                          temp.data.read.in,
                          contains(rownames(parameter.specification)[which.variable])
                        )

                        whittling <- select(
                          whittling,
                          contains(colnames(parameter.specification)[which.other.variable])
                        )

                        whittling <- select(
                          whittling,
                          contains(".WITH.")
                        )

                        # Now drop the columns that are left after whittling
                        temp.data.read.in <- dplyr::select(
                          temp.data.read.in,
                          !contains(colnames(whittling))
                        )

                      } # closes if it's a covariance

                  } # closes loop over "other" variable (columns)

                } # closes loop over each variable  (rows)

            } # closes if we are not conditioning on the (co)variances of latents

            # latent means and intercepts ('alpha' in Mplus)
            if(sum(parameters.to.condition.on == "means and intercepts of latent variables") == 0){

              # Extract the parameter specfication from Mplus as a data frame
              parameter.specification <- as.data.frame(
                antecedent.fragment.Mplus.output$tech1$parameterSpecification$alpha
              )

              names.of.latents <- colnames(parameter.specification)

              # Now whittle down the antecendent draws to hopefully leave just one column
              whittling <- select(
                temp.data.read.in,
                contains(names.of.latents)
              )

              whittling <- select(
                whittling,
                contains("MEAN.")
              )

              whittling <- select(
                whittling,
                !contains("$")
              )

              # Now drop the columns that are left after whittling
              temp.data.read.in <- dplyr::select(
                temp.data.read.in,
                !contains(colnames(whittling))
              )

            } # closes if we are not conditioning on the latent means and intercepts




            # Store those draws if this is the first antecendent fragment
            if(which.antecedent.fragment == 1){
              draws.from.antecedent.fragments.as.data.frame <- temp.data.read.in
            } # closes if this is the first antecedent fragment


          # colnames(temp.data.read.in)

          # Append the draws from this antecendent fragment
          # draws.from.antecedent.fragments.as.data.frame <- bind_cols(draws.from.antecedent.fragments.as.data.frame, temp.data.read.in)
          # colnames(draws.from.antecedent.fragments.as.data.frame)

          # Add only columns from the current antecedent that are not already present from past antecendents
          if(which.antecedent.fragment > 1){
          hope <- add_column(
              draws.from.antecedent.fragments.as.data.frame,
              !!!temp.data.read.in[setdiff(names(temp.data.read.in), names(draws.from.antecedent.fragments.as.data.frame))]
            )
          draws.from.antecedent.fragments.as.data.frame <- hope
          rm(hope)
          } # closes if this is the not the first antecedent fragment

          # select(temp.data.read.in, contains("F2.ON"))
          # select(draws.from.antecedent.fragments.as.data.frame, contains("F2.ON"))
          # select(draws.from.antecedent.fragments.as.data.frame, contains(".x"))
          # select(draws.from.antecedent.fragments.as.data.frame, contains("F1.BY.CH1IT10"))
          # select(temp.data.read.in, contains("F1.BY.CH1IT10"))
          # select(from.current.to.append, contains("F1.BY.CH1IT10"))

          } # closes loop over antecedent fragments for reading in draws

          } # closes if there are antecedent fragments


          # * * * Only keep unstandardized by discarding -----
          unstandardized.draws.from.antecedent.fragments.as.data.frame <- draws.from.antecedent.fragments.as.data.frame %>%
            select(!contains("STD"))

          # colnames(unstandardized.draws.from.antecedent.fragments.as.data.frame)
          # select(unstandardized.draws.from.antecedent.fragments.as.data.frame, contains("F2.ON"))
          # select(unstandardized.draws.from.antecedent.fragments.as.data.frame, contains("F1.BY.CH1IT10"))



          # * * * If estimating LVs, define the data to be used for estimating LVs -----
          if(fragments[[which.fragment]]$estimating.lvs == TRUE){
            # Extract the data, keeping on the distinct rows
            # These are unique response vectors

            # Keep only distinct rows
            distinct.data <- distinct(fragments[[which.fragment]]$data)
          } # closes if estimating LVs in this fragment




        # * * * ******************************************************************* ----
        current.time <- Sys.time()
        print(paste0("Sampling for Fragment ", which.fragment, " initiated ", round(current.time, units="secs")))

        sink("R Session Info.out", append = TRUE)
        print(paste0("Sampling for Fragment ", which.fragment, " initiated ", round(current.time, units="secs")))
        sink()


        # * * * Clean up existing parallel environment ------
        unregister_dopar <- function() {
          env <- foreach:::.foreachGlobals
          rm(list=ls(name=env), pos=env)
        }

        # unregister_dopar()

        # * * * Setup parallel environment ------

        # * * * Detect and set up the number of cores
        numCores <- parallel::detectCores()-1

        # strtoi works on linux
        #numCores <- strtoi(system("nproc",intern=TRUE))-1

        #future::plan(multicore, workers=numCores)
        registerDoParallel(numCores)


        # * * Set the number of cores to use
        #registerDoParallel(makeCluster(numCores, outfile=""))  # use multicore, set to the number of our cores
        #plan(multisession, workers = numCores)
        #registerDoFuture()


        # * * * Initiate object for storing draws for this fragment ------
        draws.from.MUPPET.model.this.fragment <- NULL

        # * * * Loop over estimation batches ------

        which.estimation.batch = 0
        which.estimation.batch = which.estimation.batch +1


        for(which.estimation.batch in 1:n.estimation.batches){
        # for(which.estimation.batch in 1:2){

          current.time <- Sys.time()
          print(paste0("Starting Fragment ", which.fragment, " batch ", which.estimation.batch, " of ", n.estimation.batches, " ", round(current.time, units="secs")))

          sink("R Session Info.out", append = TRUE)
          print(paste0("Starting Fragment ", which.fragment, " batch ", which.estimation.batch, " of ", n.estimation.batches, " ", round(current.time, units="secs")))
          sink()

          # Define the start and end iteration for each batch
          start.iter = (which.estimation.batch-1)*iters.per.batch + 1
          end.iter = (which.estimation.batch-1)*iters.per.batch + iters.per.batch
          # start.iter = 11
          # end.iter = start.iter + 9



        # * * * Run analyses in parallel ------
        # which.iter=0
        # which.iter=which.iter+1
        # draws.from.MUPPET.model.this.fragment <- foreach(which.iter=1:13,
        #draws.from.MUPPET.model <- foreach(which.iter=1:nrow(draws.from.antecedent.fragments.as.data.frame),
        # temp <- foreach(which.iter=1:3,



        # commenting this out to now do batches
        # draws.from.MUPPET.model.this.fragment <- foreach(which.iter=1:nrow(draws.from.antecedent.fragments.as.data.frame),

        draws.from.MUPPET.model.this.batch <- foreach(which.iter=start.iter:end.iter,
          #results.from.parallel <- foreach(which.iter=1:1000,
         .packages = c("dplyr", "MplusAutomation", "tidyverse"),
         #.combine = 'list',
         .combine = rbind
         #.options.future = list("dplyr", "MplusAutomation", "tidyverse", scheduling = Inf)
         #) %dofuture% {
        ) %dopar% {

          which.iter = as.integer(which.iter)
          #print(paste("iter", which.iter, "--start"))


          #.combine=rbind) %do% {


           # * * *  * ******************************************************************* ----
           # * * * * Define the Mplus MODEL portion of the input file -----

           # * * * * * Initiate an empty object. Will append text ------
           temp <- NULL

           # * * * * * Initiate an empty object to store list of parameter names from antecedents
           # antecedent.parameter.names <- c("")
           antecedent.parameter.names <- NULL


           # * * * * * Initiate an object capturing if there is a covariance from an antecedent fragment treated as fixed in Mplus
           # (by use of a 'with' statement)
           # In Mplus, if there's a fixed covariance s.t. the covariance matrix cannot be partitioned into uncorrelated blocks,
           # need to use a different sampler
           fixed.covariance.from.antecedent <- FALSE


           # Loop over antecedent fragments
           which.antecedent.fragment = 0
           which.antecedent.fragment = which.antecedent.fragment + 1
           for(which.antecedent.fragment in 1:length(fragments[[which.fragment]]$conditioning)){


             # Define the antecedent fragment number and name
             antecedent.fragment.number = fragments[[which.fragment]]$conditioning[which.antecedent.fragment]
             antecedent.fragment.name = fragments[[antecedent.fragment.number]]$name

             # Read in the Mplus output object from that fragment

             # Define the filename if the antecedent fragment is not conditional on other fragments
             # if(sum(fragments.conditioning[[antecedent.fragment.number]]) == 0){
             #   file.name <- paste0(
             #     fragments.names[antecedent.fragment.number], " Mplus output object.rds"
             #   )
             # } # closes if the antecedent fragment has no antecedents

             # Define the filename if the antecedent fragment is conditional on other fragments
             if(sum(fragments[[antecedent.fragment.number]]$conditioning) == 0){
               file.name <- paste0(
                 fragments[[antecedent.fragment.number]]$name, ".Mplus output object.rds"
               )
             } # closes if the antecedent fragment has no antecedents


             # if(sum(fragments.conditioning[[antecedent.fragment.number]]) != 0){
             #   file.name <- paste0(
             #     fragments.names[antecedent.fragment.number], ".shell Mplus output object.rds"
             #   )
             # } # closes if the antecedent fragment has no antecedents


             if(sum(fragments[[antecedent.fragment.number]]$conditioning) != 0){
               file.name <- paste0(
                 fragments[[antecedent.fragment.number]]$name, ".shell Mplus output object.rds"
               )
             } # closes if the antecedent fragment has no antecedents


             # file.exists(file.name)
             antecedent.fragment.Mplus.output <- readRDS(file.name)
              # antecedent.fragment.Mplus.output$tech1$parameterSpecification$lambda

             # # Read in existing draws from that fragment
             # draws.from.antecedent.fragment.as.data.frame <- read_csv(
             #   file = paste0(fragments.names[antecedent.fragment.number], " unstandardized draws.csv")
             # )


             # Define the line to append to syntax
             this.line <- paste0(
               "! Parameters from fragment ",
               antecedent.fragment.number
             )

             # Append the line
             temp <- paste0(
               temp, "\n",
               "\n",
               this.line
             )


             # * * * * * Define the parameters to condition on -------
             parameters.to.condition.on <- c(
               "loadings",
               "(co)variances of observables",
               "means and intercepts of observables",
               "thresholds",
               "structural coefficients",
               "(co)variances of latent variables",
               "means and intercepts of latent variables"
             )

             # Now exclude parameters from conditioning if requested
             # Check if there are parameters to exclude
             if(!is.null(fragments[[which.fragment]]$parameters.to.exclude.in.conditioning)){

               # If there are, extract the names of the parameter types to exclude
               temp.names <- unlist(fragments[[which.fragment]]$parameters.to.exclude.in.conditioning[which.antecedent.fragment])

               # And now drop them
               parameters.to.condition.on <- parameters.to.condition.on[! parameters.to.condition.on %in% temp.names]

             } # closes if there parameters to condition on

             # * * * * * Add code for loadings ('lambda' in Mplus)-------
             if(sum(parameters.to.condition.on == "loadings") == 1){

               # Extract the parameter specfication from Mplus as a data frame
               parameter.specification <- as.data.frame(
                 antecedent.fragment.Mplus.output$tech1$parameterSpecification$lambda
               )

               # Proceed if there are any fitted parameters
               if(any(parameter.specification > 0, na.rm = TRUE)){

                 # Go through each variable (row in the data frame)
                 which.variable = 1
                 for(which.variable in 1:nrow(parameter.specification)){

                   # Go through each other variable (column in the data frame)
                   which.other.variable = 1
                   for(which.other.variable in 1:ncol(parameter.specification)){

                     # Proceed if it's a a fitted parameter
                     if(!is.na(parameter.specification[which.variable, which.other.variable]) && parameter.specification[which.variable, which.other.variable] >0){

                       # if(parameter.specification[which.variable, which.other.variable] >0){

                       # Define the parameter name for the Mplus input code
                       antecedent.parameter.name <- paste0(
                         colnames(parameter.specification)[which.other.variable],
                         " by ",
                         rownames(parameter.specification)[which.variable]

                       )

                       # Proceed if the parameter name is new to the list of antecedent parameters, don't if it already exists among list of antecedent parameters
                       if(sum(antecedent.parameter.name==antecedent.parameter.names)==0){

                         # Define the value of the entity

                         # first define an element of the parameter name in the matrix of antecedent draws
                         antecedent.parameter.name.in.antecedent.draws <- paste0(
                           colnames(parameter.specification)[which.other.variable],
                           ".by.",
                           rownames(parameter.specification)[which.variable]
                         )

                         # Pull out the set of draws for that parameter
                         #select(draws.from.antecedent.fragment.as.data.frame, contains(antecedent.parameter.name.in.antecedent.draws))
                         # pull(draws.from.antecedent.fragment.as.data.frame, contains(antecedent.parameter.name.in.antecedent.draws))
                         # colnames(draws.from.antecedent.fragments.as.data.frame)

                         # colnames(unstandardized.draws.from.antecedent.fragments.as.data.frame)
                         # value.to.paste <- pull(
                         #   unstandardized.draws.from.antecedent.fragments.as.data.frame,
                         #   contains(antecedent.parameter.name.in.antecedent.draws)
                         # )[which.iter]


                         # The 'pull' function above sometimes returns an error
                         # says it's not a unique column
                         # I think it's because the name of the variable is a subset of others:
                         # E.g., 'Item1' is a subset of 'Item10' so it can't pick a unique column for the first
                         # Using the 'select' function below, which returns a table
                         # This seems to just return a table with 1 column, so it works
                         # Not sure why 'pull' indicates there are multiple columns, but 'select' ony yields one column
                         value.to.paste <- select(
                           unstandardized.draws.from.antecedent.fragments.as.data.frame,
                           contains(antecedent.parameter.name.in.antecedent.draws)
                         )[which.iter,1]

                         # If the value is equal to 0, perturb it a little
                         # Need to do this because if it's exactly equal to 0, Mplus will fix it properly
                         # but won't include it as one of the parameters saved in the Bayesian output
                         if(value.to.paste==0){

                           # Define a U(.0001, .0010) random variable to increment or decrement
                           epsilon <- runif(1, .0001, .0010)

                           # Add or subtract the epsilon based on a coin flip
                           should.add <- rbinom(1, 1, .5)==1
                           if(should.add) {value.to.paste = value.to.paste + epsilon }
                           if(!should.add) {value.to.paste = value.to.paste - epsilon }
                         } # closes if the value from the measurement model was exactly 0



                         # Define the line to append to syntax
                         this.line <- paste0(
                           # name.of.latent, " by ",
                           # rownames(measurement.model.Mplus.output$tech1$parameterSpecification$lambda)[which.observable],
                           #rownames(antecedent.fragment.Mplus.output$tech1$parameterSpecification$lambda)[which.observable],
                           antecedent.parameter.name,
                           " @",

                           # measurement.model.draws.as.data.frame[which.iter, measurement.model.Mplus.output$tech1$parameterSpecification$lambda[which.observable, which.latent]],
                           value.to.paste,
                           ";"

                         )

                         # Append the line
                         temp <- paste0(
                           temp, "\n",
                           this.line
                         )

                         # Add this parameter name to the collection of existing parameter names
                         antecedent.parameter.names <- c(antecedent.parameter.names, antecedent.parameter.name)


                       } # closes if the parameter is new to the list of antecedent parameters

                     } # closes if this is a fitted parameter

                   } # closes loop over "other" variable (columns)

                 } # closes loop over each variable with fitted parameters (rows)

               } # closes if there are fitted loadings

             } # closes switch for adding code for loadings

             # * * * * * Add code for (co)variances of observables ('theta' in Mplus)  -------
             if(sum(parameters.to.condition.on == "(co)variances of observables") == 1){

               # Extract the parameter specfication from Mplus as a data frame
               parameter.specification <- as.data.frame(
                 antecedent.fragment.Mplus.output$tech1$parameterSpecification$theta
               )

               # Proceed if there are any fitted parameters
               if(any(parameter.specification > 0, na.rm = TRUE)){

                 # Go through each variable (row in the data frame)
                 which.variable = 0
                 which.variable = which.variable + 1
                 for(which.variable in 1:nrow(parameter.specification)){

                   # Go through each other variable (column in the data frame)
                   which.other.variable = 0
                   which.other.variable = which.other.variable + 1
                   for(which.other.variable in 1:ncol(parameter.specification)){

                     # Proceed if it's a a fitted parameter

                     if(!is.na(parameter.specification[which.variable, which.other.variable]) && parameter.specification[which.variable, which.other.variable] >0){
                       # if(parameter.specification[which.variable, which.other.variable] >0){

                       # If it's a variance
                       if(rownames(parameter.specification)[which.variable] ==
                          colnames(parameter.specification)[which.other.variable]){

                         # Define the parameter name for the Mplus input code
                         antecedent.parameter.name <- paste0(
                           rownames(parameter.specification)[which.variable]
                         )


                         # Proceed if the parameter name is new to the list of antecedent parameters, don't if it already exists among list of antecedent parameters
                         if(sum(antecedent.parameter.name==antecedent.parameter.names)==0){

                           # Define the value of the entity

                           # first define an element of the parameter name in the matrix of antecedent draws
                           antecedent.parameter.name.in.antecedent.draws <- paste0(
                             rownames(parameter.specification)[which.variable]
                           )

                           # Now whittle down the antecendent draws to hopefully leave just one column
                           whittling <- select(
                             unstandardized.draws.from.antecedent.fragments.as.data.frame,
                             contains(antecedent.parameter.name.in.antecedent.draws)
                           )

                           whittling <- select(
                             whittling,
                             !contains("ON")
                           )

                           whittling <- select(
                             whittling,
                             !contains("BY")
                           )

                           whittling <- select(
                             whittling,
                             !contains("R-SQ")
                           )

                           whittling <- select(
                             whittling,
                             !contains("MEAN")
                           )

                           # paste the value if we've identified the isolated column
                           if(ncol(whittling) == 1) {
                             value.to.paste <- whittling[which.iter, 1]
                           } # closes if just one column in whittling


                           # If the value is equal to 0, perturb it a little
                           # Need to do this because if it's exactly equal to 0, Mplus will fix it properly
                           # but won't include it as one of the parameters saved in the Bayesian output
                           if(value.to.paste==0){

                             # Define a U(.0001, .0010) random variable to increment or decrement
                             epsilon <- runif(1, .0001, .0010)

                             # Add or subtract the epsilon based on a coin flip
                             should.add <- rbinom(1, 1, .5)==1
                             if(should.add) {value.to.paste = value.to.paste + epsilon }
                             if(!should.add) {value.to.paste = value.to.paste - epsilon }
                           } # closes if the value from the measurement model was exactly 0



                           # Define the line to append to syntax
                           this.line <- paste0(
                             antecedent.parameter.name,
                             " @",
                             value.to.paste,
                             ";"

                           )

                           # Append the line
                           temp <- paste0(
                             temp, "\n",
                             this.line
                           )

                           # Add this parameter name to the collection of existing parameter names
                           antecedent.parameter.names <- c(antecedent.parameter.names, antecedent.parameter.name)


                         } # closes if the parameter is new to the list of antecedent parameters


                       } # closes if it's a variance

                       # If it's a covariance
                       if(rownames(parameter.specification)[which.variable] !=
                          colnames(parameter.specification)[which.other.variable]){


                         # Define the parameter name for the Mplus input code
                         antecedent.parameter.name <- paste0(
                           colnames(parameter.specification)[which.variable],
                           " with ",
                           rownames(parameter.specification)[which.other.variable]

                         )

                         # Proceed if the parameter name is new to the list of antecedent parameters, don't if it already exists among list of antecedent parameters
                         if(sum(antecedent.parameter.name==antecedent.parameter.names)==0){

                           # Define the value of the entity

                           # first define an element of the parameter name in the matrix of antecedent draws
                           antecedent.parameter.name.in.antecedent.draws <- paste0(
                             rownames(parameter.specification)[which.variable],
                             ".with.",
                             colnames(parameter.specification)[which.other.variable]
                           )

                           # check if that parameter name is in the matrix of antecedent draws
                           # If it ISN'T, reverse the ordering of the variables in the name
                           # If it isn't see below
                           if(ncol(
                             select(
                               unstandardized.draws.from.antecedent.fragments.as.data.frame,
                               contains(antecedent.parameter.name.in.antecedent.draws)
                             )
                           ) == 0){
                             antecedent.parameter.name.in.antecedent.draws <- paste0(
                               colnames(parameter.specification)[which.other.variable],
                               ".with.",
                               rownames(parameter.specification)[which.variable]
                             )
                           } # closes if the parameter name is in the matrix of antecedent draws


                           # The 'pull' function above sometimes returns an error
                           # says it's not a unique column
                           # I think it's because the name of the variable is a subset of others:
                           # E.g., 'Item1' is a subset of 'Item10' so it can't pick a unique column for the first
                           # Using the 'select' function below, which returns a table
                           # This seems to just return a table with 1 column, so it works
                           # Not sure why 'pull' indicates there are multiple columns, but 'select' ony yields one column
                           value.to.paste <- select(
                             unstandardized.draws.from.antecedent.fragments.as.data.frame,
                             contains(antecedent.parameter.name.in.antecedent.draws)
                           )[which.iter,1]


                           # If the value is equal to 0, perturb it a little
                           # Need to do this because if it's exactly equal to 0, Mplus will fix it properly
                           # but won't include it as one of the parameters saved in the Bayesian output
                           if(value.to.paste==0){

                             # Define a U(.0001, .0010) random variable to increment or decrement
                             epsilon <- runif(1, .0001, .0010)

                             # Add or subtract the epsilon based on a coin flip
                             should.add <- rbinom(1, 1, .5)==1
                             if(should.add) {value.to.paste = value.to.paste + epsilon }
                             if(!should.add) {value.to.paste = value.to.paste - epsilon }
                           } # closes if the value from the measurement model was exactly 0

                           # Define the line to append to syntax
                           this.line <- paste0(
                             antecedent.parameter.name,
                             " @",
                             value.to.paste,
                             ";"

                           )



                           # Append the line
                           temp <- paste0(
                             temp, "\n",
                             this.line
                           )

                           # Add this parameter name to the collection of existing parameter names
                           antecedent.parameter.names <- c(antecedent.parameter.names, antecedent.parameter.name)





                         } # closes if the parameter is new to the list of antecedent parameters

                         # Since specifying a covariance as fixed, declare that there is a fixed covariance
                         # In Mplus, if there's a fixed covariance s.t. the covariance matrix cannot be partitioned into uncorrelated blocks,
                         # need to use a different sampler
                         fixed.covariance.from.antecedent <- TRUE

                       } # closes if it's a covariance



                     } # closes if this is a fitted parameter

                   } # closes loop over "other" variable (columns)

                 } # closes loop over each variable with fitted parameters (rows)

               } # closes if there are fitted parameters
             } # closes switch for adding code for (co)variances of observables

             # * * * * * Add code for observable means and intercepts ('nu' in Mplus)-------
             if(sum(parameters.to.condition.on == "means and intercepts of observables") == 1){

               # Extract the parameter specfication from Mplus as a data frame
               parameter.specification <- as.data.frame(
                 antecedent.fragment.Mplus.output$tech1$parameterSpecification$nu
               )

               # Proceed if there are any fitted parameters
               if(any(parameter.specification > 0, na.rm = TRUE)){

                 # Loop over observables
                 which.variable=1
                 for(which.variable in 1:ncol(parameter.specification)){

                   # proceed if this parameter is fitted (i.e., the value of the parameter specification table is > 0)
                   if(!is.na(parameter.specification[which.variable]) && parameter.specification[which.variable] >0){
                     # if(parameter.specification[1, which.variable] != 0){



                     name.of.variable <- colnames(parameter.specification)[which.variable]

                     # first define an element of the parameter name in the matrix of antecedent draws
                     antecedent.parameter.name.in.antecedent.draws <- paste0(
                       "MEAN.", name.of.variable
                     )

                     # select the draws for the intercept for this observable
                     draws.for.intercept <- select(unstandardized.draws.from.antecedent.fragments.as.data.frame, contains(antecedent.parameter.name.in.antecedent.draws))


                     # Define the parameter name
                     antecedent.parameter.name <- paste0(
                       "[", name.of.variable,"]"
                     )

                     # Proceed if the parameter name is new, don't if it already exists among list of antecedent parameters
                     if(sum(antecedent.parameter.name==antecedent.parameter.names)==0){

                       # Define the value of the entity
                       # Previous code that referred to the matrix via a column number
                       # Should work for taking values from the draws for the first antecedent
                       # But may not work in general
                       # value.to.paste <- draws.from.antecedent.fragments.as.data.frame[which.iter, param.number.thresholds.this.observable[which.threshold]]
                       # value.to.paste <- draws.from.antecedent.fragment.as.data.frame[which.iter, param.number.thresholds.this.observable[which.threshold]]

                       # this should work regardless of which antecedent

                       value.to.paste <- draws.for.intercept[which.iter, 1]


                       # If the value is equal to 0, perturb it a little
                       # Need to do this because if it's exactly equal to 0, Mplus will fix it properly
                       # but won't include it as one of the parameters saved in the Bayesian output
                       if(value.to.paste==0){

                         # Define a U(.0001, .0010) random variable to increment or decrement
                         epsilon <- runif(1, .0001, .0010)

                         # Add or subtract the epsilon based on a coin flip
                         should.add <- rbinom(1, 1, .5)==1
                         if(should.add) {value.to.paste = value.to.paste + epsilon }
                         if(!should.add) {value.to.paste = value.to.paste - epsilon }

                       } # closes if the value from the measurement model was exactly 0


                       # Define the line to append to syntax
                       # Usual way won't work b/c in Mplus the constraint is inside the brackets
                       # this.line <- paste0(
                       #   antecedent.parameter.name,
                       #   "@",
                       #   value.to.paste,
                       #   ";"
                       # )

                       this.line <- paste0(
                         "[", name.of.variable,"@",
                         value.to.paste,
                         "];"
                       )

                       # Append the line
                       temp <- paste0(
                         temp, "\n",
                         this.line
                       )

                       # Add this parameter name to the collection of existing parameter names
                       antecedent.parameter.names <- c(antecedent.parameter.names, antecedent.parameter.name)

                     } # closes if the parameter is new to the list of antecedent parameters
                   } # closes if this parameter is fitted (i.e., the value of the parameter specification table is > 0)


                 } # closes loop over variable
               } # closes if there are fitted observable means and intercepts
             } # closes switch for adding code for observable means and intercepts

             # * * * * * Add code for thresholds ('tau' in Mplus) -------
             if(sum(parameters.to.condition.on == "thresholds") == 1){

               # Extract the parameter specfication from Mplus as a data frame
               parameter.specification <- as.data.frame(
                 antecedent.fragment.Mplus.output$tech1$parameterSpecification$tau
               )

               # Proceed if there are any fitted parameters
               if(any(parameter.specification > 0, na.rm = TRUE)){

                 # Loop over observables
                 # Note the name of the observables comes from the loading matrix
                 # This is because the names in the parameter.specification table can be too long, and ends get cut off
                 which.observable=0
                 which.observable=which.observable+1
                 for(which.observable in 1:nrow(antecedent.fragment.Mplus.output$tech1$parameterSpecification$lambda)){

                   # threshold.parameter.table.as.data.frame <- as.data.frame(measurement.model.Mplus.output$tech1$parameterSpecification$tau)
                   name.of.observable <- rownames(antecedent.fragment.Mplus.output$tech1$parameterSpecification$lambda)[which.observable]

                   # first define an element of the parameter name in the matrix of antecedent draws
                   antecedent.parameter.name.in.antecedent.draws <- paste0(
                     name.of.observable, "$"
                   )

                   # select the draws for the thresholds for this observable
                   draws.for.thresholds <- select(unstandardized.draws.from.antecedent.fragments.as.data.frame, contains(antecedent.parameter.name.in.antecedent.draws))

                   # isolate the thresholds for this observable, with their parameter number
                   # do this by matching the observable name to the names of threshold parameters
                   # param.number.thresholds.this.observable <-
                   #   antecedent.fragment.Mplus.output$tech1$parameterSpecification$tau[ ,
                   #                                                                      grepl(name.of.observable, colnames(antecedent.fragment.Mplus.output$tech1$parameterSpecification$tau))
                   #   ]

                   # n.thresholds.this.observable <- length(param.number.thresholds.this.observable)
                   n.thresholds.this.observable <- ncol(draws.for.thresholds)

                   # Loop over thresholds
                   which.threshold=1
                   for(which.threshold in 1:n.thresholds.this.observable){

                     # Define the parameter name
                     antecedent.parameter.name <- paste0(
                       "[", name.of.observable, "$", which.threshold,"]"
                     )

                     # Proceed if the parameter name is new, don't if it already exists among list of antecedent parameters
                     if(sum(antecedent.parameter.name==antecedent.parameter.names)==0){

                       # Define the value of the entity
                       # Previous code that referred to the matrix via a column number
                       # Should work for taking values from the draws for the first antecedent
                       # But may not work in general
                       # value.to.paste <- draws.from.antecedent.fragments.as.data.frame[which.iter, param.number.thresholds.this.observable[which.threshold]]
                       # value.to.paste <- draws.from.antecedent.fragment.as.data.frame[which.iter, param.number.thresholds.this.observable[which.threshold]]

                       # this should work regardless of which antecedent

                       value.to.paste <- draws.for.thresholds[which.iter, which.threshold]


                       # If the value is equal to 0, perturb it a little
                       # Need to do this because if it's exactly equal to 0, Mplus will fix it properly
                       # but won't include it as one of the parameters saved in the Bayesian output
                       if(value.to.paste==0){

                         # Define a U(.0001, .0010) random variable to increment or decrement
                         epsilon <- runif(1, .0001, .0010)

                         # Add or subtract the epsilon based on a coin flip
                         should.add <- rbinom(1, 1, .5)==1
                         if(should.add) {value.to.paste = value.to.paste + epsilon }
                         if(!should.add) {value.to.paste = value.to.paste - epsilon }

                       } # closes if the value from the measurement model was exactly 0


                       # Define the line to append to syntax
                       # Usual way won't work b/c in Mplus the constraint is inside the brackets
                       # this.line <- paste0(
                       #   antecedent.parameter.name,
                       #   "@",
                       #   value.to.paste,
                       #   ";"
                       # )

                       this.line <- paste0(
                         "[", name.of.observable, "$", which.threshold, "@",
                         value.to.paste,
                         "];"
                       )

                       # Append the line
                       temp <- paste0(
                         temp, "\n",
                         this.line
                       )

                       # Add this parameter name to the collection of existing parameter names
                       antecedent.parameter.names <- c(antecedent.parameter.names, antecedent.parameter.name)

                     } # closes if the parameter is new to the list of antecedent parameters

                   } # closes loop over thresholds
                 } # closes loop over observables
               } # closes if there are fitted thresholds

             } # closes switch for adding code for thresholds

             # * * * * * Add code for structural coefficients ('beta' in Mplus)-------
             if(sum(parameters.to.condition.on == "structural coefficients") == 1){

               # Extract the parameter specfication from Mplus as a data frame
               parameter.specification <- as.data.frame(
                 antecedent.fragment.Mplus.output$tech1$parameterSpecification$beta
               )

               # Proceed if there are any fitted parameters
               if(any(parameter.specification > 0, na.rm = TRUE)){


                 # Go through each variable (row in the data frame)
                 which.variable = 1
                 for(which.variable in 1:nrow(parameter.specification)){

                   # Go through each other variable (column in the data frame)
                   which.other.variable = 1
                   for(which.other.variable in 1:ncol(parameter.specification)){

                     # Proceed if it's a a fitted parameter
                     if(!is.na(parameter.specification[which.variable, which.other.variable]) && parameter.specification[which.variable, which.other.variable] >0){
                       # if(parameter.specification[which.variable, which.other.variable] >0){

                       # Define the parameter name for the Mplus input code
                       antecedent.parameter.name <- paste0(
                         rownames(parameter.specification)[which.variable],
                         " on ",
                         colnames(parameter.specification)[which.other.variable]
                       )

                       # Proceed if the parameter name is new to the list of antecedent parameters, don't if it already exists among list of antecedent parameters
                       if(sum(antecedent.parameter.name==antecedent.parameter.names)==0){

                         # Define the value of the entity

                         # first define an element of the parameter name in the matrix of antecedent draws
                         antecedent.parameter.name.in.antecedent.draws <- paste0(
                           rownames(parameter.specification)[which.variable],
                           ".on.",
                           colnames(parameter.specification)[which.other.variable]
                         )

                         # Pull out the set of draws for that parameter
                         #select(draws.from.antecedent.fragment.as.data.frame, contains(antecedent.parameter.name.in.antecedent.draws))
                         # pull(draws.from.antecedent.fragment.as.data.frame, contains(antecedent.parameter.name.in.antecedent.draws))
                         # colnames(draws.from.antecedent.fragments.as.data.frame)

                         # colnames(unstandardized.draws.from.antecedent.fragments.as.data.frame)
                         # value.to.paste <- pull(
                         #   unstandardized.draws.from.antecedent.fragments.as.data.frame,
                         #   contains(antecedent.parameter.name.in.antecedent.draws)
                         # )[which.iter]


                         # The 'pull' function above sometimes returns an error
                         # says it's not a unique column
                         # I think it's because the name of the variable is a subset of others:
                         # E.g., 'Item1' is a subset of 'Item10' so it can't pick a unique column for the first
                         # Using the 'select' function below, which returns a table
                         # This seems to just return a table with 1 column, so it works
                         # Not sure why 'pull' indicates there are multiple columns, but 'select' ony yields one column
                         value.to.paste <- select(
                           unstandardized.draws.from.antecedent.fragments.as.data.frame,
                           contains(antecedent.parameter.name.in.antecedent.draws)
                         )[which.iter,1]

                         # If the value is equal to 0, perturb it a little
                         # Need to do this because if it's exactly equal to 0, Mplus will fix it properly
                         # but won't include it as one of the parameters saved in the Bayesian output
                         if(value.to.paste==0){

                           # Define a U(.0001, .0010) random variable to increment or decrement
                           epsilon <- runif(1, .0001, .0010)

                           # Add or subtract the epsilon based on a coin flip
                           should.add <- rbinom(1, 1, .5)==1
                           if(should.add) {value.to.paste = value.to.paste + epsilon }
                           if(!should.add) {value.to.paste = value.to.paste - epsilon }
                         } # closes if the value from the measurement model was exactly 0



                         # Define the line to append to syntax
                         this.line <- paste0(
                           # name.of.latent, " by ",
                           # rownames(measurement.model.Mplus.output$tech1$parameterSpecification$lambda)[which.observable],
                           #rownames(antecedent.fragment.Mplus.output$tech1$parameterSpecification$lambda)[which.observable],
                           antecedent.parameter.name,
                           " @",

                           # measurement.model.draws.as.data.frame[which.iter, measurement.model.Mplus.output$tech1$parameterSpecification$lambda[which.observable, which.latent]],
                           value.to.paste,
                           ";"

                         )

                         # Append the line
                         temp <- paste0(
                           temp, "\n",
                           this.line
                         )

                         # Add this parameter name to the collection of existing parameter names
                         antecedent.parameter.names <- c(antecedent.parameter.names, antecedent.parameter.name)


                         # closes if the parameter is new to the list of antecedent parameters

                       } # closes if this is a fitted parameter

                     } # closes loop over "other" variable (columns)

                   } # closes loop over each variable with fitted parameters (rows)

                 } # closes if there are fitted structural coefficients


               } # closes if there are any fitted parameters
             } # closes switch for adding code for structural coefficients

             # * * * * * Add code for (co)variances of latents ('psi' in Mplus)  -------
             if(sum(parameters.to.condition.on == "(co)variances of latent variables") == 1){

               # Extract the parameter specfication from Mplus as a data frame
               parameter.specification <- as.data.frame(
                 antecedent.fragment.Mplus.output$tech1$parameterSpecification$psi
               )

               # Proceed if there are any fitted parameters
               if(any(parameter.specification > 0, na.rm = TRUE)){

                 # Go through each variable (row in the data frame)
                 which.variable = 0
                 which.variable = which.variable + 1
                 for(which.variable in 1:nrow(parameter.specification)){

                   # Go through each other variable (column in the data frame)
                   which.other.variable = 0
                   which.other.variable = which.other.variable + 1
                   for(which.other.variable in 1:ncol(parameter.specification)){

                     # Proceed if it's a a fitted parameter

                     if(!is.na(parameter.specification[which.variable, which.other.variable]) && parameter.specification[which.variable, which.other.variable] >0){
                       # if(parameter.specification[which.variable, which.other.variable] >0){

                       # If it's a variance
                       if(rownames(parameter.specification)[which.variable] ==
                          colnames(parameter.specification)[which.other.variable]){

                         # Define the parameter name for the Mplus input code
                         antecedent.parameter.name <- paste0(
                           rownames(parameter.specification)[which.variable]
                         )


                         # Proceed if the parameter name is new to the list of antecedent parameters, don't if it already exists among list of antecedent parameters
                         if(sum(antecedent.parameter.name==antecedent.parameter.names)==0){

                           # Define the value of the entity

                           # first define an element of the parameter name in the matrix of antecedent draws
                           antecedent.parameter.name.in.antecedent.draws <- paste0(
                             rownames(parameter.specification)[which.variable]
                           )

                           # Now whittle down the antecendent draws to hopefully leave just one column
                           whittling <- select(
                             unstandardized.draws.from.antecedent.fragments.as.data.frame,
                             contains(antecedent.parameter.name.in.antecedent.draws)
                           )

                           whittling <- select(
                             whittling,
                             !contains("ON")
                           )

                           whittling <- select(
                             whittling,
                             !contains("BY")
                           )

                           whittling <- select(
                             whittling,
                             !contains("R-SQ")
                           )

                           whittling <- select(
                             whittling,
                             !contains("MEAN")
                           )

                           # paste the value if we've identified the isolated column
                           if(ncol(whittling) == 1) {
                             value.to.paste <- whittling[which.iter, 1]
                           } # closes if just one column in whittling


                           # If the value is equal to 0, perturb it a little
                           # Need to do this because if it's exactly equal to 0, Mplus will fix it properly
                           # but won't include it as one of the parameters saved in the Bayesian output
                           if(value.to.paste==0){

                             # Define a U(.0001, .0010) random variable to increment or decrement
                             epsilon <- runif(1, .0001, .0010)

                             # Add or subtract the epsilon based on a coin flip
                             should.add <- rbinom(1, 1, .5)==1
                             if(should.add) {value.to.paste = value.to.paste + epsilon }
                             if(!should.add) {value.to.paste = value.to.paste - epsilon }
                           } # closes if the value from the measurement model was exactly 0



                           # Define the line to append to syntax
                           this.line <- paste0(
                             antecedent.parameter.name,
                             " @",
                             value.to.paste,
                             ";"

                           )

                           # Append the line
                           temp <- paste0(
                             temp, "\n",
                             this.line
                           )

                           # Add this parameter name to the collection of existing parameter names
                           antecedent.parameter.names <- c(antecedent.parameter.names, antecedent.parameter.name)


                         } # closes if the parameter is new to the list of antecedent parameters


                       } # closes if it's a variance

                       # If it's a covariance
                       if(rownames(parameter.specification)[which.variable] !=
                          colnames(parameter.specification)[which.other.variable]){


                         # Define the parameter name for the Mplus input code
                         antecedent.parameter.name <- paste0(
                           colnames(parameter.specification)[which.variable],
                           " with ",
                           rownames(parameter.specification)[which.other.variable]

                         )

                         # Proceed if the parameter name is new to the list of antecedent parameters, don't if it already exists among list of antecedent parameters
                         if(sum(antecedent.parameter.name==antecedent.parameter.names)==0){

                           # Define the value of the entity

                           # first define an element of the parameter name in the matrix of antecedent draws
                           antecedent.parameter.name.in.antecedent.draws <- paste0(
                             rownames(parameter.specification)[which.variable],
                             ".with.",
                             colnames(parameter.specification)[which.other.variable]
                           )

                           # check if that parameter name is in the matrix of antecedent draws
                           # If it ISN'T, reverse the ordering of the variables in the name
                           # If it isn't see below
                           if(ncol(
                             select(
                               unstandardized.draws.from.antecedent.fragments.as.data.frame,
                               contains(antecedent.parameter.name.in.antecedent.draws)
                             )
                           ) == 0){
                             antecedent.parameter.name.in.antecedent.draws <- paste0(
                               colnames(parameter.specification)[which.other.variable],
                               ".with.",
                               rownames(parameter.specification)[which.variable]
                             )
                           } # closes if the parameter name is in the matrix of antecedent draws


                           # The 'pull' function above sometimes returns an error
                           # says it's not a unique column
                           # I think it's because the name of the variable is a subset of others:
                           # E.g., 'Item1' is a subset of 'Item10' so it can't pick a unique column for the first
                           # Using the 'select' function below, which returns a table
                           # This seems to just return a table with 1 column, so it works
                           # Not sure why 'pull' indicates there are multiple columns, but 'select' ony yields one column
                           value.to.paste <- select(
                             unstandardized.draws.from.antecedent.fragments.as.data.frame,
                             contains(antecedent.parameter.name.in.antecedent.draws)
                           )[which.iter,1]


                           # If the value is equal to 0, perturb it a little
                           # Need to do this because if it's exactly equal to 0, Mplus will fix it properly
                           # but won't include it as one of the parameters saved in the Bayesian output
                           if(value.to.paste==0){

                             # Define a U(.0001, .0010) random variable to increment or decrement
                             epsilon <- runif(1, .0001, .0010)

                             # Add or subtract the epsilon based on a coin flip
                             should.add <- rbinom(1, 1, .5)==1
                             if(should.add) {value.to.paste = value.to.paste + epsilon }
                             if(!should.add) {value.to.paste = value.to.paste - epsilon }
                           } # closes if the value from the measurement model was exactly 0

                           # Define the line to append to syntax
                           this.line <- paste0(
                             antecedent.parameter.name,
                             " @",
                             value.to.paste,
                             ";"

                           )



                           # Append the line
                           temp <- paste0(
                             temp, "\n",
                             this.line
                           )

                           # Add this parameter name to the collection of existing parameter names
                           antecedent.parameter.names <- c(antecedent.parameter.names, antecedent.parameter.name)





                         } # closes if the parameter is new to the list of antecedent parameters

                         # Since specifying a covariance as fixed, declare that there is a fixed covariance
                         # In Mplus, if there's a fixed covariance s.t. the covariance matrix cannot be partitioned into uncorrelated blocks,
                         # need to use a different sampler
                         fixed.covariance.from.antecedent <- TRUE

                       } # closes if it's a covariance



                     } # closes if this is a fitted parameter

                   } # closes loop over "other" variable (columns)

                 } # closes loop over each variable with fitted parameters (rows)

               } # closes if there are fitted parameters
             } # closes switch for adding code for (co)variances of latents

             # * * * * * Add code for latent means and intercepts ('alpha' in Mplus) -------
             if(sum(parameters.to.condition.on == "means and intercepts of latent variables") == 1){

             # Extract the parameter specfication from Mplus as a data frame
             parameter.specification <- as.data.frame(
               antecedent.fragment.Mplus.output$tech1$parameterSpecification$alpha
             )

             # Proceed if there are any fitted parameters
             if(any(parameter.specification > 0, na.rm = TRUE)){

               # Loop over observables
               which.variable=1
               for(which.variable in 1:ncol(parameter.specification)){

                 # proceed if this parameter is fitted (i.e., the value of the parameter specification table is > 0)
                 if(parameter.specification[1, which.variable] != 0){



                   name.of.variable <- colnames(parameter.specification)[which.variable]

                   # first define an element of the parameter name in the matrix of antecedent draws
                   antecedent.parameter.name.in.antecedent.draws <- paste0(
                     "MEAN.", name.of.variable
                   )

                   # select the draws for the intercept for this observable
                   draws.for.intercept <- select(unstandardized.draws.from.antecedent.fragments.as.data.frame, contains(antecedent.parameter.name.in.antecedent.draws))


                   # Define the parameter name
                   antecedent.parameter.name <- paste0(
                     "[", name.of.variable,"]"
                   )

                   # Proceed if the parameter name is new, don't if it already exists among list of antecedent parameters
                   if(sum(antecedent.parameter.name==antecedent.parameter.names)==0){

                     # Define the value of the entity
                     # Previous code that referred to the matrix via a column number
                     # Should work for taking values from the draws for the first antecedent
                     # But may not work in general
                     # value.to.paste <- draws.from.antecedent.fragments.as.data.frame[which.iter, param.number.thresholds.this.observable[which.threshold]]
                     # value.to.paste <- draws.from.antecedent.fragment.as.data.frame[which.iter, param.number.thresholds.this.observable[which.threshold]]

                     # this should work regardless of which antecedent

                     value.to.paste <- draws.for.intercept[which.iter, 1]


                     # If the value is equal to 0, perturb it a little
                     # Need to do this because if it's exactly equal to 0, Mplus will fix it properly
                     # but won't include it as one of the parameters saved in the Bayesian output
                     if(value.to.paste==0){

                       # Define a U(.0001, .0010) random variable to increment or decrement
                       epsilon <- runif(1, .0001, .0010)

                       # Add or subtract the epsilon based on a coin flip
                       should.add <- rbinom(1, 1, .5)==1
                       if(should.add) {value.to.paste = value.to.paste + epsilon }
                       if(!should.add) {value.to.paste = value.to.paste - epsilon }

                     } # closes if the value from the measurement model was exactly 0


                     # Define the line to append to syntax
                     # Usual way won't work b/c in Mplus the constraint is inside the brackets
                     # this.line <- paste0(
                     #   antecedent.parameter.name,
                     #   "@",
                     #   value.to.paste,
                     #   ";"
                     # )

                     this.line <- paste0(
                       "[", name.of.variable,"@",
                       value.to.paste,
                       "];"
                     )

                     # Append the line
                     temp <- paste0(
                       temp, "\n",
                       this.line
                     )

                     # Add this parameter name to the collection of existing parameter names
                     antecedent.parameter.names <- c(antecedent.parameter.names, antecedent.parameter.name)

                   } # closes if the parameter is new to the list of antecedent parameters
                 } # closes if this parameter is fitted (i.e., the value of the parameter specification table is > 0)


               } # closes loop over variable
             } # closes if there are fitted latent means and intercepts
             } # closes switch for adding code for latent means and intercepts


           } # closes loop over antecedent fragments

           # * * *  * ******************************************************************* ----
           # # * * * * Define the model.mplus syntax -------
           model.mplus.syntax <- temp

           # Clean up temp
           #rm(temp)


           # * * * * * Append syntax for this fragment -------
           model.mplus.syntax <- paste0(
             model.mplus.syntax,
             "\n",
             "\n",
             paste0(
               "! Parameters from fragment ",
               which.fragment
             ),
             # Mplus.MODEL.syntax.fragments[[which.fragment]]
             fragments[[which.fragment]]$model.syntax
           )


           # * * * * Define the syntax for the ANALYSIS portion of Mplus sytnax -----
           Mplus.ANALYSIS.syntax <-
             paste0("ESTIMATOR = BAYES;",  "\n",
                    paste0("CHAINS = 2;"), "\n",
                    paste0("POINT = MEAN;"), "\n",
                    paste0("PROCESSORS = 2;"),  "\n",
                    paste0("BSEED = ", sample.int(50000,1), ";")
             )


           # In Mplus, if there's a fixed covariance s.t. the covariance matrix cannot be partitioned into uncorrelated blocks,
           # need to use a different sampler
           if(fixed.covariance.from.antecedent == TRUE){
             Mplus.ANALYSIS.syntax <-
               paste0("ESTIMATOR = BAYES;",  "\n",
                      "ALGORITHM=GIBBS(RW);",  "\n",
                      paste0("CHAINS = 2;"), "\n",
                      paste0("POINT = MEAN;"), "\n",
                      paste0("PROCESSORS = 2;"),  "\n",
                      paste0("BSEED = ", sample.int(50000,1), ";")
               )
           }

           # * * * * Define the title for Mplus model -----
            Mplus.TITLE <- paste0(fragments[[which.fragment]]$name, ";")

           # * * * * Define the name of the BPARAMETER file saved from Mplus, which has the draws -----
            BPARAMETER.line <- paste0("BPARAMETER = ", fragments[[which.fragment]]$name, ".iter.", which.iter, ".draws.out;")


           # * * * * Define the mplus object for the fragment ------------

           # Two ways below
           # 1. If not estimating latent variables, will save the BPARAMETER file
           # 2. If  estimating latent variables, will not save the BPARAMETER file
           #    but will save the draws for the latent variables.
           #    And will change the BSEED for generating draws (plausible values in Mplus) for latent variables

           # Define the object if NOT estimating latent variables in this fragment
           if(fragments[[which.fragment]]$estimating.lvs == FALSE){

            fragment.mplus.object <- mplusObject(
              TITLE = Mplus.TITLE,
              VARIABLE = fragments[[which.fragment]]$variable.syntax,
              ANALYSIS = Mplus.ANALYSIS.syntax,
              MODEL = model.mplus.syntax,
              MODELPRIORS = fragments[[which.fragment]]$priors.syntax,
              SAVEDATA = BPARAMETER.line,


              OUTPUT = "
STANDARDIZED TECH1 TECH8;
",
              PLOT = "
TYPE = PLOT3;
",

              rdata = fragments[[which.fragment]]$data
            )
          } # closes object if NOT estimating latent variables in this fragment

           # Define the object if estimating latent variables in this fragment
           if(fragments[[which.fragment]]$estimating.lvs == TRUE){


             # * * * * Define the syntax for the ANALYSIS portion of Mplus sytnax -----
             # Mplus.ANALYSIS.syntax.previous <- Mplus.ANALYSIS.syntax
             # Mplus.ANALYSIS.syntax <- paste0(
             #   Mplus.ANALYSIS.syntax.previous,
             #   "\n",
             #   paste0("!Set seed for Markov chains"), "\n",
             #   paste0("BSEED = ", sample.int(50000,1), ";")
             # )
             # rm(Mplus.ANALYSIS.syntax.previous)

             # * * * * Define the syntax for the SAVEDATA portion of Mplus sytnax -----
             Mplus.SAVEDATA.syntax <-
               paste0("FILE = ", fragments[[which.fragment]]$name, " lvs", ".iter.", which.iter, ".dat;", "\n",
                      paste0("SAVE = FSCORES(20);"), "\n",
                      paste0("FACTORS = ", fragments[[which.fragment]]$lvs.to.estimate, ";")
               )

             # Define the Mplus object
             fragment.mplus.object <- mplusObject(
               TITLE = Mplus.TITLE,
               VARIABLE = fragments[[which.fragment]]$variable.syntax,
               ANALYSIS = Mplus.ANALYSIS.syntax,
               MODEL = model.mplus.syntax,
               MODELPRIORS = fragments[[which.fragment]]$priors.syntax,
               # SAVEDATA = BPARAMETER.line,
               SAVEDATA = Mplus.SAVEDATA.syntax,
               OUTPUT = "
STANDARDIZED TECH1 TECH8;
",
               PLOT = "
TYPE = PLOT3;
",

               # rdata = fragments[[which.fragment]]$data
               rdata = distinct.data # only using the distinct row response vectors
             )
           } # closes object if estimating latent variables in this fragment




           # * * * * Fit the model in Mplus and read in (some) results ------------

           fitted.model.mplus <- mplusModeler(
             object = fragment.mplus.object,

             # dataout = paste0(fragments.names[which.fragment], ".iter.", which.iter, " data.dat"),
             dataout = paste0(fragments[[which.fragment]]$name, ".iter.", which.iter, ".data.dat"),

             writeData = "always", # write out the data if it's not there
             hashfilename = FALSE, # do not add a hash to the data file name

             # modelout = paste0(fragments.names[which.fragment], ".iter.", which.iter, ".inp"),
             modelout = paste0(fragments[[which.fragment]]$name, ".iter.", which.iter, ".inp"),

             run = 1L
           )

           # Extract results if NOT estimating latent variables in this fragment
           if(fragments[[which.fragment]]$estimating.lvs == FALSE){

             # * * * * Extract the draws from the fitted  model -----
             # * * * * * Read in the Mplus output ---------

             # See comments above for case sensitivity in Windows and Linux
             # Mplus.output <- readModels(
             #   paste0(tolower(fragments[[which.fragment]]$name),
             #
             #          ".iter.",
             #          which.iter,
             #          ".out"))

             Mplus.output <- readModels(
               paste0(fragments[[which.fragment]]$name,

                      ".iter.",
                      which.iter,
                      ".out"))

             # Save the Mplus.output from the fragment as a separate object
             fragment.Mplus.output <- Mplus.output

             # Save the Mplus output object
             file.name <- paste0(fragments[[which.fragment]]$name, ".iter.", which.iter, ".Mplus output object.rds")
             saveRDS(
               object = fragment.Mplus.output,
               file = file.name
             )


             # Get the names of the file with MCMC draws
             Mplus.MCMC.file.name <- Mplus.output$savedata_info$bayesFile

             # Get the names of the columns of the file with MCMC draws
             # Mplus.output$savedata_info$bayesVarNames

             # * * * * * Read in the MCMC draws ------------
             Mplus.bparameters <- read_table(Mplus.MCMC.file.name, col_names = FALSE)
             colnames(Mplus.bparameters) <- Mplus.output$savedata_info$bayesVarNames


             # * * * * * Extract a random row (after convergence) as the saved draw -------

             # define a temp object
             temp1 <- Mplus.bparameters

             # define a temp object for each of the 2 chains
             # Stores the draws from the second halves of each chain
             temp2.1 <- temp1[(1+1*nrow(temp1)/4):(nrow(temp1)/2) ,]
             temp2.2 <- temp1[(1+3*nrow(temp1)/4):nrow(temp1) ,]

             # Put those second halves together
             temp2 <- rbind(temp2.1, temp2.2)

             # Take a sample from those second halves, to be returned as the draw
             temp3 <- temp2[sample(nrow(temp2), 1), ]
             # temp3
           } # closes if NOT estimating latent variables in this fragment

           # Extract results if estimating latent variables in this fragment
           if(fragments[[which.fragment]]$estimating.lvs == TRUE){

             # * * * * Extract the draws from the fitted  model -----
             # * * * * * Read in the Mplus output ---------

             # See comments above on case sensitivity in Linux and Window
             # Mplus.output <- readModels(
             #   paste0(tolower(fragments[[which.fragment]]$name),
             #
             #          ".iter.",
             #          which.iter,
             #          ".out"))

             Mplus.output <- readModels(
               paste0(fragments[[which.fragment]]$name,

                      ".iter.",
                      which.iter,
                      ".out"))

             # Save the Mplus.output from the fragment as a separate object
             fragment.Mplus.output <- Mplus.output

             # Save the Mplus output object
             file.name <- paste0(fragments[[which.fragment]]$name, ".iter.", which.iter, ".Mplus output object.rds")
             saveRDS(
               object = fragment.Mplus.output,
               file = file.name
             )



             # Read in the savedata from Mplus
             # See comments above about case sensitivity in Linux and Windows
             # savedat <- readModels(paste0(tolower(fragments[[which.fragment]]$name),
             #
             #                              ".iter.",
             #                              which.iter,
             #                              ".out"),
             #                       what = "savedata")


             savedat <- readModels(paste0(fragments[[which.fragment]]$name,

                                          ".iter.",
                                          which.iter,
                                          ".out"),
                                   what = "savedata")


             # Extract a draw for the latent variables
             # Effectively looks at what Mplus generates as the savedata file
             # Then looks for the name(s) of the latent variable(s)
             # Then takes the 5th draw for that value
             temp2 <- savedat$savedata %>%
               dplyr::select(contains(fragments[[which.fragment]]$lvs.to.estimate)) %>%
               dplyr::select(contains("020"))


             # Convert the results to a row matrix with 1 column per person
             temp3 <- t(temp2)
             # temp3


           } # closes if estimating latent variables in this fragment

          # Return the result
           #print(paste("iter", which.iter, "--end"))
          return(temp3)

           # clean up some files here, if desired
           # if(!retain.iteration.files){
           #
           #   # For now, keep all from first iteration, delete rest
           #   if(which.iter != 1){
           #
           #     # file names to delete
           #     to_be_deleted <- list.files(getwd(), pattern = paste0(fragments[[which.fragment]]$name, ".iter.", which.iter), ignore.case = TRUE)
           #     # file.paths.to.be.deleted <- file.path(getwd(), to_be_deleted)
           #     # file.remove(to_be_deleted)
           #     # file.remove(file.paths.to.be.deleted)
           #
           #     # Try a more robust file removal approach
           #     for (file in to_be_deleted) {
           #       if (file.exists(file)) {
           #         removal_result <- try(file.remove(file), silent = TRUE)
           #         if (inherits(removal_result, "try-error") || !removal_result) {
           #           warning(paste("Failed to remove file:", file))
           #         }
           #       }
           #     }
           #
           #
           #   } # closes if which.iter is not equal to 1
           #
           # }


         } # closes foreach

          # Convert list to dataframe
          #draws.from.MUPPET.model.this.batch <- as.data.frame(data.table::rbindlist(draws.from.MUPPET.model.this.batch.list))

          # Append draws from this batch
          draws.from.MUPPET.model.this.fragment <- rbind(draws.from.MUPPET.model.this.fragment, draws.from.MUPPET.model.this.batch)

          current.time <- Sys.time()
          print(paste0("Cleaning up files for Fragment ", which.fragment, " batch ", which.estimation.batch, " of ", n.estimation.batches,  " ", round(current.time, units="secs")))

          sink("R Session Info.out", append = TRUE)
          print(paste0("Cleaning up files for Fragment ", which.fragment, " batch ", which.estimation.batch, " of ", n.estimation.batches,  " ", round(current.time, units="secs")))
          sink()

          # Save the draws from all the batches so far
          # A temp file, to be used in case of need of troubleshooting
          # Will remove later

          # convert to data frame if it's not
          if(!is.data.frame(draws.from.MUPPET.model.this.fragment)) {
            draws.from.MUPPET.model.this.fragment <- as.data.frame(draws.from.MUPPET.model.this.fragment)
          }


          write_csv(
            draws.from.MUPPET.model.this.fragment,
            paste0(fragments[[which.fragment]]$name, " draws up through batch ", which.estimation.batch, ".csv")
          )

          # Delete file with draws up through previous branch
          if(which.estimation.batch > 1){
            file.remove(paste0(fragments[[which.fragment]]$name, " draws up through batch ", which.estimation.batch-1, ".csv"))
          }


          # * * * Copy and rename the file with Mplus output object from the first iteration ----
          # To be a "shell"
          # Needed if this fragment is antecedent to other fragments
          if(which.estimation.batch == 1){
            original.file.name <-  paste0(
              fragments[[which.fragment]]$name, ".iter.1.Mplus output object.rds"
            )
            new.file.name <- paste0(
              fragments[[which.fragment]]$name, ".shell Mplus output object.rds"
            )

            file.copy(original.file.name, new.file.name)

          }


          # * * * Also copy files for iteration 1 to the fragment folder Copy and rename the file with Mplus output object from the first iteration ----
          if(which.estimation.batch == 1){

            # The .inp file
            original.file.name <-  paste0(
              fragments[[which.fragment]]$name, ".iter.1.inp"
            )
            new.file.name <- paste0(
              fragment.folder,
              fragments[[which.fragment]]$name, ".iter.1.inp"
            )
            file.copy(original.file.name, new.file.name)

            # The data.dat file
            original.file.name <-  paste0(
              fragments[[which.fragment]]$name, ".iter.1.data.dat"
            )
            new.file.name <- paste0(
              fragment.folder,
              fragments[[which.fragment]]$name, ".iter.1.data.dat"
            )
            file.copy(original.file.name, new.file.name)

            # The .out file
            original.file.name <-  paste0(
              fragments[[which.fragment]]$name, ".iter.1.out"
            )
            new.file.name <- paste0(
              fragment.folder,
              fragments[[which.fragment]]$name, ".iter.1.out"
            )
            file.copy(original.file.name, new.file.name)

            # The draws.out file
            original.file.name <-  paste0(
              fragments[[which.fragment]]$name, ".iter.1.draws.out"
            )
            new.file.name <- paste0(
              fragment.folder,
              fragments[[which.fragment]]$name, ".iter.1.draws.out"
            )
            file.copy(original.file.name, new.file.name)


            # The .gh5 file
            original.file.name <-  paste0(
              fragments[[which.fragment]]$name, ".iter.1.gh5"
            )
            new.file.name <- paste0(
              fragment.folder,
              fragments[[which.fragment]]$name, ".iter.1.gh5"
            )
            file.copy(original.file.name, new.file.name)


            # The .rds file
            original.file.name <-  paste0(
              fragments[[which.fragment]]$name, ".iter.1.Mplus output object.rds"
            )
            new.file.name <- paste0(
              fragment.folder,
              fragments[[which.fragment]]$name, ".iter.1.Mplus output object.rds"
            )
            file.copy(original.file.name, new.file.name)

          }



          # * * * Delete iteration-specific files if desired -----
          if(!retain.iteration.files){


            # make a copy of the first iteration, to be saved
            file.name <- list.files(getwd(), pattern = ".iter.1.inp")

            # only keep the file for this fragment
            file.name <- file.name[grep(fragments[[which.fragment]]$name, file.name)]

            # copy the file
            file.copy(file.name, "temp.inp")

            # delete files
            to_be_deleted <- list.files(getwd(), pattern = ".iter")
            file.remove(to_be_deleted)

            # rename copy of first iteration
            file.rename("temp.inp", file.name)
          }


          # tidy up
          rm(draws.from.MUPPET.model.this.batch)

          current.time <- Sys.time()
          print(paste0("Completed Fragment ", which.fragment, " batch ", which.estimation.batch, " of ", n.estimation.batches, " ", round(current.time, units="secs")))

          sink("R Session Info.out", append = TRUE)
          print(paste0("Completed Fragment ", which.fragment, " batch ", which.estimation.batch, " of ", n.estimation.batches, " ", round(current.time, units="secs")))
          sink()

        } # closes loop over estimation batches


        # * * * Close parallel environment ------
        stopImplicitCluster()
        current.time <- Sys.time()
        print(paste0("Sampling for Fragment ", which.fragment, " completed ", round(current.time, units="secs")))

        sink("R Session Info.out", append = TRUE)
        print(paste0("Sampling for Fragment ", which.fragment, " completed ", round(current.time, units="secs")))
        sink()


        # * * * ******************************************************************* ----
        # * * * Assemble the draws from the MUPPET model if not estimating latent variables from ----
        # draws for this fragment
        # unstandardized draws from previous fragments not present in the current fragment
        # older code using 'matches' is commented out below

        if(fragments[[which.fragment]]$estimating.lvs == FALSE){

          draws.from.MUPPET.model <- add_column(
            draws.from.MUPPET.model.this.fragment,

            !!!unstandardized.draws.from.antecedent.fragments.as.data.frame[
              setdiff(
                names(unstandardized.draws.from.antecedent.fragments.as.data.frame), names(draws.from.MUPPET.model.this.fragment)
              )
            ],
            .before = 1
          )
          # colnames(draws.from.MUPPET.model)
          # colnames(draws.all.parameters.Mplus)

          # * * * Save the draws from MUPPET for this fragment -----
          if(save.draws.from.MUPPET) write_csv(
            # draws.all.parameters.Mplus,
            draws.from.MUPPET.model,

            # paste0(fragments.names[which.fragment], " draws.csv")
            paste0(fragments[[which.fragment]]$name, " draws.csv")
          )

          # Delete the last temp file with draws from doing things in batches
          file.name <- paste0(fragments[[which.fragment]]$name, " draws up through batch ", n.estimation.batches, ".csv")
          if(file.exists(file.name)){
            file.remove(file.name)
          }


          # * * * Define a folder for storing output
          current.time <- Sys.time()
          print(paste0("Posterior summary for Fragment ", which.fragment, " initiated ", round(current.time, units="secs")))

          sink("R Session Info.out", append = TRUE)
          print(paste0("Posterior summary for Fragment ", which.fragment, " initiated ", round(current.time, units="secs")))
          sink()

          posterior.summary.folder <- paste0(fragment.folder, "/Summary of Posterior/")
          if(!dir.exists(posterior.summary.folder)) dir.create(posterior.summary.folder)
          setwd(posterior.summary.folder)

          # Compute the summary statistics for the draws
          draws.to.analyze <- as.mcmc(draws.from.MUPPET.model)
          # dim(draws.to.analyze)
          # dim(na.omit(draws.from.MUPPET.model[, 66:75]))
          summary.statistics.MUPPET <- MCMCsummary(
            draws.to.analyze,
            HPD=TRUE,
            n.eff=FALSE,
            Rhat=FALSE,
            round=8,
            func=median,
            func_name = "median"
          )

          # * * * Write out the summary statistics for the MUPPET model ------
          if(save.summary.stats.from.MUPPET){
            file.name=paste0(
              # fragments.names[which.fragment],
              fragments[[which.fragment]]$name,
              " summary statistics",
              ".csv"
            )
            write.csv(
              x=summary.statistics.MUPPET,
              file=file.name
            )
          } # closes if saving summary statistics from MUPPET model

          # * * * Write out the summary statistics for the MUPPET model in Word ------
          if(save.summary.stats.from.MUPPET.in.Word){

            # Function to remove text up until and including a character
            remove_prefix_function <- function(
    text_vector,
    delimiter # will replace up to (and including) first instance of this
            ) {
              sub(paste0("^.*?", delimiter), "", text_vector)
            }

            # Create folder for output
            # fragment.folder <- paste0(getwd(), "/Fragment ", which.fragment, " ", fragments[[which.fragment]]$name, "/")
            # if(!dir.exists(fragment.folder)) dir.create(fragment.folder)


            # Prepare the summary statistics table for publishing ------
            # digits.to.round.for.Word = 2

            # which.fragment = 3
            # fragments[[which.fragment]]$name

            # * Read in the summary statistics table ----
            # file.name <- paste0(fragments[[which.fragment]]$name, " summary statistics.csv")
            # raw.table <-  read_csv(
            #   file = file.name
            # )
            raw.table <- summary.statistics.MUPPET

            # Set the working table
            working.table <- raw.table
            working.table <- rownames_to_column(working.table, var = "parameter")
            # colnames(working.table)[1] = "parameter"

            # * Discarding chain, iterations, and R-squared rows -----
            working.table <- working.table %>%
              filter(
                !grepl("chain", parameter, ignore.case = TRUE)
              )  %>%
              filter(
                !grepl("iteration", parameter, ignore.case = TRUE)
              )   %>%
              filter(
                !grepl("R-SQUARE", parameter, ignore.case = TRUE)
              )



            # * Round to selected digits -----
            working.table <- working.table %>%
              mutate_if(is.numeric, round, digits.to.round.for.Word)

            # * Create interval in single column -----
            working.table$`95% HPDI` <-
              paste0(
                "(",
                str_trim(format(working.table$`95%_HPDL`, nmall=digits.to.round.for.Word)), # remove whitespace
                ", ",
                str_trim(format(working.table$`95%_HPDU`, nmall=digits.to.round.for.Word)), # remove whitespace
                ")"
              )


            # * Drop columns not needed -----
            working.table <- select(working.table, -median)
            working.table <- select(working.table, -`95%_HPDL`)
            working.table <- select(working.table, -`95%_HPDU`)

            # * Capitalize columns -----
            colnames(working.table) <- capitalize(colnames(working.table))


            # * Create table ----
            # working.table <- tt(
            #   working.table,
            #   notes = list(
            #     a = list(i = 0, j=ncol(working.table), text = "Highest posterior density interval")
            #   )
            # )


            # * Store the table with all parameters -----
            working.table.all.parameters <- working.table

            # Write out the table
            table.to.write <- working.table.all.parameters

            if(nrow(table.to.write) >0){
              # file.name <- paste0(fragments[[which.fragment]]$name, " summary statistics.docx")
              # tt(table.to.write) |> save_tt(file.name, overwrite=TRUE)

              file.name <- paste0(fragments[[which.fragment]]$name, " summary statistics reduced.csv")
              write.csv(table.to.write, paste0(file.name), row.names = FALSE, quote = TRUE)

            }


            # * * unstandardized solution -----
            if(1==1){
              working.table <- working.table.all.parameters %>%
                filter(
                  !grepl("STD", Parameter, ignore.case = TRUE)
                )

              # * * Remove text before a character from the parameter names -----
              working.table$Parameter <- remove_prefix_function(
                text_vector = working.table$Parameter,
                delimiter = "_" # will remove text before (and including) the first instance of this character
              )


              # * * * loadings -----
              loadings.table <- working.table %>%
                filter(
                  grepl("\\.BY\\.", Parameter, ignore.case = TRUE)
                )

              # Write out the table
              table.to.write <- loadings.table

              # sort it
              table.to.write <- dplyr::arrange(table.to.write, Parameter)

              # Write out the table
              if(nrow(table.to.write) >0){
                file.name <- paste0(fragments[[which.fragment]]$name, " unstandardized loadings summary statistics.docx")
                tt(table.to.write) |> save_tt(paste0(file.name), overwrite=TRUE)


              }


              # * * * intercepts and means -----
              intercepts.and.means.table <- working.table %>%
                filter(
                  grepl("MEAN\\.", Parameter, ignore.case = TRUE)
                ) %>%
                filter(
                  !grepl("\\$", Parameter, ignore.case = TRUE)

                )

              intercepts.and.means.table$Parameter <- remove_prefix_function(
                text_vector = intercepts.and.means.table$Parameter,
                delimiter = "MEAN\\." # will remove text before (and including) the first instance of this character
              )


              # Write out the table
              table.to.write <- intercepts.and.means.table

              # sort it
              # table.to.write <- dplyr::arrange(table.to.write, Parameter)

              if(nrow(table.to.write) >0){
                file.name <- paste0(fragments[[which.fragment]]$name, " unstandardized intercepts and means summary statistics.docx")
                tt(table.to.write) |> save_tt(paste0(file.name), overwrite=TRUE)

              }

              # * * * thresholds -----
              thresholds.table <- working.table %>%
                filter(
                  grepl("MEAN\\.", Parameter, ignore.case = TRUE)
                ) %>%
                filter(
                  grepl("\\$", Parameter, ignore.case = TRUE)
                )

              # Write out the table
              table.to.write <- thresholds.table

              # sort it
              # table.to.write <- dplyr::arrange(table.to.write, Parameter)

              if(nrow(table.to.write) >0){
                file.name <- paste0(fragments[[which.fragment]]$name, " unstandardized thresholds summary statistics.docx")
                tt(table.to.write) |> save_tt(paste0(file.name), overwrite=TRUE)

              }





              # * * * covariances -----
              covariances.table <- working.table %>%
                filter(
                  grepl("\\.WITH\\.", Parameter, ignore.case = TRUE)
                )

              # Write out the table
              table.to.write <- covariances.table

              # sort it
              table.to.write <- dplyr::arrange(table.to.write, Parameter)

              if(nrow(table.to.write) >0){
                file.name <- paste0(fragments[[which.fragment]]$name, " unstandardized covariances summary statistics.docx")
                tt(table.to.write) |> save_tt(paste0(file.name), overwrite=TRUE)

              }




              # * * *  structural coefficents -----
              coefficients.table <- working.table %>%
                filter(
                  grepl("\\.ON\\.", Parameter, ignore.case = TRUE)
                )

              # Write out the table
              table.to.write <- coefficients.table

              # sort it
              table.to.write <- dplyr::arrange(table.to.write, Parameter)

              if(nrow(table.to.write) >0){
                file.name <- paste0(fragments[[which.fragment]]$name, " unstandardized coefficients summary statistics.docx")
                tt(table.to.write) |> save_tt(paste0(file.name), overwrite=TRUE)

              }



              # * * *  variances and other parameters  -----
              variances.and.others.table <- working.table %>%
                filter(
                  !grepl("\\.BY\\.", Parameter, ignore.case = TRUE)
                ) %>%
                filter(
                  !grepl("MEAN\\.", Parameter, ignore.case = TRUE)
                ) %>%
                filter(
                  !grepl("\\.WITH\\.", Parameter, ignore.case = TRUE)
                ) %>%
                filter(
                  !grepl(".ON\\.", Parameter, ignore.case = TRUE)
                )

              # Write out the table
              table.to.write <- variances.and.others.table

              # sort it
              # table.to.write <- dplyr::arrange(table.to.write, Parameter)

              if(nrow(table.to.write) >0){
                file.name <- paste0(fragments[[which.fragment]]$name, " unstandardized variances and remaining parameters summary statistics.docx")
                tt(table.to.write) |> save_tt(paste0(file.name), overwrite=TRUE)

              }

            } # closes switch for unstandardized solution



            # * * STDYX solution -----
            if(1==1){
              working.table <- working.table.all.parameters %>%
                filter(
                  grepl("STDYX", Parameter, ignore.case = TRUE)
                )

              # * * Remove text before a character from the parameter names -----
              working.table$Parameter <- remove_prefix_function(
                text_vector = working.table$Parameter,
                delimiter = "_" # will remove text before (and including) the first instance of this character
              )


              # * * * loadings -----
              loadings.table <- working.table %>%
                filter(
                  grepl("\\.BY\\.", Parameter, ignore.case = TRUE)
                )

              # Write out the table
              table.to.write <- loadings.table

              # sort it
              table.to.write <- dplyr::arrange(table.to.write, Parameter)

              # Write out the table
              if(nrow(table.to.write) >0){
                file.name <- paste0(fragments[[which.fragment]]$name, " standardized loadings summary statistics.docx")
                tt(table.to.write) |> save_tt(paste0(file.name), overwrite=TRUE)

              }


              # * * * intercepts and means -----
              intercepts.and.means.table <- working.table %>%
                filter(
                  grepl("MEAN\\.", Parameter, ignore.case = TRUE)
                ) %>%
                filter(
                  !grepl("\\$", Parameter, ignore.case = TRUE)

                )

              intercepts.and.means.table$Parameter <- remove_prefix_function(
                text_vector = intercepts.and.means.table$Parameter,
                delimiter = "MEAN\\." # will remove text before (and including) the first instance of this character
              )


              # Write out the table
              table.to.write <- intercepts.and.means.table

              # sort it
              # table.to.write <- dplyr::arrange(table.to.write, Parameter)

              if(nrow(table.to.write) >0){
                file.name <- paste0(fragments[[which.fragment]]$name, " standardized intercepts and means summary statistics.docx")
                tt(table.to.write) |> save_tt(paste0(file.name), overwrite=TRUE)

              }

              # * * * thresholds -----
              thresholds.table <- working.table %>%
                filter(
                  grepl("MEAN\\.", Parameter, ignore.case = TRUE)
                ) %>%
                filter(
                  grepl("\\$", Parameter, ignore.case = TRUE)
                )

              # Write out the table
              table.to.write <- thresholds.table

              # sort it
              # table.to.write <- dplyr::arrange(table.to.write, Parameter)

              if(nrow(table.to.write) >0){
                file.name <- paste0(fragments[[which.fragment]]$name, " standardized thresholds summary statistics.docx")
                tt(table.to.write) |> save_tt(paste0(file.name), overwrite=TRUE)

              }





              # * * * covariances -----
              covariances.table <- working.table %>%
                filter(
                  grepl("\\.WITH\\.", Parameter, ignore.case = TRUE)
                )

              # Write out the table
              table.to.write <- covariances.table

              # sort it
              table.to.write <- dplyr::arrange(table.to.write, Parameter)

              if(nrow(table.to.write) >0){
                file.name <- paste0(fragments[[which.fragment]]$name, " standardized covariances summary statistics.docx")
                tt(table.to.write) |> save_tt(paste0(file.name), overwrite=TRUE)

              }




              # * * *  structural coefficents -----
              coefficients.table <- working.table %>%
                filter(
                  grepl("\\.ON\\.", Parameter, ignore.case = TRUE)
                )

              # Write out the table
              table.to.write <- coefficients.table

              # sort it
              table.to.write <- dplyr::arrange(table.to.write, Parameter)

              if(nrow(table.to.write) >0){
                file.name <- paste0(fragments[[which.fragment]]$name, " standardized coefficients summary statistics.docx")
                tt(table.to.write) |> save_tt(paste0(file.name), overwrite=TRUE)

              }



              # * * *  variances and other parameters  -----
              variances.and.others.table <- working.table %>%
                filter(
                  !grepl("\\.BY\\.", Parameter, ignore.case = TRUE)
                ) %>%
                filter(
                  !grepl("MEAN\\.", Parameter, ignore.case = TRUE)
                ) %>%
                filter(
                  !grepl("\\.WITH\\.", Parameter, ignore.case = TRUE)
                ) %>%
                filter(
                  !grepl(".ON\\.", Parameter, ignore.case = TRUE)
                )

              # Write out the table
              table.to.write <- variances.and.others.table

              # sort it
              # table.to.write <- dplyr::arrange(table.to.write, Parameter)

              if(nrow(table.to.write) >0){
                file.name <- paste0(fragments[[which.fragment]]$name, " standardized variances and remaining parameters summary statistics.docx")
                tt(table.to.write) |> save_tt(paste0(file.name), overwrite=TRUE)

              }

            } # closes switch for STDYX solution
          } # closes if saving summary statistics in Word

          # * * * Write out the summary plots for the MUPPET model ------
          if(save.summary.plots.from.MUPPET != "none"){


            # Create a folder for plots
            plots.folder <- paste0(posterior.summary.folder, "/Plots/")
            if(!dir.exists(plots.folder)) dir.create(plots.folder)
            setwd(plots.folder)

            # Append/Replace columns indicating the chain and iteration number
            draws.to.analyze <- draws.from.MUPPET.model
            draws.to.analyze$Chain.number <- NA
            draws.to.analyze$Iteration.number <- NA

            for (which.chain in 1:n.chains) {
              start_row <- (which.chain-1) * n.iters.per.chain.after.warmup.and.burnin + 1
              end_row <- which.chain * n.iters.per.chain.after.warmup.and.burnin
              draws.to.analyze$Chain.number[start_row:end_row] = which.chain
              draws.to.analyze$Iteration.number[start_row:end_row] = seq(1:n.iters.per.chain.after.warmup.and.burnin)

            }

            # Select only unstandardized values if desired
            if(save.summary.plots.from.MUPPET == "unstandardized"){
              draws = draws.to.analyze

              draws <- draws %>%
                select(
                  !contains("STD")
                )
            } # closes selecting only unstandardized values if desired

            # Need to convert to an MCMC.list
            draws.to.analyze <- df_to_mcmc_list(draws, chain_col = "Chain.number", iter_col = "Iteration.number")

            # Tidy up
            rm(draws)


            # Plot the results
            # trace, densities, and autocorrelations
            # mcmcplots::mcmcplot(
            mcmcplot.mod(
              mcmcout = draws.to.analyze,
              #parms = parameters.for.convergence,
              dir = plots.folder,
              style="plain",
              filename = "MCMC Plots"
            )

          } # closes if saving summary plots from MUPPET model
          current.time <- Sys.time()
          print(paste0("Posterior summary for Fragment ", which.fragment, " completed ", round(current.time, units="secs")))

          sink("R Session Info.out", append = TRUE)
          print(paste0("Posterior summary for Fragment ", which.fragment, " completed ", round(current.time, units="secs")))
          sink()


          # Move back to home folder
          setwd(home.folder)


        } # closes if NOT estimating latent variables in this fragment


        # * * * Assemble the draws from the MUPPET model if estimating latent variables ----
        if(fragments[[which.fragment]]$estimating.lvs == TRUE){

          # Store the draws from the MUPPET model as being *per response vector*
          draws.from.MUPPET.model.this.fragment.per.response.vector <- draws.from.MUPPET.model.this.fragment

          # remove draws from MUPPET model for this fragment object
          # Will be reconstituted later per each case in the dataset
          rm(draws.from.MUPPET.model.this.fragment)

          # Append the iteration number as a column
          # Iteration.number <- seq(1,nrow(draws.from.MUPPET.model.this.fragment.per.response.vector))
          # draws.from.MUPPET.model.this.fragment.per.response.vector <- cbind(draws.from.MUPPET.model.this.fragment.per.response.vector, Iteration.number)

          # Remove rownames
          # row.names(draws.from.MUPPET.model.this.fragment.per.response.vector) <- NULL


          # function to find matching row in a data frame
          find_matching_row2 <- function(main_df, other_df, row_index) {
            # Extract the target row from the other data frame
            target_row <- other_df[row_index, ]

            # Check if the column names match between the data frames
            if(!all(names(target_row) == names(main_df))) {
              warning("Column names don't match between data frames")
              return(NA)
            }

            # Look for matching rows
            matches <- which(apply(main_df, 1, function(row) {
              all(row == target_row)
            }))

            return(matches)
          }


          # Create object to store the draws *per case*
          draws.from.MUPPET.model.this.fragment <- matrix(NA, nrow=nrow(draws.from.MUPPET.model.this.fragment.per.response.vector), ncol=nrow(fragments[[which.fragment]]$data))


          # Loop over people
          i=0
          i=i+1
          for(i in 1:nrow(fragments[[which.fragment]]$data)){


            # Find the row in the matrix of distinct data (i.e., unique response patterns)
            # that matches this person's data
            which.response.pattern <- find_matching_row2(
              main_df = distinct.data,
              other_df = fragments[[which.fragment]]$data,
              row_index = i)


            # Store the draws for this person
            # Plucks out appropriate column of draws pertaining to the correct distinct response pattern
            draws.from.MUPPET.model.this.fragment[,i] <- draws.from.MUPPET.model.this.fragment.per.response.vector[,which.response.pattern]
          }


          # Append the iteration number as a column
          Iteration.number <- seq(1,nrow(draws.from.MUPPET.model.this.fragment))
          draws.from.MUPPET.model.this.fragment <- cbind(draws.from.MUPPET.model.this.fragment, Iteration.number)


          # Convert to a dataframe
          draws.from.MUPPET.model.this.fragment <- as.data.frame(draws.from.MUPPET.model.this.fragment)

          # * * * Save the draws from MUPPET for just this fragment -----
          if(save.draws.from.MUPPET) write_csv(
            draws.from.MUPPET.model.this.fragment,
            paste0(fragments[[which.fragment]]$name, " draws.csv")
          )

          # * * * Define a folder for storing output
          current.time <- Sys.time()
          print(paste0("Posterior summary for Fragment ", which.fragment, " initiated ", round(current.time, units="secs")))

          posterior.summary.folder <- paste0(fragment.folder, "/Summary of Posterior/")
          if(!dir.exists(posterior.summary.folder)) dir.create(posterior.summary.folder)
          setwd(posterior.summary.folder)

          # Compute the summary statistics for the draws
          draws.to.analyze <- as.mcmc(draws.from.MUPPET.model.this.fragment)
          summary.statistics.MUPPET <- MCMCsummary(
            draws.to.analyze,
            HPD=TRUE,
            n.eff=FALSE,
            Rhat=FALSE,
            round=8,
            func=median,
            func_name = "median"
          )

          # * * * Write out the summary statistics for the MUPPET model ------
          if(save.summary.stats.from.MUPPET){
            file.name=paste0(
              fragments[[which.fragment]]$name,
              " summary statistics",
              ".csv"
            )
            write.csv(
              x=summary.statistics.MUPPET,
              file=file.name
            )
          } # closes if saving summary statistics from MUPPET model
          print(paste0("Posterior summary for Fragment ", which.fragment, " completed"))


          # Move back to home folder
          setwd(home.folder)


        } # closes if estimating latent variables in this fragment


        } # closes if the fragment is conditional on other fragments

      } # closes if this fragment is to be fit


        # * Move fragment specific files to fragment folder -----

        # Select files with the fragment name in them
        files_to_copy <- list.files(getwd(), pattern = paste0(fragments[[which.fragment]]$name))

        # Select files with a . in the filename, to avoid selecting the subfolder
        files_to_copy <- grep("\\.", files_to_copy, value=TRUE)

        # Move them
        file.copy(files_to_copy, file.path(fragment.folder, files_to_copy))
        rm(files_to_copy)


        # * End the timer -------
        end.time <- Sys.time()
        time.taken <- round(difftime(end.time, start.time, units='mins'), 2)
        time.taken
        print(paste0("Fragment ", which.fragment, " took ", time.taken, " minutes"))

        # * Sink out file some session info -----
        sink("R Session Info.out", append = TRUE)
        print(paste0("Fragment ", which.fragment, " took ", time.taken, " minutes"))
        sink()

      } # closes loop over fragments (starting at 1)


     #} # closes if software.environment == Mplus

  #} # closes if modular == TRUE

  # * End the timer -------
  end.time <- Sys.time()
  time.taken <- round(difftime(end.time, start.time.fragment.1, units='mins'), 2)
  time.taken
  print(paste0("MUPPET took ", time.taken, " minutes"))

  # * Sink out file some session info -----
  sink("R Session Info.out", append = TRUE)
  print(paste0("MUPPET took ", time.taken, " minutes"))
  sink()

  # Define what is returned
  to.return <- paste0("MUPPET took ", time.taken, " minutes")
  return(to.return)

} # closes MUPPET.modular.function

