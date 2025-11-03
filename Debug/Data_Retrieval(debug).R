############################ 1. General Information ############################

# This file contains the updated code for replicating exact datasets from completed 
# SEM simulation studies by restoring the original RNG state used during data 
# generation. Updated to reflect new methods (upi, lsam) and VITA-based generation.

#' @param condition_id Integer. Identifier for the experimental condition to replicate. Must be within range of available conditions in results file.
#' @param rep_id Integer. Replication number within the condition to reproduce. Must be within range of completed replications for that condition.
#' @param study Integer. Study identifier (1 or 2). Determines which models and default results file to use. Default is 1.
#' @param results_file Character. Path to results file containing simulation data. If NULL, uses study-specific default path. Default is NULL.
#' @param base_path Character. Base directory path for simulation files. Default is ".".
#' @return A data.frame containing the replicated dataset with attributes:
#'   - condition: The experimental condition parameters
#'   - condition_id: Integer identifier for the condition
#'   - rep_id: Replication number within the condition
#'   - study: Study identifier (1 or 2)
#'   - observed_R2: Observed R-squared values for endogenous variables
#'   - observed_reliabilities: Observed reliability values for latent variables

############################### 2. Function ####################################

replicate_condition_data <- function(condition_id, rep_id, study = 1, 
                                     results_file = NULL, base_path = ".") {
  
  # required packages (!)
  required_packages <- c("lavaan", "covsim", "rvinecopulib")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop("packages not installed: ", paste(missing_packages, collapse = ", "))
  }
  
  library(lavaan); library(covsim); library(rvinecopulib)
  
  # load appropriate models based on study
  if (study == 1) {
    models_path <- file.path(base_path, "Simulations/Study_1/Simulation/Models(1).RData")
    default_file <- file.path(base_path, "Simulations/Study_1/Data/Data_Study_1_final.RData")
  } else if (study == 2) {
    models_path <- file.path(base_path, "Simulations/Study_2/Simulation/Models(2).RData")
    default_file <- file.path(base_path, "Simulations/Study_2/Data/Data_Study_2_final.RData")
  } else {
    stop("simulation must be 1 or 2")
  }
  
  # check if models file exists
  if (!file.exists(models_path)) {
    stop("Models(X).RData file not found: ", models_path)
  }
  load(models_path) 
  
  # load GenerateData function
  generate_data_path <- file.path(base_path, "Simulations/GenerateData.R")
  if (!file.exists(generate_data_path)) {
    stop("GenerateData.R not found: ", generate_data_path)
  }
  source(generate_data_path)
  
  # load results file
  if (is.null(results_file)) {
    results_file <- default_file
  }
  
  if (!file.exists(results_file)) {
    stop("results file not found: ", results_file)
  }
  
  load(results_file) 
  
  # validate condition_id
  if (length(all_results) < condition_id) {
    stop("condition ", condition_id, " not found in results (max: ", 
         length(all_results), ")")
  }
  
  # extract condition data
  cond_data <- all_results[[condition_id]]
  condition <- cond_data$condition
  results <- cond_data$results
  
  # check if rep_id exists based on available methods
  if (study == 1) {
    # study 1: lms, qml, upi, lsam
    max_reps <- max(
      length(results$lms_tables),
      length(results$qml_tables),
      length(results$upi_tables),
      length(results$lsam_tables),
      na.rm = TRUE
    )
  } else {
    # study 2: lsam, qml, upi
    max_reps <- max(
      length(results$lsam_tables),
      length(results$qml_tables),
      length(results$upi_tables),
      na.rm = TRUE
    )
  }
  
  if (rep_id > max_reps) {
    stop("replication ", rep_id, " not found (max: ", max_reps, ")")
  }
  
  # RNG state - this is essential for replication
  if (is.null(results$rng_states)) {
    stop("RNG states not found in results. Cannot replicate exact data.")
  }
  
  rng_state <- results$rng_states[[rep_id]]
  if (is.null(rng_state)) {
    stop("RNG state for replication ", rep_id, " not found. Cannot replicate exact data.")
  }
  
  # restore RNG state
  assign(".Random.seed", rng_state, envir = .GlobalEnv)
  
  # get population model
  model_name <- condition$model_name
  population_model <- all_models[[model_name]]
  
  if (is.null(population_model)) {
    stop("Model ", model_name, " not found")
  }
  
  # distribution parameters for VITA-based generation
  distributions <- list(
    normal = list(
      distr.exo = "normal",
      nonnormal.shape = NULL,
      nonnormal.rate = NULL
    ),
    nonnormal = list(
      distr.exo = "nonnormal",
      nonnormal.shape = if (study == 1) c(1, 1) else c(1, 1, 1),  # 2 exog for Sim 1, 3 for Sim 2
      nonnormal.rate = if (study == 1) c(1, 1) else c(1, 1, 1)
    ),
    uniform = list(
      distr.exo = "uniform", 
      nonnormal.shape = NULL,
      nonnormal.rate = NULL
    )
  )
  
  # distribution name from condition
  # handle both old naming (Distribution) and potential new naming (generation_distribution)
  dist_name <- if (!is.null(condition$generation_distribution)) {
    tolower(condition$generation_distribution)
  } else if (!is.null(condition$Distribution)) {
    tolower(condition$Distribution)
  } else {
    stop("distribution information not found in condition")
  }
  
  dist_params <- distributions[[dist_name]]
  
  if (is.null(dist_params)) {
    stop("unknown distribution: ", dist_name)
  }
  
  # generate data using VITA-based approach
  tryCatch({
    Data <- GenerateData(
      model           = population_model,
      N               = condition$N,
      distr.exo       = dist_params$distr.exo,
      nonnormal.shape = dist_params$nonnormal.shape,
      nonnormal.rate  = dist_params$nonnormal.rate,
      distr.epsilon   = "normal",
      distr.zeta      = "normal",
      add.eta         = FALSE,
      return.info     = TRUE
    )
  }, error = function(e) {
    stop("failed to generate data for condition ", condition_id, 
         ", rep ", rep_id, ": ", e$message)
  })
  
  # add attributes with condition info
  attr(Data, "condition") <- condition
  attr(Data, "condition_id") <- condition_id
  attr(Data, "rep_id") <- rep_id
  attr(Data, "study") <- study
  attr(Data, "model_name") <- model_name
  attr(Data, "distribution") <- dist_name
  
  # return the data with all attributes preserved
  Data
}

# example:
data_replicate <- replicate_condition_data(
  condition_id = 1, 
  rep_id = 1,
  study = 1
)

# following this one now may proceed to load the analysis model and fit with whatever the estimation desired
# you should obtain the same results as reported in the specific .RData file

