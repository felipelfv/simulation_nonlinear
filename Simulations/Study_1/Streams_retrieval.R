replicate_condition_data <- function(condition_id, rep_id, study = 1, 
                                     results_file = NULL) {
  
  library(lavaan); library(covsim); library(copula)
  
  # appropriate models based on study
  if (study == 1) {
    load("Simulations/Study_1/Simulation/Models(1).RData")  # all_models
    default_file <- "Simulations/Study_1/Data/Results_Study_1_final.RData"
  } else if (study == 2) {
    load("Simulations/Study_2/Simulation/Models(2).RData")  # all_models
    default_file <- "Simulations/Study_2/Data/Results_Study_2_final.RData"
  } else {
    stop("Study must be 1 or 2")
  }
  
  source("Simulations/GenerateData.R")  
  
  # load results file
  if (is.null(results_file)) {
    results_file <- default_file
  }
  
  if (!file.exists(results_file)) {
    stop("Results file not found: ", results_file)
  }
  
  load(results_file)
  
  if (length(all_results) < condition_id) {
    stop("Condition ", condition_id, " not found in results (max: ", 
         length(all_results), ")")
  }
  
  # condition data
  cond_data <- all_results[[condition_id]]
  condition <- cond_data$condition
  results <- cond_data$results
  
  # if rep_id exists
  max_reps <- if (study == 1) {
    max(length(results$lms_tables), 
        length(results$qml_tables),
        length(results$dblcent_tables),
        length(results$sam_tables))
  } else {
    max(length(results$qml_tables),
        length(results$dblcent_tables),
        length(results$sam_tables))
  }
  
  if (rep_id > max_reps) {
    stop("Replication ", rep_id, " not found (max: ", max_reps, ")")
  }
  
  # RNG state as this is essential for replication
  if (is.null(results$rng_states)) {
    stop("RNG states not found in results. Cannot replicate exact data.")
  }
  
  rng_state <- results$rng_states[[rep_id]]
  if (is.null(rng_state)) {
    stop("RNG state for replication ", rep_id, " not found. Cannot replicate exact data.")
  }
  
  assign(".Random.seed", rng_state, envir = .GlobalEnv)
  
  # population model
  model_name <- condition$model_name
  population_model <- all_models[[model_name]]
  
  if (is.null(population_model)) {
    stop("Model ", model_name, " not found")
  }
  
  # distribution parameters
  distributions <- list(
    normal    = list(skewness = c(0, 0, 0), 
                     excesskurtosis = c(0, 0, 0), 
                     distr.exo = "normal.rIG"),
    nonnormal = list(skewness = c(2, 2, 2), 
                     excesskurtosis = c(7, 7, 7), 
                     distr.exo = "nonnormal.rIG"),
    uniform   = list(skewness = c(0, 0, 0), 
                     excesskurtosis = c(0, 0, 0), 
                     distr.exo = "unif")
  )
  
  # study 1, use 2-element vectors; for Study 2, use 3-element vectors
  if (study == 1) {
    distributions$normal$skewness <- c(0, 0)
    distributions$normal$excesskurtosis <- c(0, 0)
    distributions$nonnormal$skewness <- c(2, 2)
    distributions$nonnormal$excesskurtosis <- c(7, 7)
    distributions$uniform$skewness <- c(0, 0)
    distributions$uniform$excesskurtosis <- c(0, 0)
  }
  
  dist_params <- distributions[[tolower(condition$Distribution)]]
  
  Data <- GenerateData(
    model          = population_model,
    N              = condition$N,
    skewness       = dist_params$skewness,
    excesskurtosis = dist_params$excesskurtosis,
    distr.exo      = dist_params$distr.exo,
    distr.epsilon  = "normal",
    distr.zeta     = "normal",
    add.eta        = FALSE,
    return.info    = TRUE
  )
  
  # data with condition info as attributes
  attr(Data, "condition") <- condition
  attr(Data, "condition_id") <- condition_id
  attr(Data, "rep_id") <- rep_id
  attr(Data, "study") <- study
  
  Data
}

# example:
data_to_replicate <- replicate_condition_data(
  condition_id = 1, 
  rep_id = 15,
  study = 1
)
