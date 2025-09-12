replicate_condition_data <- function(condition_id, rep_id, results_dir, 
                                     exo.mean = rep(0, 2), R2 = list("eta3" = 0.30)) {
  
  library(covsim)
  
  # necessary scripts
  source("Simulation_Scripts/Simulation_Study_1/Simulation_1/GenerateData.R")
  source("Simulation_Scripts/Simulation_Study_1/Simulation_1/Models(1).R")
  
  # population variances
  pop_var_file <- file.path(results_dir, "population_variances.RData")
  if (!file.exists(pop_var_file)) {
    stop("Population variances file not found in ", results_dir)
  }
  load(pop_var_file)  # Loads 'pop_variances_list'
  
  # full results to get condition info and RNG states
  checkpoint_files <- list.files(results_dir, pattern = "all_results_upto_condition_.*\\.RData", 
                                 full.names = TRUE)
  if (length(checkpoint_files) == 0) {
    final_results_file <- file.path(results_dir, "final_results.RData")
    if (file.exists(final_results_file)) {
      load(final_results_file)
    } else {
      stop("No checkpoint or final results files found in ", results_dir)
    }
  } else {
    for (file in rev(sort(checkpoint_files))) {
      load(file)  # Loads 'all_results'
      if (length(all_results) >= condition_id) break
    }
  }
  
  if (length(all_results) < condition_id) {
    stop("Condition ", condition_id, " not found in results files")
  }
  
  # RNG state and condition parameters
  cond_data <- all_results[[condition_id]]
  rng_state <- cond_data$rng_states[[rep_id]]
  cond_row <- cond_data$condition
  
  # RNG state
  assign(".Random.seed", rng_state, envir = .GlobalEnv)
  
  # distribution parameters
  is_normal <- cond_row$Distribution == "normal"
  skewness <- rep(ifelse(is_normal, 0, 2), 2)
  excesskurtosis <- rep(ifelse(is_normal, 0, 7), 2)
  
  # pre-computed population variances
  key <- paste(cond_row$Population, 
               cond_row$Distribution, 
               cond_row$Exo_method, 
               cond_row$Rel, 
               sep = "_")
  
  pop_vars <- pop_variances_list[[key]]
  if (is.null(pop_vars)) {
    stop("Population variances not found for key: ", key)
  }
  
  pop_var_nozeta <- pop_vars$pop_var_nozeta
  eta_pop_vars <- pop_vars$eta_pop_vars
  
  # original data with pre-computed population variances
  Data <- GenerateData(
    model = get(as.character(cond_row$Population)),
    N = cond_row$N,
    skewness = skewness,
    excesskurtosis = excesskurtosis,
    exo.mean = exo.mean,
    distr.exo = as.character(cond_row$Exo_method),
    distr.zeta = "normal",
    distr.epsilon = as.character(cond_row$Epsilon),
    rel = cond_row$Rel,
    R2 = R2,
    add.eta = FALSE,
    compute_pop_vars = FALSE,  # use provided
    pop_var_nozeta = pop_var_nozeta,  # pre-computed
    eta_pop_vars = eta_pop_vars  # pre-computed
  )
  
  Data
}

Data <- replicate_condition_data(
  condition_id = 1, 
  rep_id = 15, 
  results_dir = "sim_results/run_20250626_1435" 
)
