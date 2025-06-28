#### 1. General Information ####

# note:
# This code is not processing multiple conditions in parallel. 
# It processes the conditions one after another. 
# The parallelization happens only within each condition, 
# where multiple replications are run simultaneously across X cores

# note:
# Additionally, starting parameters are estimated 
# using the double-centering approach, 
# and the means of the observed variables are used to 
# generate good starting parameters for faster convergence


# Other scripts needed
source("Simulation_Scripts/Simulation_Study_1/Simulation_1/GenerateData.R") # for generating data
source("Simulation_Scripts/Simulation_Study_1/Simulation_1/Methods(1).R") # for the estimation approaches
source("Simulation_Scripts/Simulation_Study_1/Simulation_1/Models(1).R") # for fit and pop models
source("Simulation_Scripts/Simulation_Study_1/Simulation_1/Design(1).R") # for the conditions

RNGkind("Mersenne-Twister", "Inversion", "Rejection")

dir.create("sim_results", showWarnings = FALSE)
timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
results_dir <- paste0("sim_results/run_", timestamp)
dir.create(results_dir, showWarnings = FALSE)

# ============================================================================
# PRE-COMPUTE POPULATION VARIANCES FOR ALL UNIQUE CONDITIONS
# ============================================================================

cat("\n=== Pre-computing population variances for all unique conditions ===\n")

# unique combinations that affect population variances
unique_pop_conditions <- unique(conditions[, c("Population", "Distribution", "Exo_method", "Rel")])
unique_pop_conditions <- unique_pop_conditions[order(
  unique_pop_conditions$Population, 
  unique_pop_conditions$Distribution,
  unique_pop_conditions$Exo_method,
  unique_pop_conditions$Rel
), ]

cat(sprintf("Found %d unique population conditions\n", nrow(unique_pop_conditions)))

pop_variances_list <- list()

# population variances for each unique combination
for(i in 1:nrow(unique_pop_conditions)) {
  cond <- unique_pop_conditions[i, ]
  
  # key for this combination
  key <- paste(cond$Population, cond$Distribution, cond$Exo_method, cond$Rel, sep = "_")
  
  cat(sprintf("\nComputing population variances for: %s\n", key))
  
  # distribution parameters
  skewness <- rep(ifelse(cond$Distribution == "normal", 0, 2), 2)
  excesskurtosis <- rep(ifelse(cond$Distribution == "normal", 0, 7), 2)
  
  # population variances
  pop_data <- GenerateData(
    model = get(cond$Population),
    N = 10,  # this does not matter now: we only need variances and for that N.pop matters, but:
    # 10 to avoid "Error in apply(EXO, 2, stats::var) : dim(X) must have a positive length"
    compute_pop_vars = TRUE,
    N.pop = 5000000,  # !!
    skewness = skewness,
    excesskurtosis = excesskurtosis,
    exo.mean = exo.mean, 
    distr.exo = cond$Exo_method,
    distr.zeta = "normal",
    distr.epsilon = "normal", 
    rel = cond$Rel,
    R2 = R2, 
    seed = 123
  )
  
  pop_variances_list[[key]] <- list(
    pop_var_nozeta = attr(pop_data, "pop_var_nozeta"),
    eta_pop_vars = attr(pop_data, "eta_pop_vars")
  )
  
  cat("  eta_pop_vars:\n")
  print(pop_variances_list[[key]]$eta_pop_vars)
}

save(pop_variances_list, file = paste0(results_dir, "/population_variances.RData"))
cat("\n=== Population variance computation complete ===\n")

# parallel backend
n_cores <- detectCores() - 2
cl <- makeCluster(n_cores)
registerDoParallel(cl)

clusterExport(cl, c("GenerateData", "method_analytic", "method_dblcent", "method_sam", 
                    "population.linear.model", "population.full.model", 
                    "fit.full.model", "pop_variances_list",
                    "exo.mean", "R2", "conditions"))

# ============================================================================
# SIMULATION PART
# ============================================================================

# ============================================================================
# FUNCTION FOR SIMULATION - LMS
# ============================================================================

# For now (10/04), process_replication only for LMS in parallel
process_lms_only <- function(i, cond, n_params, skewness, excesskurtosis, pop_var_nozeta, eta_pop_vars) {
  # container for LMS 
  lms_res <- list(
    estimates = rep(NA, n_params),
    se = rep(NA, n_params),
    pvals = rep(NA, n_params),
    ci_lower = rep(NA, n_params),  
    ci_upper = rep(NA, n_params),  
    timing = NA
  )
  
  Data <- try(GenerateData(
    model = get(conditions$Population[cond]),
    N = conditions$N[cond],
    skewness = skewness,
    excesskurtosis = excesskurtosis,
    exo.mean = exo.mean,
    distr.exo = conditions$Exo_method[cond],
    distr.zeta = "normal",
    distr.epsilon = conditions$Epsilon[cond],
    rel = conditions$Rel[cond],
    R2 = R2,
    add.eta = FALSE,
    compute_pop_vars = FALSE,  # don't compute, use provided
    pop_var_nozeta = pop_var_nozeta,  # use pre-computed
    eta_pop_vars = eta_pop_vars  # use pre-computed
  ), silent = TRUE)
  
  if(!inherits(Data, "try-error")) {
    analysis_model <- get(conditions$Analysis_model[cond])
    
    start_time_method <- Sys.time()
    result <- try(method_analytic(Data = Data, model.fit = analysis_model, method = "lms"))
    lms_res$timing <- as.numeric(difftime(Sys.time(), start_time_method, units = "secs"))
    
    if(!inherits(result, "try-error")) {
      lms_res$estimates <- result$Estimates
      lms_res$se <- result$`Standard Errors`
      lms_res$pvals <- result$`P-values`
      lms_res$ci_lower <- result$CI_lower  
      lms_res$ci_upper <- result$CI_upper  
    }
  }
  
  # So, here we return only LMS results as well as the data for each iteration
  # Then I use for the remaining methods accordingly
  list(data = if(!inherits(Data, "try-error")) Data else NULL,
       lms_res = lms_res)
}

# ============================================================================
# MAIN SIMULATION LOOP (MODIFIED)
# ============================================================================

start_time <- Sys.time()
all_results <- list()

# simulation for all conditions
for(cond in 1:nrow(conditions)) {
  condition_start <- Sys.time()
  cat("\nRunning condition", cond, "of", nrow(conditions), "\n")
  cat("Parameters: N =", conditions$N[cond], 
      ", rel =", conditions$Rel[cond], 
      ", pop_model =", conditions$Population[cond], 
      ", analy_model =", conditions$Analysis_model[cond],
      ", distribution =", conditions$Distribution[cond], 
      ", exo_method =", conditions$Exo_method[cond],
      ", epsilon =", conditions$Epsilon[cond], "\n")
  
  skewness <- rep(ifelse(conditions$Distribution[cond] == "normal", 0, 2), 2)
  excesskurtosis <- rep(ifelse(conditions$Distribution[cond] == "normal", 0, 7), 2)
  
  # RETRIEVE PRE-COMPUTED POPULATION VARIANCES
  key <- paste(conditions$Population[cond], 
               conditions$Distribution[cond], 
               conditions$Exo_method[cond], 
               conditions$Rel[cond], 
               sep = "_")
  
  pop_vars <- pop_variances_list[[key]]
  pop_var_nozeta <- pop_vars$pop_var_nozeta
  eta_pop_vars <- pop_vars$eta_pop_vars
  
  cat("Using pre-computed population variances for:", key, "\n")
  
  # N of parameters as a function of the analysis model
  n_params <- ifelse(conditions$Analysis_model[cond] == "fit.full.model", 5, 3)
  
  res <- list(
    lms = array(NA, dim = c(rep, n_params, 5), 
                dimnames = list(NULL, NULL, c("beta", "se", "pval", "ci_lower", "ci_upper"))),  
    qml = array(NA, dim = c(rep, n_params, 5), 
                dimnames = list(NULL, NULL, c("beta", "se", "pval", "ci_lower", "ci_upper"))),  
    dblcent = array(NA, dim = c(rep, n_params, 5),  
                dimnames = list(NULL, NULL, c("beta", "se", "pval", "ci_lower", "ci_upper"))),  
    sam = array(NA, dim = c(rep, n_params, 5),  
                dimnames = list(NULL, NULL, c("beta", "se", "pval", "ci_lower", "ci_upper"))),  
    timing = data.frame(
      lms = numeric(rep),
      qml = numeric(rep),
      dblcent = numeric(rep),
      sam = numeric(rep)
    )
  )
  
  # parallel processing ONLY for LMS 
  cat("\nRunning LMS in parallel\n")
  lms_parallel_results <- foreach(i = seq_len(rep), 
                                  .packages = c("modsem", "lavaan", "covsim", "faux"), 
                                  .errorhandling = "pass",
                                  .options.RNG = 123 + cond * 1000) %dorng% {
                                    process_lms_only(i, cond, n_params, skewness, excesskurtosis, 
                                                     pop_var_nozeta, eta_pop_vars)
                                  }
  
  # RNG states for this condition
  rng_states <- attr(lms_parallel_results, "rng")
  
  # other methods sequentially
  cat("\nRunning QML, UCA, and SAM sequentially\n")
  for(i in 1:rep) {
    # skip if data generation failed
    if(is.null(lms_parallel_results[[i]]$data)) {
      next
    }
    
    # LMS results
    # order is preserved - good design from foreach
    res$lms[i, , 1] <- lms_parallel_results[[i]]$lms_res$estimates
    res$lms[i, , 2] <- lms_parallel_results[[i]]$lms_res$se
    res$lms[i, , 3] <- lms_parallel_results[[i]]$lms_res$pvals
    res$lms[i, , 4] <- lms_parallel_results[[i]]$lms_res$ci_lower  
    res$lms[i, , 5] <- lms_parallel_results[[i]]$lms_res$ci_upper  
    res$timing$lms[i] <- lms_parallel_results[[i]]$lms_res$timing
    
    # dataset from parallel results (!) - make sure we get the same exact data used in lms
    Data <- lms_parallel_results[[i]]$data
    analysis_model <- get(conditions$Analysis_model[cond])
    
    # QML method
    start_time_method <- Sys.time()
    result <- suppressWarnings(try(method_analytic(Data = Data, model.fit = analysis_model, method = "qml")))
    res$timing$qml[i] <- as.numeric(difftime(Sys.time(), start_time_method, units = "secs"))
    
    if(!inherits(result, "try-error")) {
      res$qml[i, , 1] <- result$Estimates
      res$qml[i, , 2] <- result$`Standard Errors`
      res$qml[i, , 3] <- result$`P-values`
      res$qml[i, , 4] <- result$CI_lower  
      res$qml[i, , 5] <- result$CI_upper  
    }
    
    # UCA method
    start_time_method <- Sys.time()
    result <- suppressWarnings(try(method_dblcent(Data = Data, model.fit = analysis_model)))
    res$timing$dblcent[i] <- as.numeric(difftime(Sys.time(), start_time_method, units = "secs"))
    
    if(!inherits(result, "try-error")) {
      res$dblcent[i, , 1] <- result$Estimates
      res$dblcent[i, , 2] <- result$`Standard Errors`
      res$dblcent[i, , 3] <- result$`P-values`
      res$dblcent[i, , 4] <- result$CI_lower  
      res$dblcent[i, , 5] <- result$CI_upper 
    }
    
    # SAM method
    start_time_method <- Sys.time()
    result <- suppressWarnings(try(method_sam(Data = Data, model.fit = analysis_model)))
    res$timing$sam[i] <- as.numeric(difftime(Sys.time(), start_time_method, units = "secs"))
    
    if(!inherits(result, "try-error")) {
      res$sam[i, , 1] <- result$Estimates
      res$sam[i, , 2] <- result$`Standard Errors`
      res$sam[i, , 3] <- result$`P-values`
      res$sam[i, , 4] <- result$CI_lower  
      res$sam[i, , 5] <- result$CI_upper  
    }
  }
  
  all_results[[cond]] <- list(
    condition = conditions[cond, ],
    results = res,
    rng_states = rng_states,
    population_variances = pop_vars  # population variances used
  )
  
  # results for (this) condition
  condition_filename <- sprintf(
    "%s/condition_%d_N%d_Rel%s_%s_%s_%s_%s.RData",
    results_dir, 
    cond,
    conditions$N[cond],
    conditions$Rel[cond],
    substr(conditions$Population[cond], 11, 20),
    conditions$Distribution[cond],
    conditions$Exo_method[cond],
    conditions$Epsilon[cond]
  )
  
  tryCatch({
    save(res, file = condition_filename)
  }, error = function(e) {
    warning(paste("Failed to save condition results:", e$message))
  })
  
  # all results until 10th condition
  if(cond %% 10 == 0 || cond == nrow(conditions)) {
    tryCatch({
      save(all_results, rng_states, conditions, 
           file = sprintf("%s/all_results_upto_condition_%d.RData", results_dir, cond))
    }, error = function(e) {
      warning(paste("Failed to save checkpoint:", e$message))
    })
  }
  
  # removing large objects from memory
  rm(lms_parallel_results)  
  gc()    
  
  # progress info
  condition_time <- difftime(Sys.time(), condition_start, units = "mins")
  cat(sprintf("\nCondition %d completed in %.2f minutes\n", cond, condition_time))
  cat(sprintf("Progress: %.1f%% complete (%d of %d conditions)\n", 
              100 * cond / nrow(conditions), cond, nrow(conditions)))
}

stopCluster(cl)

total_time <- difftime(Sys.time(), start_time, units = "hours")
cat(sprintf("\nTotal simulation completed in %.2f hours\n", total_time))
save(all_results, file = paste0(results_dir, "/final_results.RData"))