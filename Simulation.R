#### 1. General Information ####

# This code is not processing multiple conditions in parallel. 
# It processes the conditions one after another. 
# The parallelization happens only within each condition, 
# where multiple replications are run simultaneously across X cores.

# Other scripts needed
source("GenerateData.R")
source("Methods.R")
source("Models.R")
source("Design.R")

# parallel backend
n_cores <- detectCores() - 2
cl <- makeCluster(n_cores)
registerDoParallel(cl)

clusterExport(cl, c("GenerateData", "method_analytic", "method_uca", "method_sam", "population.interaction.model",
                    "population.linear.model", "population.full.model", "fit.interaction.model", "fit.full.model"))

dir.create("sim_results", showWarnings = FALSE)
timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
results_dir <- paste0("sim_results/run_", timestamp)
dir.create(results_dir, showWarnings = FALSE)

start_time <- Sys.time()
all_results <- list()

# For now (10/04), process_replication only for LMS in parallel
process_lms_only <- function(i, cond, n_params, skewness, excesskurtosis) {
  # Container just for LMS results 
  lms_res <- list(
    estimates = rep(NA, n_params),
    se = rep(NA, n_params),
    pvals = rep(NA, n_params),
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
    #target.var = target.var,
    R2 = R2,
    add.eta = FALSE), silent = TRUE)
  
  if(!inherits(Data, "try-error")) {
    analysis_model <- get(conditions$Analysis_model[cond])
    
    start_time_method <- Sys.time()
    result <- try(method_analytic(Data = Data, model.fit = analysis_model, method = "lms"))
    lms_res$timing <- as.numeric(difftime(Sys.time(), start_time_method, units = "secs"))
    
    if(!inherits(result, "try-error")) {
      lms_res$estimates <- result$Estimates
      lms_res$se <- result$`Standard Errors`
      lms_res$pvals <- result$`P-values`
    }
  }
  
  # So, here we return only LMS results as well as the data for each iteration
  # Then I use for the remaining methods accordingly
  list(data = if(!inherits(Data, "try-error")) Data else NULL,
       lms_res = lms_res)
}


# Simulation for all conditions
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
  
  # N of parameters as a function of the analysis model
  n_params <- ifelse(conditions$Analysis_model[cond] == "fit.full.model", 5, 3)
  
  res <- list(
    lms = array(NA, dim = c(rep, n_params, 3),
                dimnames = list(NULL, NULL, c("beta", "se", "pval"))),
    qml = array(NA, dim = c(rep, n_params, 3),
                dimnames = list(NULL, NULL, c("beta", "se", "pval"))),
    uca = array(NA, dim = c(rep, n_params, 3),
                dimnames = list(NULL, NULL, c("beta", "se", "pval"))),
    sam = matrix(NA, rep, n_params),
    timing = data.frame(
      lms = numeric(rep),
      qml = numeric(rep),
      uca = numeric(rep),
      sam = numeric(rep)
    )
  )
  
  # Parallel processing ONLY for LMS 
  cat("\nRunning LMS in parallel\n")
  lms_parallel_results <- foreach(i = seq_len(rep), 
                                  .packages = c("modsem", "lavaan", "covsim"), 
                                  .errorhandling = "pass",
                                  .options.RNG = conditions$Seed[cond]) %dorng% {
                                    process_lms_only(i, cond, n_params, skewness, excesskurtosis)
                                  }
  
  # RNG states for this condition
  rng_states <- attr(lms_parallel_results, "rng")
  
  # Process other methods sequentially
  cat("\nRunning QML, UCA, and SAM sequentially\n")
  for(i in 1:rep) {
    # skip if data generation failed
    if(is.null(lms_parallel_results[[i]]$data)) {
      next
    }
    
    # LMS results
    res$lms[i, , 1] <- lms_parallel_results[[i]]$lms_res$estimates
    res$lms[i, , 2] <- lms_parallel_results[[i]]$lms_res$se
    res$lms[i, , 3] <- lms_parallel_results[[i]]$lms_res$pvals
    res$timing$lms[i] <- lms_parallel_results[[i]]$lms_res$timing
    
    # Dataset from parallel results
    Data <- lms_parallel_results[[i]]$data
    analysis_model <- get(conditions$Analysis_model[cond])
    
    # QML method
    start_time_method <- Sys.time()
    result <- try(method_analytic(Data = Data, model.fit = analysis_model, method = "qml"))
    res$timing$qml[i] <- as.numeric(difftime(Sys.time(), start_time_method, units = "secs"))
    
    if(!inherits(result, "try-error")) {
      res$qml[i, , 1] <- result$Estimates
      res$qml[i, , 2] <- result$`Standard Errors`
      res$qml[i, , 3] <- result$`P-values`
    }
    
    # UCA method
    start_time_method <- Sys.time()
    result <- try(method_uca(Data = Data, model.fit = analysis_model))
    res$timing$uca[i] <- as.numeric(difftime(Sys.time(), start_time_method, units = "secs"))
    
    if(!inherits(result, "try-error")) {
      res$uca[i, , 1] <- result$Estimates
      res$uca[i, , 2] <- result$`Standard Errors`
      res$uca[i, , 3] <- result$`P-values`
    }
    
    # SAM method
    start_time_method <- Sys.time()
    result <- try(method_sam(Data = Data, model.fit = analysis_model))
    res$timing$sam[i] <- as.numeric(difftime(Sys.time(), start_time_method, units = "secs"))
    
    if(!inherits(result, "try-error")) {
      # As of now, still no standard errors; no SE and p-value
      res$sam[i, ] <- result$Estimates
    }
  }
  
  all_results[[cond]] <- list(
    condition = conditions[cond, ],
    results = res,
    rng_states = rng_states  # rng states
  )
  
  # RResults for this condition
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
  
  # All results periodically (probably stop after the 50th)
  if(cond %% 10 == 0 || cond == nrow(conditions)) {
    tryCatch({
      save(all_results, conditions, 
           file = sprintf("%s/all_results_upto_condition_%d.RData", results_dir, cond))
    }, error = function(e) {
      warning(paste("Failed to save checkpoint:", e$message))
    })
  }
  
  # removing large objects from memory
  rm(lms_parallel_results)  
  gc()    
  
  # Progress information
  condition_time <- difftime(Sys.time(), condition_start, units = "mins")
  cat(sprintf("\nCondition %d completed in %.2f minutes\n", cond, condition_time))
  cat(sprintf("Progress: %.1f%% complete (%d of %d conditions)\n", 
              100 * cond / nrow(conditions), cond, nrow(conditions)))
}

stopCluster(cl)

total_time <- difftime(Sys.time(), start_time, units = "hours")
cat(sprintf("\nTotal simulation completed in %.2f hours\n", total_time))
save(all_results, file = paste0(results_dir, "/final_results.RData"))