# scripts needed
source("Simulation_Scripts/GenerateData(2).R")
source("Simulation_Scripts/Simulation_Study_2/Method(2).R")
source("Simulation_Scripts/Simulation_Study_2/Models(2).R")
source("Simulation_Scripts/Simulation_Study_2/Design(2).R")

rep <- 100

dir.create("sim_results", showWarnings = FALSE)
timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
results_dir <- paste0("sim_results_2/run_", timestamp)
dir.create(results_dir, showWarnings = FALSE)

start_time <- Sys.time()
all_results <- list()

param_names <- c("eta4~eta1", "eta4~eta2", "eta4~eta3", "eta4~norm", 
                 "eta4~eta1:eta3", "eta4~eta2:eta3", "eta4~eta3:eta3", 
                 "eta4~eta5:eta2", "eta4~eta5:eta5", "eta5~eta1:eta2", 
                 "eta3~gender", "eta3~age", "eta3~norm", 
                 "norm~gender", "norm~age")
n_params <- length(param_names)

# for all conditions
for(cond in 1:nrow(conditions)) {
  condition_start <- Sys.time()
  cat("\nRunning condition", cond, "of", nrow(conditions), "\n")
  cat("Parameters: N =", conditions$N[cond], 
      ", rel =", conditions$Rel[cond], 
      ", pop_model =", conditions$Population[cond], 
      ", analy_model =", conditions$Analysis_model[cond],
      ", distribution =", conditions$Distribution[cond], 
      ", epsilon =", conditions$Epsilon[cond], "\n")
  
  skewness <- rep(ifelse(conditions$Distribution[cond] == "normal", 0, 2), 4)
  excesskurtosis <- rep(ifelse(conditions$Distribution[cond] == "normal", 0, 7), 4)
  
  res <- list(
    sam_estimates = matrix(NA, nrow = rep, ncol = n_params,
                           dimnames = list(NULL, param_names)),
    timing = numeric(rep)
  )
  
  analysis_model <- get(conditions$Analysis_model[cond])
  
  # sequential processing for replications
  cat("\nRunning SAM sequentially\n")
  for(i in 1:rep) {
    
    Data <- try(GenerateData(
      model = get(conditions$Population[cond]),
      N = conditions$N[cond],
      skewness = skewness,
      excesskurtosis = excesskurtosis,
      exo.mean = exo.mean,
      distr.zeta = "normal",
      distr.epsilon = conditions$Epsilon[cond],
      rel = conditions$Rel[cond],
      R2 = R2,
      add.eta = FALSE), silent = TRUE)
    
    if(!inherits(Data, "try-error")) {
      # SAM method
      start_time_method <- Sys.time()
      result <- suppressWarnings(try(method_sam(Data = Data, model.fit = analysis_model)))
      res$timing[i] <- as.numeric(difftime(Sys.time(), start_time_method, units = "secs"))
      
      if(!inherits(result, "try-error")) {
        res$sam_estimates[i, ] <- result$Estimates
      }
    }
    
    # progress info within condition
    if(i %% 10 == 0 || i == rep) {
      cat(sprintf("\rProgress: %d/%d replications", i, rep))
    }
  }
  cat("\n")
  
  all_results[[cond]] <- list(
    condition = conditions[cond, ],
    results = res
  )
  
  # result for single cond
  condition_filename <- sprintf(
    "%s/condition_%d_N%d_Rel%s_%s_%s_%s.RData",
    results_dir, 
    cond,
    conditions$N[cond],
    conditions$Rel[cond],
    substr(conditions$Population[cond], 11, 20),
    conditions$Distribution[cond],
    conditions$Epsilon[cond]
  )
  
  tryCatch({
    save(res, file = condition_filename)
  }, error = function(e) {
    warning(paste("Failed to save condition results:", e$message))
  })
  
  if(cond %% 10 == 0 || cond == min(10, nrow(conditions))) {
    tryCatch({
      save(all_results, conditions, 
           file = sprintf("%s/all_results_upto_condition_%d.RData", results_dir, cond))
    }, error = function(e) {
      warning(paste("Failed to save checkpoint:", e$message))
    })
  }
  
  condition_time <- difftime(Sys.time(), condition_start, units = "mins")
  cat(sprintf("\nCondition %d completed in %.2f minutes\n", cond, condition_time))
  cat(sprintf("Progress: %.1f%% complete (%d of %d conditions)\n", 
              100 * cond / nrow(conditions), cond, nrow(conditions)))
}

total_time <- difftime(Sys.time(), start_time, units = "hours")
cat(sprintf("\nTotal simulation completed in %.2f hours\n", total_time))
save(all_results, file = paste0(results_dir, "/final_results.RData"))