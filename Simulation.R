#library(parallel)
library(foreach)
library(doParallel)

# Other scripts needed
source("GenerateData.R")
source("Methods.R")
source("Models.R")
source("Design.R")

# Setup parallel backend
n_cores <- detectCores() - 2 # 8 in total, we use 6; change this when using the semlab pc
cl <- makeCluster(n_cores)
registerDoParallel(cl)

clusterExport(cl, c("GenerateData", "method_analytic", "method_uca", "method_sam", "population.interaction.model",
                    "population.linear.model", "population.full.model", "fit.interaction.model", "fit.full.model"))

# Directories for results
dir.create("sim_results", showWarnings = FALSE)
timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
results_dir <- paste0("sim_results/run_", timestamp)
dir.create(results_dir, showWarnings = FALSE)

start_time <- Sys.time()
all_results <- list()

# Single replication. This is used below within foreach
process_replication <- function(i, cond, n_params, skewness, excesskurtosis) {
  set.seed(1234 + i + (cond * 1000)) # Important seed like this (not sure how to use streams yet)
  
  # Local result structure
  # Later we aggregate all into one (combining the different processing levels)
  local_res <- list(
    lms = array(NA, dim = c(1, n_params, 3)),
    qml = array(NA, dim = c(1, n_params, 3)),
    uca = array(NA, dim = c(1, n_params, 3)),
    sam = rep(NA, n_params),
    timing = data.frame(
      lms = NA,
      qml = NA,
      uca = NA,
      sam = NA
    )
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
    target.var = target.var,
    R2 = R2,
    add.eta = FALSE), silent = TRUE)
  
  if(!inherits(Data, "try-error")) {
    analysis_model <- get(conditions$Analysis_model[cond])
    
    # LMS and QML methods
    for(m in c("lms", "qml")) {
      start_time_method <- Sys.time()
      result <- try(method_analytic(Data = Data, model.fit = analysis_model, method = m))
      local_res$timing[[m]] <- as.numeric(difftime(Sys.time(), start_time_method, units = "secs"))
      
      if(!inherits(result, "try-error")) {
        local_res[[m]][1, , 1] <- result$Estimates
        local_res[[m]][1, , 2] <- result$`Standard Errors`
        local_res[[m]][1, , 3] <- result$`P-values`
      }
    }
    
    # UCA method
    start_time_method <- Sys.time()
    result <- try(method_uca(Data = Data, model.fit = analysis_model))
    local_res$timing$uca <- as.numeric(difftime(Sys.time(), start_time_method, units = "secs"))
    
    if(!inherits(result, "try-error")) {
      local_res$uca[1, , 1] <- result$Estimates
      local_res$uca[1, , 2] <- result$`Standard Errors`
      local_res$uca[1, , 3] <- result$`P-values`
    }
    
    # SAM method
    start_time_method <- Sys.time()
    result <- try(method_sam(Data = Data, model.fit = analysis_model))
    local_res$timing$sam <- as.numeric(difftime(Sys.time(), start_time_method, units = "secs"))
    
    if(!inherits(result, "try-error")) {
      # As of now, still no standard errors; hence, no SE and p-value
      local_res$sam <- result$Estimates
    }
  }
  local_res
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
  
  # Distribution parameter assignment
  skewness <- rep(ifelse(conditions$Distribution[cond] == "normal", 0, 2), 2)
  excesskurtosis <- rep(ifelse(conditions$Distribution[cond] == "normal", 0, 7), 2)
  
  # Number of parameters as a function of the analysis model
  n_params <- ifelse(conditions$Analysis_model[cond] == "fit.full.model", 5, 3)
  
  # Results structure
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
  
  # Parallel processing of replications
  results <- foreach(i = seq_len(rep), 
                     .packages = c("modsem", "lavaan", "covsim"), 
                     .errorhandling = "pass") %dopar% {
                       process_replication(i, cond, n_params, skewness, excesskurtosis)
                     }
  
  # Get valid results indices
  valid_indices <- which(!sapply(results, inherits, "try-error"))
  
  # Process results for valid indices
  for(i in valid_indices) {
    # Array methods (lms, qml, uca)
    for(method in c("lms", "qml", "uca")) {
      res[[method]][i, , ] <- results[[i]][[method]][1, , ]
    }
    
    # SAM results
    res$sam[i, ] <- results[[i]]$sam
    
    # Timing data
    res$timing[i, ] <- results[[i]]$timing
  }
  
  all_results[[cond]] <- list(
    condition = conditions[cond, ],
    results = res
  )
  
  # Save results for this condition
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
  
  save(res, file = condition_filename)
  
  # Save all results periodically
  if(cond %% 10 == 0 || cond == nrow(conditions)) {
    save(all_results, conditions, 
         file = sprintf("%s/all_results_upto_condition_%d.RData", results_dir, cond))
  }
  
  # Update progress information
  condition_time <- difftime(Sys.time(), condition_start, units = "mins")
  cat(sprintf("\nCondition %d completed in %.2f minutes\n", cond, condition_time))
  cat(sprintf("Progress: %.1f%% complete (%d of %d conditions)\n", 
              100 * cond / nrow(conditions), cond, nrow(conditions)))
}

stopCluster(cl)

total_time <- difftime(Sys.time(), start_time, units = "hours")
cat(sprintf("\nTotal simulation completed in %.2f hours\n", total_time))
# Save final results
save(all_results, conditions, file = paste0(results_dir, "/final_results.RData"))
