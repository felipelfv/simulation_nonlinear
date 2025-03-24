# Source required functions
source("GenerateData.R")
source("Methods.R")
source("Models.R")
source("Design.R")

# Setup parallel backend
n_cores <- detectCores() - 2  # 8 in total; 6 used
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Necessary functions that need to be exported to the cluster
clusterExport(cl, c("GenerateData", "method_analytic", "method_uca", "method_sam", 
                    "population.interaction.model", "population.linear.model", 
                    "population.full.model", "fit.interaction.model", "fit.full.model"))

# Directories for results
dir.create("sim_results", showWarnings = FALSE)
timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
results_dir <- paste0("sim_results/run_", timestamp)
dir.create(results_dir, showWarnings = FALSE)

start_time <- Sys.time()
all_results <- list()

# Debug mode for diagnostics
debug_mode <- FALSE  # Set to TRUE to see detailed debug output

# Process a single replication - can be used within foreach
process_replication <- function(i, cond, n_params, skewness, excesskurtosis) {
  set.seed(1234 + i + (cond * 1000)) # Important seed
  
  # Create local result structure with single replication
  local_res <- list(
    lms = array(NA, dim = c(1, n_params, 3)),
    qml = array(NA, dim = c(1, n_params, 3)),
    uca = array(NA, dim = c(1, n_params, 3)),
    sam = rep(NA, n_params),
    timing = data.frame(
      lms = NA_real_,
      qml = NA_real_,
      uca = NA_real_,
      sam = NA_real_
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
      local_res$sam <- result$Estimates
    }
  }
  
  return(local_res)
}

# Simulation for all conditions
for(cond in 1:nrow(conditions)) {
  condition_start <- Sys.time()
  cat("\nRunning condition", cond, "of", nrow(conditions), "\n")
  cat("Parameters: N =", conditions$N[cond], 
      ", rel =", conditions$Rel[cond], 
      ", pop_model =", conditions$Population[cond], 
      ", distribution =", conditions$Distribution[cond], 
      ", exo_method =", conditions$Exo_method[cond],
      ", epsilon =", conditions$Epsilon[cond], "\n")
  
  # Vectorized distribution parameter assignment
  skewness <- rep(ifelse(conditions$Distribution[cond] == "normal", 0, 2), 2)
  excesskurtosis <- rep(ifelse(conditions$Distribution[cond] == "normal", 0, 7), 2)
  
  # Determine number of parameters based on analysis model
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
  
  # Debug message
  if(debug_mode) {
    cat("Starting parallel processing with", rep, "replications\n")
  }
  
  # Parallel processing of replications - using the function now
  results <- foreach(i = seq_len(rep), 
                     .packages = c("modsem", "lavaan", "covsim"), 
                     .errorhandling = "pass") %dopar% {
                       process_replication(i, cond, n_params, skewness, excesskurtosis)
                     }
  
  # Debug information about results
  if(debug_mode) {
    cat("Parallel processing completed. Checking results:\n")
    cat("Number of results:", length(results), "\n")
    
    # Check if all results are errors
    error_count <- sum(sapply(results, function(x) inherits(x, "try-error")))
    cat("Number of error results:", error_count, "\n")
    
    # Sample a non-error result to check structure
    valid_idx <- which(!sapply(results, inherits, "try-error"))
    if(length(valid_idx) > 0) {
      sample_idx <- valid_idx[1]
      cat("Structure of sample result (index", sample_idx, "):\n")
      print(names(results[[sample_idx]]))
      for(method in c("lms", "qml", "uca")) {
        cat("Dimensions of", method, "array:", 
            paste(dim(results[[sample_idx]][[method]]), collapse="×"), "\n")
        cat("Sample values in", method, "array:", 
            paste(results[[sample_idx]][[method]][1,1,], collapse=", "), "\n")
      }
      cat("SAM result has", length(results[[sample_idx]]$sam), "elements\n")
      cat("Sample values in SAM:", paste(head(results[[sample_idx]]$sam), collapse=", "), "\n")
    } else {
      cat("No valid results to examine structure\n")
      # If all are errors, print the first error
      if(error_count > 0) {
        cat("First error message:", as.character(results[[1]]), "\n")
      }
    }
  }
  
  # Identify valid results using vectorization
  valid_results <- !sapply(results, inherits, "try-error")
  valid_indices <- which(valid_results)
  
  # Combine parallel results - use for loop for safety, but with vectorized valid indices
  for(i in valid_indices) {
    # Check if results[[i]] contains expected elements
    expected_items <- c("lms", "qml", "uca", "sam", "timing")
    missing_items <- expected_items[!expected_items %in% names(results[[i]])]
    
    if(length(missing_items) > 0) {
      if(debug_mode) {
        cat("Warning: Result", i, "is missing items:", 
            paste(missing_items, collapse=", "), "\n")
      }
      next  # Skip this iteration
    }
    
    # Process array methods (lms, qml, uca) with error checking
    for(method in c("lms", "qml", "uca")) {
      tryCatch({
        # If there's valid data, copy it
        if(!is.null(results[[i]][[method]]) && 
           is.array(results[[i]][[method]]) && 
           all(dim(results[[i]][[method]]) == c(1, n_params, 3))) {
          res[[method]][i, , ] <- results[[i]][[method]][1, , ]
          if(debug_mode && i == valid_indices[1]) {
            cat("Successfully copied", method, "data for first valid rep\n")
          }
        }
      }, error = function(e) {
        if(debug_mode) {
          cat("Error copying", method, "data for rep", i, ":", conditionMessage(e), "\n")
          cat("Dimensions - Source:", paste(dim(results[[i]][[method]]), collapse="×"), 
              "Target:", paste(dim(res[[method]][i,,]), collapse="×"), "\n")
        }
      })
    }
    
    # Process SAM results with error checking
    tryCatch({
      if(!is.null(results[[i]]$sam) && length(results[[i]]$sam) == n_params) {
        res$sam[i, ] <- results[[i]]$sam
        if(debug_mode && i == valid_indices[1]) {
          cat("Successfully copied SAM data for first valid rep\n")
        }
      }
    }, error = function(e) {
      if(debug_mode) {
        cat("Error copying SAM data for rep", i, ":", conditionMessage(e), "\n")
      }
    })
    
    # Process timing data with error checking
    tryCatch({
      if(!is.null(results[[i]]$timing) && 
         all(c("lms", "qml", "uca", "sam") %in% names(results[[i]]$timing))) {
        res$timing[i, ] <- results[[i]]$timing
        if(debug_mode && i == valid_indices[1]) {
          cat("Successfully copied timing data for first valid rep\n")
        }
      }
    }, error = function(e) {
      if(debug_mode) {
        cat("Error copying timing data for rep", i, ":", conditionMessage(e), "\n")
      }
    })
  }
  
  # Verify that res contains at least some data
  if(debug_mode) {
    # Check each method to see if there's any non-NA data
    for(method in c("lms", "qml", "uca")) {
      na_count <- sum(is.na(res[[method]]))
      total_elements <- prod(dim(res[[method]]))
      cat(method, "array: ", total_elements - na_count, "/", total_elements, 
          "elements filled (", 
          round(100 * (total_elements - na_count) / total_elements, 1), "%)\n", sep="")
    }
    
    # Check SAM matrix
    na_count <- sum(is.na(res$sam))
    total_elements <- prod(dim(res$sam))
    cat("SAM matrix: ", total_elements - na_count, "/", total_elements, 
        "elements filled (", 
        round(100 * (total_elements - na_count) / total_elements, 1), "%)\n", sep="")
  }
  
  all_results[[cond]] <- list(
    condition = conditions[cond, ],
    results = res
  )
  
  # Generate filename using sprintf
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
  
  # Also, periodically save all results
  if(cond %% 10 == 0 || cond == nrow(conditions)) {
    save(all_results, conditions, 
         file = sprintf("%s/all_results_upto_condition_%d.RData", results_dir, cond))
  }
  
  # Update progress information with condition completion time 
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