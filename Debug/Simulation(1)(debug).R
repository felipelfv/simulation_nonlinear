# Other scripts needed
source("Simulation_Study_1/Simulation/GenerateData(1).R")
source("Simulation_Study_1/Simulation/Methods(1).R")
source("Simulation_Study_1/Simulation/Models(1).R")
source("Simulation_Study_1/Simulation/Design(1).R")

# data hash:
# (see: https://www.rdocumentation.org/packages/digest/versions/0.6.37/topics/digest)
if(!require(digest)) {
  install.packages("digest")
  library(digest)
}

# for debug: function to get a data snapshot for debugging
get_data_snapshot <- function(data, n_rows=5, n_cols=3) {
  if(is.null(data) || !is.data.frame(data)) return("NULL or invalid data")
  
  n_cols <- min(n_cols, ncol(data))
  n_rows <- min(n_rows, nrow(data))
  # subset
  subset_data <- data[1:n_rows, 1:n_cols, drop=FALSE]
  
  # summary stats for first column
  if(ncol(data) > 0 && nrow(data) > 0) {
    col1 <- data[[1]]
    col_stats <- sprintf("Col1 stats: min=%.3f, max=%.3f, mean=%.3f", 
                         min(col1, na.rm=TRUE), max(col1, na.rm=TRUE), mean(col1, na.rm=TRUE))
  } else {
    col_stats <- "No data for column stats"
  }
  
  list(
    data = subset_data,
    stats = col_stats
  )
}

# parallel backend
n_cores <- detectCores() - 2
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Debug: parallel setup info
cat("Debug: Using", n_cores, "cores for parallel processing\n")

clusterExport(cl, c("GenerateData", "method_analytic", "method_uca", "method_sam", "population.interaction.model",
                    "population.linear.model", "population.full.model", "fit.interaction.model", "fit.full.model"))

# Debug: Also export digest function and snapshot function
clusterExport(cl, c("digest", "get_data_snapshot"))

dir.create("sim_results", showWarnings = FALSE)
timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
results_dir <- paste0("sim_results/run_", timestamp)
dir.create(results_dir, showWarnings = FALSE)

# Debug: Create debug log file
debug_log <- file(paste0(results_dir, "/debug_log.txt"), "w")
cat("Debug log: ", as.character(Sys.time()), "\n", file=debug_log)
cat("Configuration: cores =", n_cores, "\n", file=debug_log)

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
  
  # Debug: Create data hash and snapshot for verification
  data_hash <- NULL
  data_snapshot <- NULL
  
  if(!inherits(Data, "try-error")) {
    data_hash <- digest(Data)
    data_snapshot <- get_data_snapshot(Data)
    
    # Debug: Print data snapshot in worker
    cat(sprintf("Debug [Worker %d]: Data snapshot for condition %d, rep %d:\n", i, cond, i))
    print(data_snapshot$data)
    cat(sprintf("Debug [Worker %d]: %s\n", i, data_snapshot$stats))
    
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
  
  # Return both LMS results and data with hash and snapshot
  list(data = if(!inherits(Data, "try-error")) Data else NULL,
       lms_res = lms_res,
       data_hash = data_hash,
       data_snapshot = data_snapshot)  # Debug: data snapshot
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
  
  # Debug: condition details to file
  cat("\nDEBUG: Starting condition", cond, "at", as.character(Sys.time()), "\n", file=debug_log)
  cat("Parameters: N =", conditions$N[cond], 
      ", rel =", conditions$Rel[cond], 
      ", pop_model =", conditions$Population[cond], 
      ", analy_model =", conditions$Analysis_model[cond],
      ", distribution =", conditions$Distribution[cond], 
      ", exo_method =", conditions$Exo_method[cond],
      ", epsilon =", conditions$Epsilon[cond], "\n", file=debug_log)
  
  skewness <- rep(ifelse(conditions$Distribution[cond] == "normal", 0, 2), 2)
  excesskurtosis <- rep(ifelse(conditions$Distribution[cond] == "normal", 0, 7), 2)
  
  # N of parameters as a function of the analysis model
  n_params <- ifelse(conditions$Analysis_model[cond] == "fit.full.model", 5, 3)
  
  # Debug: number of parameters
  cat("Debug: Using n_params =", n_params, "for", conditions$Analysis_model[cond], "\n")
  
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
    ),
    # Debug: data tracking fields
    data_hashes = character(rep),
    data_snapshots = vector("list", rep)
  )
  
  # Parallel processing ONLY for LMS 
  cat("\nRunning LMS in parallel\n")
  lms_parallel_results <- foreach(i = seq_len(rep), 
                                  .packages = c("modsem", "lavaan", "covsim", "digest"), 
                                  .errorhandling = "pass",
                                  .options.RNG = conditions$Seed[cond]) %dorng% {
                                    process_lms_only(i, cond, n_params, skewness, excesskurtosis)
                                  }
  
  # Debug: log summary of parallel results
  cat("Debug: LMS parallel processing completed\n", file=debug_log)
  cat(sprintf("Debug: Got %d results from parallel processing\n", length(lms_parallel_results)))
  
  # RNG states for this condition
  rng_states <- attr(lms_parallel_results, "rng")
  
  # Debug: counter for successful/failed replications
  success_count <- 0
  failed_count <- 0
  
  # Process other methods sequentially
  cat("\nRunning QML, UCA, and SAM sequentially\n")
  for(i in 1:rep) {
    # Debug: replication processing
    cat(sprintf("\n=== Debug: Processing replication %d/%d for condition %d ===\n", i, rep, cond))
    
    # skip if data generation failed
    if(is.null(lms_parallel_results[[i]]$data)) {
      failed_count <- failed_count + 1
      cat(sprintf("Debug: Skipping rep %d due to NULL data\n", i))
      next
    }
    
    # Store data hash and snapshot
    res$data_hashes[i] <- lms_parallel_results[[i]]$data_hash
    res$data_snapshots[[i]] <- lms_parallel_results[[i]]$data_snapshot
    
    # Debug: LMS data snapshot from parallel processing
    cat("Debug: LMS DATA SNAPSHOT from parallel processing:\n")
    print(lms_parallel_results[[i]]$data_snapshot$data)
    cat(lms_parallel_results[[i]]$data_snapshot$stats, "\n")
    
    # LMS results
    res$lms[i, , 1] <- lms_parallel_results[[i]]$lms_res$estimates
    res$lms[i, , 2] <- lms_parallel_results[[i]]$lms_res$se
    res$lms[i, , 3] <- lms_parallel_results[[i]]$lms_res$pvals
    res$timing$lms[i] <- lms_parallel_results[[i]]$lms_res$timing
    
    # Debug: LMS results were stored correctly
    cat(sprintf("Debug: LMS est[1] = %.3f, se[1] = %.3f, pval[1] = %.3f\n", 
                res$lms[i, 1, 1], res$lms[i, 1, 2], res$lms[i, 1, 3]))

    cat(sprintf("Debug: LMS est[2] = %.3f, se[2] = %.3f, pval[2] = %.3f\n", 
                res$lms[i, 2, 1], res$lms[i, 2, 2], res$lms[i, 2, 3]))
    
    # Dataset from parallel results
    Data <- lms_parallel_results[[i]]$data
    
    # Debug: snapshot of the data to verify consistency
    current_snapshot <- get_data_snapshot(Data)
    cat("Debug: CURRENT DATA SNAPSHOT before QML/UCA/SAM:\n")
    print(current_snapshot$data)
    cat(current_snapshot$stats, "\n")
    
    # Debug: data integrity using hash
    current_hash <- digest(Data)
    if(current_hash != res$data_hashes[i]) {
      warning(sprintf("DATA INTEGRITY ERROR: Hash mismatch for rep %d", i))
      cat(sprintf("DATA INTEGRITY ERROR: Hash mismatch for rep %d\n", i), file=debug_log)
      cat(sprintf("  Original: %s\n  Current: %s\n", 
                  res$data_hashes[i], current_hash), file=debug_log)
    } else {
      cat(sprintf("Debug: Data hash verified for rep %d: %s\n", i, current_hash))
    }
    
    analysis_model <- get(conditions$Analysis_model[cond])
    
    # QML method
    start_time_method <- Sys.time()
    result <- suppressWarnings(try(method_analytic(Data = Data, model.fit = analysis_model, method = "qml")))
    res$timing$qml[i] <- as.numeric(difftime(Sys.time(), start_time_method, units = "secs"))
    
    # Debug: data remains consistent after QML
    qml_snapshot <- get_data_snapshot(Data)
    cat("Debug: DATA SNAPSHOT after QML:\n")
    print(qml_snapshot$data)
    cat(qml_snapshot$stats, "\n")
    
    if(!inherits(result, "try-error")) {
      res$qml[i, , 1] <- result$Estimates
      res$qml[i, , 2] <- result$`Standard Errors`
      res$qml[i, , 3] <- result$`P-values`
      
      # Debug: Log QML results
      cat(sprintf("Debug: QML successful, first est: %.3f\n", result$Estimates[1]))
    } else {
      cat(sprintf("Deubg: QML failed: %s\n", as.character(result)))
      cat(sprintf("Debug: QML failed for rep %d: %s\n", i, as.character(result)), file=debug_log)
    }
    
    # UCA method
    start_time_method <- Sys.time()
    result <- suppressWarnings(try(method_uca(Data = Data, model.fit = analysis_model)))
    res$timing$uca[i] <- as.numeric(difftime(Sys.time(), start_time_method, units = "secs"))
    
    # Debug: data remains consistent after UCA
    uca_snapshot <- get_data_snapshot(Data)
    cat("Debug: DATA SNAPSHOT after UCA:\n")
    print(uca_snapshot$data)
    cat(uca_snapshot$stats, "\n")
    
    if(!inherits(result, "try-error")) {
      res$uca[i, , 1] <- result$Estimates
      res$uca[i, , 2] <- result$`Standard Errors`
      res$uca[i, , 3] <- result$`P-values`
      
      # Debug: UCA results
      cat(sprintf("Debug: UCA successful, first est: %.3f\n", result$Estimates[1]))
    } else {
      cat(sprintf("Debug: UCA failed: %s\n", as.character(result)))
      cat(sprintf("Debug: UCA failed for rep %d: %s\n", i, as.character(result)), file=debug_log)
    }
    
    # SAM method
    start_time_method <- Sys.time()
    result <- suppressWarnings(try(method_sam(Data = Data, model.fit = analysis_model)))
    res$timing$sam[i] <- as.numeric(difftime(Sys.time(), start_time_method, units = "secs"))
    
    # Debug: data remains consistent after SAM
    sam_snapshot <- get_data_snapshot(Data)
    cat("Debug: DATA SNAPSHOT after SAM:\n")
    print(sam_snapshot$data)
    cat(sam_snapshot$stats, "\n")
    
    if(!inherits(result, "try-error")) {
      # As of now, still no standard errors; no SE and p-value
      res$sam[i, ] <- result$Estimates
      
      # Debug: SAM results
      cat(sprintf("Debug: SAM successful, first est: %.3f\n", result$Estimates[1]))
    } else {
      cat(sprintf("Debug: SAM failed: %s\n", as.character(result)))
      cat(sprintf("Debug: SAM failed for rep %d: %s\n", i, as.character(result)), file=debug_log)
    }
    
    # Debug: data snapshots across methods
    cat("\nDebug: COMPARING DATA SNAPSHOTS ACROSS METHODS:\n")
    if(identical(current_snapshot$data, qml_snapshot$data) && 
       identical(qml_snapshot$data, uca_snapshot$data) && 
       identical(uca_snapshot$data, sam_snapshot$data)) {
      cat("✓ DATA CONSISTENT across all methods for this replication\n")
    } else {
      cat("✗ DATA INCONSISTENCY detected between methods!\n")
      
      # Check each method pair
      if(!identical(current_snapshot$data, qml_snapshot$data)) 
        cat("- Difference between initial data and QML data\n")
      if(!identical(current_snapshot$data, uca_snapshot$data)) 
        cat("- Difference between initial data and UCA data\n")
      if(!identical(current_snapshot$data, sam_snapshot$data)) 
        cat("- Difference between initial data and SAM data\n")
      if(!identical(qml_snapshot$data, uca_snapshot$data)) 
        cat("- Difference between QML data and UCA data\n")
      if(!identical(qml_snapshot$data, sam_snapshot$data)) 
        cat("- Difference between QML data and SAM data\n")
      if(!identical(uca_snapshot$data, sam_snapshot$data)) 
        cat("- Difference between UCA data and SAM data\n")
      
      # Log this severe issue to debug file
      cat("\nDATA INCONSISTENCY detected in replication", i, "\n", file=debug_log)
    }
    
    success_count <- success_count + 1
  }
  
  # Debug: success/failure for this condition
  cat(sprintf("\nDebug: Condition %d summary: %d successful, %d failed replications\n", 
              cond, success_count, failed_count))
  cat(sprintf("Debug: Condition %d summary: %d successful, %d failed replications\n", 
              cond, success_count, failed_count), file=debug_log)
  
  # Debug: Check array dimensions and NA patterns
  cat("\nDebg: Results arrays check:\n")
  cat("- LMS dimension:", paste(dim(res$lms), collapse="x"), "\n")
  cat("- QML dimension:", paste(dim(res$qml), collapse="x"), "\n")
  cat("- UCA dimension:", paste(dim(res$uca), collapse="x"), "\n")
  cat("- SAM dimension:", paste(dim(res$sam), collapse="x"), "\n\n")
  
  lms_na_pct <- sum(is.na(res$lms)) / length(res$lms) * 100
  qml_na_pct <- sum(is.na(res$qml)) / length(res$qml) * 100
  uca_na_pct <- sum(is.na(res$uca)) / length(res$uca) * 100
  sam_na_pct <- sum(is.na(res$sam)) / length(res$sam) * 100
  
  cat(sprintf("NA percentages: LMS=%.1f%%, QML=%.1f%%, UCA=%.1f%%, SAM=%.1f%%\n", 
              lms_na_pct, qml_na_pct, uca_na_pct, sam_na_pct))
  
  all_results[[cond]] <- list(
    condition = conditions[cond, ],
    results = res,
    rng_states = rng_states  # rng states
  )
  
  # Results for this condition
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
    cat(sprintf("Debug: Saved results to %s\n", condition_filename))
  }, error = function(e) {
    warning(paste("Failed to save condition results:", e$message))
    cat("Debug: Failed to save results:", e$message, "\n", file=debug_log)
  })
  
  # All results until 10th condition
  if(cond %% 10 == 0 || cond == min(10, nrow(conditions))) {
    tryCatch({
      save(all_results, conditions, 
           file = sprintf("%s/all_results_upto_condition_%d.RData", results_dir, cond))
      cat(sprintf("Debug: Saved all results up to condition %d\n", cond))
    }, error = function(e) {
      warning(paste("Failed to save checkpoint:", e$message))
      cat("Debug: Failed to save checkpoint:", e$message, "\n", file=debug_log)
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
  
  # Debug: condition completion
  cat(sprintf("Debug: Condition %d completed in %.2f minutes\n", 
              cond, condition_time), file=debug_log)
}

stopCluster(cl)

# Debug: Close the log file
cat("Debug log ended: ", as.character(Sys.time()), "\n", file=debug_log)
close(debug_log)