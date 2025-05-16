# Set up libraries
library(modsem)
library(lavaan)
library(covsim)

# replications
rep <- 1

dir.create("sim_results", showWarnings = FALSE)
timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
results_dir <- paste0("sim_results/run_", timestamp)
dir.create(results_dir, showWarnings = FALSE)

source("Simulation_Scripts/GenerateData(2).R")
source("Simulation_Scripts/Simulation_Study_2/Simulation_2/Methods(2).R")
source("Simulation_Scripts/Simulation_Study_2/Simulation_2/Models(2).R")
source("Simulation_Scripts/Simulation_Study_2/Simulation_2/Design(2).R")

if(!"Seed" %in% names(conditions_2)) {
  conditions_2$Seed <- 123:(123 + nrow(conditions_2) - 1) # change latter to avoid dependencies 
}

start_time <- Sys.time()
all_results <- list()

# Process each condition
for(cond in 1:nrow(conditions_2)) {
  condition_start <- Sys.time()
  
  # progress information
  cat("\nRunning condition", cond, "of", nrow(conditions_2), "\n")
  cat("Parameters: N =", conditions_2$N[cond], 
      ", rel =", conditions_2$Rel[cond], 
      ", pop_model =", conditions_2$Population[cond], 
      ", analy_model =", conditions_2$Analysis_model[cond],
      ", distribution =", conditions_2$Distribution[cond], 
      ", exo_method =", conditions_2$Exo_method[cond],
      ", epsilon =", conditions_2$Epsilon[cond],
      ", num_exo_vars =", conditions_2$Num_exo_vars[cond], "\n")
  
  # model number for file naming
  model_num <- as.numeric(gsub(".*\\.", "", conditions_2$Population[cond]))
  
  # Delete this later; just for checking as of now
  pop_model <- tryCatch({
    get(conditions_2$Population[cond])
  }, error = function(e) {
    cat("Failed to get population model:", e$message, "\n")
    return(NULL)
  })
  
  fit_model <- tryCatch({
    get(conditions_2$Analysis_model[cond])
  }, error = function(e) {
    cat("Failed to get fit model:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(pop_model) || is.null(fit_model)) {
    cat("Models not found. Skipping condition.\n")
    next
  }
  
  # proper lengths for this model
  num_vars <- conditions_2$Num_exo_vars[cond]
  skewness_vec <- rep(0, num_vars)
  excesskurtosis_vec <- rep(0, num_vars)
  exo_mean_vec <- rep(0, num_vars)
  
  R2_val <- list("eta4" = conditions_2$R2[cond])
  
  # test dataset first to determine parameters
  cat("Generating test dataset to determine parameters...\n")
  test_data <- tryCatch({
    set.seed(conditions_2$Seed[cond])
    GenerateData(
      model = pop_model,
      N = 100, # small sample for testing
      skewness = skewness_vec,
      excesskurtosis = excesskurtosis_vec,
      exo.mean = exo_mean_vec,
      center.exogenous.latent = TRUE,
      center.exogenous.manifest = FALSE,
      center.lv.dependent = FALSE,
      center.lv.prod = FALSE,
      center.indicators = FALSE,
      distr.exo = conditions_2$Exo_method[cond],
      distr.zeta = "normal",
      distr.epsilon = conditions_2$Epsilon[cond],
      rel = conditions_2$Rel[cond],
      R2 = R2_val,
      add.eta = FALSE
    )
  }, error = function(e) {
    cat("Error generating test data:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(test_data)) {
    cat("Failed to generate test data. Skipping condition.\n")
    next
  }
  
  # Try to run LMS on test data to get parameter names
  test_lms_result <- tryCatch({
    method_analytic(Data = test_data, model.fit = fit_model, method = "lms")
  }, error = function(e) {
    cat("Error running test LMS:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(test_lms_result)) {
    cat("Failed to run test LMS. Skipping condition.\n")
    next
  }
  
  # Extract parameter names from test run
  params <- names(test_lms_result$Estimates)
  n_params <- length(params)
  
  cat("Successfully determined", n_params, "parameters:", paste(params, collapse=", "), "\n")
  
  # results containers
  res <- list(
    lms = array(NA, dim = c(rep, n_params, 3),
                dimnames = list(NULL, params, c("beta", "se", "pval"))),
    qml = array(NA, dim = c(rep, n_params, 3),
                dimnames = list(NULL, params, c("beta", "se", "pval"))),
    uca = array(NA, dim = c(rep, n_params, 3),
                dimnames = list(NULL, params, c("beta", "se", "pval"))),
    sam = matrix(NA, rep, n_params, dimnames = list(NULL, params)),
    timing = data.frame(
      lms = numeric(rep),
      qml = numeric(rep),
      uca = numeric(rep),
      sam = numeric(rep)
    )
  )
  
  # replication sequentially
  successful_reps <- c()
  
  for(i in 1:rep) {
    cat("Processing replication", i, "of", rep, "\n")
    
    # data for this replication
    set.seed(conditions_2$Seed[cond] + i)
    
    Data <- tryCatch({
      GenerateData(
        model = pop_model,
        N = conditions_2$N[cond],
        skewness = skewness_vec,
        excesskurtosis = excesskurtosis_vec,
        exo.mean = exo_mean_vec,
        center.exogenous.latent = TRUE,
        center.exogenous.manifest = FALSE,
        center.lv.dependent = FALSE,
        center.lv.prod = FALSE,
        center.indicators = FALSE,
        distr.exo = conditions_2$Exo_method[cond],
        distr.zeta = "normal",
        distr.epsilon = conditions_2$Epsilon[cond],
        rel = conditions_2$Rel[cond],
        R2 = R2_val,
        add.eta = FALSE
      )
    }, error = function(e) {
      cat("  Error generating data for replication", i, ":", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(Data)) {
      cat("  Skipping replication", i, "due to data generation failure\n")
      next
    }
    
    # LMS method
    start_time_method <- Sys.time()
    lms_result <- tryCatch({
      method_analytic(Data = Data, model.fit = fit_model, method = "lms")
    }, error = function(e) {
      cat("  Error with LMS method for replication", i, ":", e$message, "\n")
      return(NULL)
    })
    res$timing$lms[i] <- as.numeric(difftime(Sys.time(), start_time_method, units = "secs"))
    
    if (!is.null(lms_result)) {
      res$lms[i, , 1] <- lms_result$Estimates[params]
      res$lms[i, , 2] <- lms_result$`Standard Errors`[params]
      res$lms[i, , 3] <- lms_result$`P-values`[params]
      successful_reps <- c(successful_reps, i)
      cat("  LMS successful for replication", i, "\n")
    }
    
    # QML method
    start_time_method <- Sys.time()
    qml_result <- tryCatch({
      method_analytic(Data = Data, model.fit = fit_model, method = "qml")
    }, error = function(e) {
      cat("  Error with QML method for replication", i, ":", e$message, "\n")
      return(NULL)
    })
    res$timing$qml[i] <- as.numeric(difftime(Sys.time(), start_time_method, units = "secs"))
    
    if (!is.null(qml_result)) {
      res$qml[i, , 1] <- qml_result$Estimates[params]
      res$qml[i, , 2] <- qml_result$`Standard Errors`[params]
      res$qml[i, , 3] <- qml_result$`P-values`[params]
      cat("  QML successful for replication", i, "\n")
    }
    
    # UCA method
    start_time_method <- Sys.time()
    uca_result <- tryCatch({
      method_uca(Data = Data, model.fit = fit_model)
    }, error = function(e) {
      cat("  Error with UCA method for replication", i, ":", e$message, "\n")
      return(NULL)
    })
    res$timing$uca[i] <- as.numeric(difftime(Sys.time(), start_time_method, units = "secs"))
    
    if (!is.null(uca_result)) {
      res$uca[i, , 1] <- uca_result$Estimates[params]
      res$uca[i, , 2] <- uca_result$`Standard Errors`[params]
      res$uca[i, , 3] <- uca_result$`P-values`[params]
      cat("  UCA successful for replication", i, "\n")
    }
    
    # SAM method
    start_time_method <- Sys.time()
    sam_result <- tryCatch({
      method_sam(Data = Data, model.fit = fit_model)
    }, error = function(e) {
      cat("  Error with SAM method for replication", i, ":", e$message, "\n")
      return(NULL)
    })
    res$timing$sam[i] <- as.numeric(difftime(Sys.time(), start_time_method, units = "secs"))
    
    if (!is.null(sam_result)) {
      res$sam[i, ] <- sam_result$Estimates[params]
      cat("  SAM successful for replication", i, "\n")
    }
  }
  
  # results for this condition
  if (length(successful_reps) > 0) {
    all_results[[cond]] <- list(
      condition = conditions_2[cond, ],
      results = res,
      successful_reps = successful_reps
    )
    
    condition_filename <- sprintf(
      "%s/model_%d_N%d_Rel%.1f.RData",
      results_dir, 
      model_num,
      conditions_2$N[cond],
      conditions_2$Rel[cond]
    )
    
    tryCatch({
      save(res, file = condition_filename)
      cat("Saved results for condition", cond, "to", condition_filename, "\n")
    }, error = function(e) {
      warning(paste("Failed to save condition results:", e$message))
    })
    
    cat("\nSummary for condition", cond, ":\n")
    cat("Successful replications:", length(successful_reps), "out of", rep, "\n")
    
  }
  
  # condition completion and timing
  condition_time <- difftime(Sys.time(), condition_start, units = "mins")
  cat(sprintf("\nCondition %d completed in %.2f minutes\n", cond, condition_time))
  cat(sprintf("Progress: %.1f%% complete (%d of %d conditions)\n", 
              100 * cond / nrow(conditions_2), cond, nrow(conditions_2)))
}

total_time <- difftime(Sys.time(), start_time, units = "mins")
cat(sprintf("\nTotal simulation time: %.2f minutes (%.2f hours)\n", 
            total_time, total_time/60))

final_filename <- sprintf("%s/complete_results.RData", results_dir)
save(all_results, conditions_2, file = final_filename)
cat("Saved complete results to", final_filename, "\n")