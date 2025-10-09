############################ 1. General Information ############################
# See README file for more information concerning this file. 

# This file contains the code necessary to run the simulation study. 
# It is dependent on the file "Models.RData" where we store the lavaan-based 
# syntax models for generating the data. It is also dependent on the file
# "Methods.R" where we specify the functions for estimating the different
# approaches

# Relevant to re-start after running this script once to default Mer-Twi.
# RNGkind("Mersenne-Twister", "Inversion", "Rejection")

############################### 2. Documentation ################################

#' Simulation Study Parameters and Settings
#' 
#' @param N_REPLICATIONS    Integer. Number of Monte Carlo replications per condition (default = 1000).
#' @param SAMPLE_SIZES      Integer vector. Sample sizes to simulate. Default is c(400, 1000).
#' @param SEED_START        Integer. Starting seed for reproducibility (default = 123)
#' @param USE_ROBUST_SE     Logical. Whether to use robust standard errors for LMS, QML, and DBLCENT (default = FALSE)
#' 
#' @param analysis.model    Character. Lavaan syntax for the analysis model with interaction 
#'                          and quadratic effects. Used for fitting all methods.
#' 
#' @param distributions     List. Distribution specifications with three types:
#'   - normal:    skewness = c(0,0), excesskurtosis = c(0,0), distr.exo = "normal.rIG".
#'   - nonnormal: skewness = c(2,2), excesskurtosis = c(7,7), distr.exo = "nonnormal.rIG".
#'   - uniform:   skewness = c(0,0), excesskurtosis = c(0,0), distr.exo = "unif".
#' 
#' @param conditions        Data.frame. Full factorial design with:
#'   - N:             Sample sizes (400, 1000).
#'   - Rel:           Reliability levels (0.4, 0.6, 0.8).
#'   - Distribution:  Distribution types (normal, nonnormal, uniform).
#'   - Model_Type:    "full" (with interaction/quadratic) or "linear" (without).
#' 
#' @param n_cores          Integer. Number of parallel cores (default = detectCores() - 6).
#' 
#' @return all_results     List. Contains for each condition:
#'   - condition:         Row from conditions data.frame
#'   - results:           List with:
#'     * lms_tables:      Parameter tables from LMS estimation
#'     * qml_tables:      Parameter tables from QML estimation  
#'     * dblcent_tables:  Parameter tables from DBLCENT estimation
#'     * sam_tables:      Parameter tables from SAM estimation
#'     * timing:          Data.frame with computation times for each method
#'     * warnings:        List of warning messages for each method
#'     * observed_r2:     Observed R² values for eta3
#'     * observed_rel:    Matrix of observed reliabilities (9 indicators)
#'     * rng_states:      RNG states for reproducibility
#'   - true_parameters:   True population parameters based on Model_Type
#' 
#' @note Output Files:
#'   - Checkpoint files: Results_Study_1_checkpoint_[n].RData (every 5 conditions)
#'   - Final file:       Results_Study_1_final.RData (all conditions)
#'   - Directory:        Simulations/Study_1/Data/
#' 
#' @note Dependencies:
#'   - Models(1).RData:   Contains all_models with population models.
#'   - Methods.R:         Contains method_analytic(), method_dblcent(), method_sam().
#'   - GenerateData.R:    Contains GenerateData() function.
#' 
#' @note Error Handling:
#'   - Data generation errors: Skip iteration and return NULL.
#'   - Method estimation errors: Store NULL for table, NA for timing, preserve warnings.
#'   - Warnings tracked separately without stopping execution.

############################### 3. Simulation ##################################

library(lavaan); library(modsem); library(doParallel); 
library(doRNG); library(covsim); library(copula); 

# required files
load("Simulations/Study_1/Simulation/Models(1).RData")  # all_models
source("Simulations/Methods.R")  # methods file
source("Simulations/GenerateData.R") # generate data function

# SIMULATION PARAMETERS

N_REPLICATIONS <- 1000
SAMPLE_SIZES   <- c(400, 1000)
SEED_START     <- 123
USE_ROBUST_SE  <- FALSE  # relevant for supplementary materials: TRUE for robust SE

# analysis model (for fitting - no fixed values)
analysis.model <- "
eta1 =~ x1 + x2 + x3
eta2 =~ x4 + x5 + x6
eta3 =~ x7 + x8 + x9

eta3 ~ eta1 + eta2 + eta1:eta2 + eta1:eta1 + eta2:eta2
"

# distribution parameters
distributions <- list(
  normal    = list(skewness = c(0, 0), excesskurtosis = c(0, 0), distr.exo = "normal.rIG"),
  nonnormal = list(skewness = c(2, 2), excesskurtosis = c(7, 7), distr.exo = "nonnormal.rIG"),
  uniform   = list(skewness = c(0, 0), excesskurtosis = c(0, 0), distr.exo = "unif")
)

# output directory
base_dir <- "Simulations/Study_1/Data"
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
file_suffix <- ifelse(USE_ROBUST_SE, "_robustse", "") # _robustse suffix to filename only if using robust SE
results_base <- file.path(base_dir, paste0("Results_Study_1", file_suffix))

# conditions for sim study 1
conditions <- expand.grid(
  N             = SAMPLE_SIZES,
  Rel           = c(0.4, 0.6, 0.8),
  Distribution  = names(distributions),
  Model_Type    = c("full", "linear"),  # linear = no interaction/quadratic terms
  stringsAsFactors = FALSE
)

# model names based on model type and reliability
conditions$model_name <- ifelse(
  conditions$Model_Type == "linear",
  paste0("null_normal_rel", gsub("\\.", "", as.character(conditions$Rel))),
  paste0("normal_rel", gsub("\\.", "", as.character(conditions$Rel)))
)

# track actual distribution used for generation
conditions$generation_distribution <- conditions$Distribution

# track robust SE setting in conditions
conditions$robust_se <- USE_ROBUST_SE

# parallel processing setup
n_cores <- max(1, detectCores() - 6)
cl <- makeCluster(n_cores)
registerDoParallel(cl)
clusterExport(cl, c(
  "GenerateData", "method_analytic", "method_dblcent", "method_sam",
  "analysis.model", "all_models", "distributions", "USE_ROBUST_SE"
))

# main sim study 1 loop
all_results <- list()
start_time  <- Sys.time()

for (cond in 1:nrow(conditions)) {
  cat("\n========================================")
  cat("\nCondition", cond, "of", nrow(conditions))
  cat("\n- Sample size:", conditions$N[cond])
  cat("\n- Target reliability:", conditions$Rel[cond])
  cat("\n- Model type:", conditions$Model_Type[cond])
  cat("\n- Actual distribution:", conditions$Distribution[cond])
  cat("\n- Using model:", conditions$model_name[cond])
  cat("\n- Robust SE:", USE_ROBUST_SE)  # ADDED: Display robust SE setting
  cat("\n========================================\n")
  
  # get the appropriate population model
  population_model <- all_models[[conditions$model_name[cond]]]
  if (is.null(population_model)) {
    cat("WARNING: Model", conditions$model_name[cond], "not found. Skipping...\n")
    next
  }
  
  # distribution parameters for data generation
  dist_params <- distributions[[conditions$Distribution[cond]]]
  
  # push condition-specific into cluster
  clusterExport(cl, c("population_model", "dist_params", "conditions", "cond"), envir = environment())
  
  res <- list(
    lms_tables     = vector("list", N_REPLICATIONS),
    qml_tables     = vector("list", N_REPLICATIONS),
    dblcent_tables = vector("list", N_REPLICATIONS),
    sam_tables     = vector("list", N_REPLICATIONS),
    
    timing = data.frame(
      lms     = numeric(N_REPLICATIONS),
      qml     = numeric(N_REPLICATIONS),
      dblcent = numeric(N_REPLICATIONS),
      sam     = numeric(N_REPLICATIONS)
    ),
    
    warnings = list(
      lms     = vector("list", N_REPLICATIONS),
      qml     = vector("list", N_REPLICATIONS),
      dblcent = vector("list", N_REPLICATIONS),
      sam     = vector("list", N_REPLICATIONS)
    ),
    
    observed_r2  = numeric(N_REPLICATIONS),
    observed_rel = matrix(NA, nrow = N_REPLICATIONS, ncol = 9),
    robust_se_used = USE_ROBUST_SE  # track robust SE setting in results
  )
  
  # run replications in parallel
  parallel_results <- foreach(
    i = 1:N_REPLICATIONS,
    .packages = c("lavaan", "modsem", "covsim", "copula"),
    .errorhandling = "pass",
    .options.RNG = SEED_START + cond * 1000
  ) %dorng% {
    
    # data generation (skip iteration on hard error)
    Data <- try(GenerateData(
      model           = population_model,
      N               = conditions$N[cond],
      skewness        = dist_params$skewness,
      excesskurtosis  = dist_params$excesskurtosis,
      distr.exo       = dist_params$distr.exo,
      distr.epsilon   = "normal",
      distr.zeta      = "normal",
      add.eta         = FALSE,
      return.info     = TRUE
    ), silent = TRUE)
    
    if (inherits(Data, "try-error")) return(NULL)
    
    observed_metrics <- list(
      r2  = attr(Data, "observed_R2")$eta3,
      rel = unlist(attr(Data, "observed_reliabilities"))
    )
    
    data_clean <- as.data.frame(Data)
    attributes(data_clean) <- attributes(data_clean)[c("names", "row.names", "class")]
    
    results <- list(observed_metrics = observed_metrics)
    
    # helper: run a method and capture warnings only
    run_with_warnings <- function(expr) {
      warns <- NULL
      t0 <- Sys.time()
      out <- withCallingHandlers(
        try(expr, silent = TRUE),
        warning = function(w) {
          warns <<- c(warns, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
      if (!inherits(out, "try-error")) {
        list(table = out,
             timing = as.numeric(difftime(Sys.time(), t0, units = "secs")),
             warnings = warns)
      } else {
        # hard error: no table/timing; warnings (if any) already captured
        list(table = NULL, timing = NA_real_, warnings = warns)
      }
    }
    
    # LMS 
    lms_res <- run_with_warnings(
      method_analytic(Data = data_clean, model.fit = analysis.model, 
                      method = "lms", robust.se = USE_ROBUST_SE)
    )
    results$lms_table    <- lms_res$table
    results$lms_timing   <- lms_res$timing
    results$lms_warnings <- lms_res$warnings
    
    # QML 
    qml_res <- run_with_warnings(
      method_analytic(Data = data_clean, model.fit = analysis.model, 
                      method = "qml", robust.se = USE_ROBUST_SE)
    )
    results$qml_table    <- qml_res$table
    results$qml_timing   <- qml_res$timing
    results$qml_warnings <- qml_res$warnings
    
    # UPI/DBLCENT
    dblcent_res <- run_with_warnings(
      method_dblcent(Data = data_clean, model.fit = analysis.model, 
                     robust.se = USE_ROBUST_SE)
    )
    results$dblcent_table    <- dblcent_res$table
    results$dblcent_timing   <- dblcent_res$timing
    results$dblcent_warnings <- dblcent_res$warnings
    
    # LSAM 
    # Note: SAM doesn't have robust.se option in the same way
    sam_res <- run_with_warnings(
      method_sam(Data = data_clean, model.fit = analysis.model)
    )
    results$sam_table    <- sam_res$table
    results$sam_timing   <- sam_res$timing
    results$sam_warnings <- sam_res$warnings
    
    results
  }
  
  # RNG states (for reproducibility if needed)
  rng_states_for_condition <- attr(parallel_results, "rng")
  
  # collect results (skip NULLs and foreach error objects) 
  for (i in 1:N_REPLICATIONS) {
    iter <- parallel_results[[i]]
    if (is.null(iter) || inherits(iter, "error")) next
    
    if (!is.null(iter$observed_metrics)) {
      res$observed_r2[i]    <- iter$observed_metrics$r2
      res$observed_rel[i, ] <- iter$observed_metrics$rel[1:9]
    }
    
    # LMS
    if (!is.null(iter$lms_table))  res$lms_tables[[i]]   <- iter$lms_table
    if (!is.null(iter$lms_timing)) res$timing$lms[i]     <- iter$lms_timing
    if (!is.null(iter$lms_warnings)) res$warnings$lms[[i]] <- iter$lms_warnings
    
    # QML
    if (!is.null(iter$qml_table))  res$qml_tables[[i]]   <- iter$qml_table
    if (!is.null(iter$qml_timing)) res$timing$qml[i]     <- iter$qml_timing
    if (!is.null(iter$qml_warnings)) res$warnings$qml[[i]] <- iter$qml_warnings
    
    # DBLCENT/UPI
    if (!is.null(iter$dblcent_table))  res$dblcent_tables[[i]] <- iter$dblcent_table
    if (!is.null(iter$dblcent_timing)) res$timing$dblcent[i]   <- iter$dblcent_timing
    if (!is.null(iter$dblcent_warnings)) res$warnings$dblcent[[i]] <- iter$dblcent_warnings
    
    # LSAM
    if (!is.null(iter$sam_table))  res$sam_tables[[i]]   <- iter$sam_table
    if (!is.null(iter$sam_timing)) res$timing$sam[i]     <- iter$sam_timing
    if (!is.null(iter$sam_warnings)) res$warnings$sam[[i]] <- iter$sam_warnings
  }
  
  res$rng_states <- rng_states_for_condition
  
  # summary during sim study 1
  cat("\nObserved metrics across replications:")
  cat("\n- Mean R²:", mean(res$observed_r2, na.rm = TRUE))
  cat("\n- Mean reliabilities:", round(colMeans(res$observed_rel, na.rm = TRUE), 3))
  
  cat("\n\nWarnings encountered:")
  cat("\n- LMS:",     sum(lengths(res$warnings$lms) > 0),     "iterations with warnings")
  cat("\n- QML:",     sum(lengths(res$warnings$qml) > 0),     "iterations with warnings")
  cat("\n- DBLCENT:", sum(lengths(res$warnings$dblcent) > 0), "iterations with warnings")
  cat("\n- SAM:",     sum(lengths(res$warnings$sam) > 0),     "iterations with warnings")
  cat("\n")
  
  # store condition results
  all_results[[cond]] <- list(
    condition = conditions[cond, ],
    results   = res,
    true_parameters = if (conditions$Model_Type[cond] == "linear") {
      c(0.316, 0.316, 0, 0, 0)   # linear model: interaction & quadratics 0
    } else {
      c(0.316, 0.316, 0.139, 0.101, 0.101)  # full model
    }
  )
  
  # checkpoint every 5 conditions (and final)
  if (cond %% 5 == 0 || cond == nrow(conditions)) {
    save(all_results, conditions, file = sprintf("%s_checkpoint_%d.RData", results_base, cond))
  }
  
  gc()
}

stopCluster(cl)

# save with appropriate suffix if or no robust SE
save(all_results, conditions, file = paste0(results_base, "_final.RData"))

total_time <- difftime(Sys.time(), start_time, units = "hours")
cat(sprintf("\n\nSimulation completed in %.2f hours\n", total_time))
cat(sprintf("\nResults saved with suffix: %s\n", file_suffix))