############################ 1. General Information ############################

# See README file for more information concerning this file. 

# This file contains the code necessary to run the simulation study. 
# It is dependent on the file "Models.RData" where we store the lavaan-based 
# syntax models for generating the data. It is also dependent on the file
# "Methods.R" where we specify the functions for estimating the different
# approaches

# Relevant to re-start after running this script once to default Mer-Twi.
# RNGkind("Mersenne-Twister", "Inversion", "Rejection")

############################### 2. Documentation ###############################

#' Simulation 1 Parameters and Settings
#' 
#' @param N_REPLICATIONS    Integer. Number of Monte Carlo replications per condition (default = 1000).
#' @param SAMPLE_SIZES      Integer vector. Sample sizes to simulate. Default is c(400, 1000).
#' @param SEED_START        Integer. Starting seed for reproducibility (default = 123)
#' @param USE_ROBUST_SE     Logical. Whether to use robust standard errors for LMS, QML, and UPI (default = FALSE)
#' 
#' @param analysis.model    Character. Lavaan syntax for the analysis model with interaction 
#'                          and quadratic effects. Used for fitting all methods.
#' 
#' @param conditions  Data.frame. Full factorial design with:
#'   - N:             Sample sizes (400, 1000).
#'   - Rel:           Reliability levels (0.4, 0.6, 0.8).
#'   - Distr_Exo:     Distribution for exogenous latent variables (normal, nonnormal, uniform).
#'   - Distr_Epsilon: Distribution for measurement errors (normal, exp.rate1).
#'   - Distr_Zeta:    Distribution for structural disturbances (normal, exp.rate1).
#'   - Model_Type:    "full" (with interaction/quadratic) or "linear" (without).
#' 
#' @param n_cores         Integer. Number of parallel cores (default = detectCores() - 4).
#' 
#' @return all_results    List. Contains for each condition:
#'   - condition:         Row from conditions data.frame
#'   - results:           List with:
#'     * lms_tables:      Parameter tables from LMS estimation
#'     * qml_tables:      Parameter tables from QML estimation  
#'     * upi_tables:      Parameter tables from UPI estimation
#'     * lsam_tables:     Parameter tables from LSAM estimation
#'     * timing:          Data.frame with computation times for each method
#'     * warnings:        List of warning messages for each method
#'     * observed_r2:     Observed r squared values for eta3
#'     * observed_rel:    Matrix of observed reliabilities (9 indicators)
#'     * rng_states:      RNG states for reproducibility
#'   - true_parameters:   True population parameters based on Model_Type
#' 
#' @note Output Files:
#'   - Checkpoint files: Data_Study_1_checkpoint_[n].RData (every 5 conditions)
#'   - Final file:       Data_Study_1_final.RData (all conditions)
#'   - Directory:        Simulations/Study_1/Data/
#' 
#' @note Files dependencies:
#'   - Models(1).RData:   Contains all_models with population models.
#'   - Methods.R:         Contains method_analytic(), method_upi(), method_lsam().
#'   - GenerateData.R:    Contains GenerateData() function (VITA-based generation).
#' 
#' @note Error Handling:
#'   - Data generation errors: Skip iteration and return NULL.
#'   - Method estimation errors: Store NULL for table, NA for timing, preserve warnings.
#'   - Warnings tracked separately without stopping execution.
#'   
#' @note Dependencies:
#'   Required packages: lavaan, modsem, doParallel, doRNG, covsim, and rvinecopulib
#'   

############################### 3. Simulation ##################################

# required files
source("Simulations/Study_1/Simulation/Models(1).R")  # all_models
source("Simulations/Methods.R")  # methods file
source("Simulations/GenerateData.R") 

# simulation parameters
N_REPLICATIONS <- 1000
SAMPLE_SIZES   <- c(400, 1000)
SEED_START     <- 123
USE_ROBUST_SE  <- FALSE  # relevant for supplemental materials: TRUE for robust SE

# analysis model (for fitting - no fixed values)
analysis.model <- "
eta1 =~ x1 + x2 + x3
eta2 =~ x4 + x5 + x6
eta3 =~ x7 + x8 + x9

eta3 ~ eta1 + eta2 + eta1:eta2 + eta1:eta1 + eta2:eta2
"

# output directory
base_dir <- "Simulations/Study_1/Data"
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
file_suffix <- ifelse(USE_ROBUST_SE, "_robustse", "") 
results_base <- file.path(base_dir, paste0("Data_Study_1", file_suffix))

# full factorial design: 2 × 3 × 3 × 2 × 2 × 2 = 144 conditions
conditions <- expand.grid(
  N              = SAMPLE_SIZES,
  Rel            = c(0.4, 0.6, 0.8),
  Distr_Exo      = c("normal", "nonnormal", "uniform"),
  Distr_Epsilon  = c("normal", "exp.rate1"),
  Distr_Zeta     = c("normal", "exp.rate1"),
  Model_Type     = c("full", "linear"),
  stringsAsFactors = FALSE
)

# model names based on model type and reliability
conditions$model_name <- ifelse(
  conditions$Model_Type == "linear",
  paste0("null_normal_rel", gsub("\\.", "", as.character(conditions$Rel))),
  paste0("normal_rel", gsub("\\.", "", as.character(conditions$Rel)))
)

# just to track robust SE setting in conditions
conditions$robust_se <- USE_ROBUST_SE

cat("Total conditions:", nrow(conditions), "\n")

# parallel processing setup
n_cores <- max(1, detectCores() - 4)
cl <- makeCluster(n_cores)
registerDoParallel(cl)
clusterExport(cl, c(
  "GenerateData", "method_analytic", "method_upi", "method_lsam",
  "analysis.model", "all_models", "USE_ROBUST_SE"
))

# main simulation study 1 loop
all_results <- list()
start_time  <- Sys.time()

for (cond in 1:nrow(conditions)) {
  cat("\nCondition", cond, "of", nrow(conditions))
  cat("\n- Sample size:", conditions$N[cond])
  cat("\n- Target reliability:", conditions$Rel[cond])
  cat("\n- Model type:", conditions$Model_Type[cond])
  cat("\n- Exogenous distribution:", conditions$Distr_Exo[cond])
  cat("\n- Epsilon distribution:", conditions$Distr_Epsilon[cond])
  cat("\n- Zeta distribution:", conditions$Distr_Zeta[cond])
  cat("\n- Using model:", conditions$model_name[cond])
  cat("\n- Robust SE:", USE_ROBUST_SE)
  
  # get the appropriate population model
  population_model <- all_models[[conditions$model_name[cond]]]
  if (is.null(population_model)) {
    cat("\nWARNING: Model", conditions$model_name[cond], "not found. Skipping...\n")
    next
  }
  
  # build distribution parameters dynamically
  dist_params <- list(
    distr.exo       = conditions$Distr_Exo[cond],
    distr.epsilon   = conditions$Distr_Epsilon[cond],
    distr.zeta      = conditions$Distr_Zeta[cond],
    nonnormal.shape = if (conditions$Distr_Exo[cond] == "nonnormal") c(1, 1) else NULL,
    nonnormal.rate  = if (conditions$Distr_Exo[cond] == "nonnormal") c(1, 1) else NULL
  )
  
  # push condition-specific into cluster
  clusterExport(cl, c("population_model", "dist_params", "conditions", "cond"), envir = environment())
  
  res <- list(
    lms_tables  = vector("list", N_REPLICATIONS),
    qml_tables  = vector("list", N_REPLICATIONS),
    upi_tables  = vector("list", N_REPLICATIONS),
    lsam_tables = vector("list", N_REPLICATIONS),
    
    timing = data.frame(
      lms  = numeric(N_REPLICATIONS),
      qml  = numeric(N_REPLICATIONS),
      upi  = numeric(N_REPLICATIONS),
      lsam = numeric(N_REPLICATIONS)
    ),
    
    warnings = list(
      lms  = vector("list", N_REPLICATIONS),
      qml  = vector("list", N_REPLICATIONS),
      upi  = vector("list", N_REPLICATIONS),
      lsam = vector("list", N_REPLICATIONS)
    ),
    
    observed_r2  = numeric(N_REPLICATIONS),
    observed_rel = matrix(NA, nrow = N_REPLICATIONS, ncol = 9),
    robust_se_used = USE_ROBUST_SE
  )
  
  # run replications in parallel
  parallel_results <- foreach(
    i = 1:N_REPLICATIONS,
    .packages = c("lavaan", "modsem", "covsim", "rvinecopulib"),  
    .errorhandling = "pass",
    .options.RNG = SEED_START + cond * 1000
  ) %dorng% {
    
    # data generation 
    Data <- try(GenerateData(
      model           = population_model,
      N               = conditions$N[cond],
      distr.exo       = dist_params$distr.exo,
      nonnormal.shape = dist_params$nonnormal.shape,
      nonnormal.rate  = dist_params$nonnormal.rate,
      distr.epsilon   = dist_params$distr.epsilon,
      distr.zeta      = dist_params$distr.zeta,
      add.eta         = FALSE, # needs to be false; otherwise triggers error in estimation
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
    
    # UPI
    upi_res <- run_with_warnings(
      method_upi(Data = data_clean, model.fit = analysis.model, 
                 robust.se = USE_ROBUST_SE)
    )
    results$upi_table    <- upi_res$table
    results$upi_timing   <- upi_res$timing
    results$upi_warnings <- upi_res$warnings
    
    # LSAM 
    lsam_res <- run_with_warnings(
      method_lsam(Data = data_clean, model.fit = analysis.model, joint = TRUE)
    )
    results$lsam_table    <- lsam_res$table
    results$lsam_timing   <- lsam_res$timing
    results$lsam_warnings <- lsam_res$warnings
    
    results
  }
  
  # RNG states (for reproducibility if needed for any specific iteration)
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
    
    # UPI
    if (!is.null(iter$upi_table))  res$upi_tables[[i]] <- iter$upi_table
    if (!is.null(iter$upi_timing)) res$timing$upi[i]   <- iter$upi_timing
    if (!is.null(iter$upi_warnings)) res$warnings$upi[[i]] <- iter$upi_warnings
    
    # LSAM
    if (!is.null(iter$lsam_table))  res$lsam_tables[[i]]   <- iter$lsam_table
    if (!is.null(iter$lsam_timing)) res$timing$lsam[i]     <- iter$lsam_timing
    if (!is.null(iter$lsam_warnings)) res$warnings$lsam[[i]] <- iter$lsam_warnings
  }
  
  res$rng_states <- rng_states_for_condition
  
  # summary during simulation study 1
  cat("\n\nObserved metrics across replications:")
  cat("\nMean R^2:", mean(res$observed_r2, na.rm = TRUE))
  cat("\nMean reliabilities:", round(colMeans(res$observed_rel, na.rm = TRUE), 3))
  
  cat("\n\nWarnings encountered:")
  cat("\nLMS:", sum(lengths(res$warnings$lms) > 0),   "iterations with warnings")
  cat("\nQML:", sum(lengths(res$warnings$qml) > 0),   "iterations with warnings")
  cat("\nUPI:", sum(lengths(res$warnings$upi) > 0),   "iterations with warnings")
  cat("\nLSAM:", sum(lengths(res$warnings$lsam) > 0), "iterations with warnings")
  cat("\n")
  
  # store condition results 
  # note we save to disk immediately and clear from memory 
  condition_result <- list(
    condition = conditions[cond, ],
    results   = res,
    true_parameters = if (conditions$Model_Type[cond] == "linear") {
      c(0.316, 0.316, 0, 0, 0)   # linear model: interaction & quadratics = 0
    } else {
      c(0.316, 0.316, 0.139, 0.101, 0.101)  # full model
    }
  )
  
  # save condition to its own file
  save(condition_result, file = sprintf("%s_condition_%03d.RData", results_base, cond))
  cat("\nCondition", cond, "saved to disk\n")
  # clear from memory
  rm(condition_result, res, parallel_results, rng_states_for_condition)
  
  gc()
}

stopCluster(cl)

total_time <- difftime(Sys.time(), start_time, units = "hours")
cat(sprintf("\nSimulation completed in %.2f hours\n", total_time))

# note to myself (15/1/26): maybe move this to results file.
# once the simulation lop is done, one may:
# load conditions back into all_results for final save
all_results <- vector("list", nrow(conditions))
for (cond in 1:nrow(conditions)) {
  load(sprintf("%s_condition_%03d.RData", results_base, cond))
  all_results[[cond]] <- condition_result
}

# if needed all results from all conditions combined:
# save with appropriate suffix if or no robust SE
save(all_results, conditions, file = paste0(results_base, "_final.RData"))