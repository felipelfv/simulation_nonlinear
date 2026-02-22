############################ 1. General Information ############################

# See README file for more information concerning this file. 

# This file contains the code necessary to run the simulation study 2. 
# It is dependent on the file "Models(2).RData" where we store the lavaan-based 
# syntax models for generating the data. It is also dependent on the file
# "Methods.R" where we specify the functions for estimating the different
# approaches

# Relevant to re-start after running this script once to default Mer-Twi.
# RNGkind("Mersenne-Twister", "Inversion", "Rejection")

############################### 2. Documentation ###############################

#' Simulation 2 Parameters and Settings
#' 
#' @param N_REPLICATIONS    Integer. Number of Monte Carlo replications per condition (default = 1000).
#' @param SAMPLE_SIZES      Integer vector. Sample sizes to simulate. Default is c(400, 1000).
#' @param RELIABILITIES     Numeric vector. Reliability levels (0.4, 0.6, 0.8).
#' @param SEED_START        Integer. Starting seed for reproducibility (default = 123)
#' @param USE_ROBUST_SE     Logical. Whether to use robust standard errors for QML and UPI (default = FALSE)
#' 
#' @param analysis.model    Character. Lavaan syntax for the analysis model with 5 factors.

#' @param conditions        Data.frame. Full factorial design with:
#'   - N:             Sample sizes (400, 1000).
#'   - Rel:           Reliability levels (0.4, 0.6, 0.8).
#'   - Distr_Exo:     Distribution for exogenous latent variables (normal, nonnormal, uniform).
#'   - Distr_Epsilon: Distribution for measurement errors (normal, exp.rate1).
#'   - Distr_Zeta:    Distribution for structural disturbances (normal, exp.rate1).
#'   - Model_Type:    "full" (with interaction/quadratic) or "linear" (without).
#'   - model_name:    Generated names for population models ("normal_rel[X]" or "null_model_rel[X]").
#'   - robust_se:     Logical indicator of robust SE setting for this run.
#' 
#' @param n_cores          Integer. Number of parallel cores (default = detectCores() - 4).
#' 
#' @return all_results     List. Contains for each condition:
#'   - condition:         Row from conditions data.frame
#'   - results:           List with:
#'     * lsam_tables:     Parameter tables from LSAM estimation
#'     * qml_tables:      Parameter tables from QML estimation  
#'     * upi_tables:      Parameter tables from UPI estimation
#'     * timing:          Data.frame with computation times for each method
#'     * warnings:        List of warning messages for each method (QML filters bias warnings)
#'     * observed_r2:     Matrix of observed r squared values for eta4 and eta5 (N_REPLICATIONS x 2)
#'     * observed_rel:    Matrix of observed reliabilities for 5 factors (N_REPLICATIONS x 5)
#'     * robust_se_used:  Logical indicator of robust SE setting
#'     * rng_states:      RNG states for reproducibility
#'   - true_parameters:   List of true population parameters based on Model_Type:
#'     * Intercepts:      eta4_intercept (0.1), eta5_intercept (0.1)
#'     * Main effects:    eta4 predictors (0.20 each), eta5 predictors (0.16 each)
#'     * Interactions:    eta4_eta1eta2 (0.11), eta4_eta1eta3 (0.11), eta5_eta1eta4 (0.08), eta5_eta2eta4 (0.08)
#'     * Quadratics:      eta4_eta1eta1 (0.08), eta4_eta2eta2 (0.08), eta5_eta1eta1 (0.06), eta5_eta3eta3 (0.06)
#'     * For linear model: all interactions and quadratics set to 0
#' 
#' @note Output Files:
#'   - Individual files: Data_Study_2[_robustse]_condition_[NNN].RData (one per condition)
#'   - Final file:       Data_Study_2[_robustse]_final.RData (all conditions combined)
#'   - Directory:        Simulations/Study_2/Data/
#' 
#' @note File dependencies:
#'   - Models(2).RData:   Contains all_models with population models for Study 2.
#'   - Methods.R:         Contains method_analytic(), method_upi(), method_lsam().
#'   - GenerateData.R:    Contains GenerateData() function (VITA-based generation).
#' 
#' @note Error Handling:
#'   - Data generation errors: Skip iteration and return NULL.
#'   - Method estimation errors: Store NULL for table, NA for timing, preserve warnings.
#'   - QML warnings about exogenous-endogenous interactions are filtered (known bias).
#'   - Warnings tracked separately without stopping execution.
#'   
#' @note Dependencies:
#'   Required packages: lavaan, modsem, doParallel, doRNG, covsim, and rvinecopulib
#'    

############################### 3. Simulation ##################################

# required packages
library(lavaan)
library(modsem)
library(doParallel)
library(doRNG)
library(covsim)
library(rvinecopulib)

source("Simulations/Study_2/Simulation/Models(2).R")
source("Simulations/Methods.R")  # methods file
source("Simulations/GenerateData.R")

# SIMULATION PARAMETERS

N_REPLICATIONS <- 1000
SAMPLE_SIZES <- c(400, 1000)
RELIABILITIES <- c(0.4, 0.6, 0.8)
SEED_START <- 123
USE_ROBUST_SE  <- FALSE  # relevant for supplemental materials: TRUE for robust SE

# output directory
base_dir <- "Simulations/Study_2/Data"
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
file_suffix <- ifelse(USE_ROBUST_SE, "_robustse", "") 
results_base <- file.path(base_dir, paste0("Data_Study_2", file_suffix))

# analysis model for 5 factors with new interaction structure
analysis.model <- "
# Measurement model
eta1 =~ x1 + x2 + x3
eta2 =~ x4 + x5 + x6
eta3 =~ x7 + x8 + x9
eta4 =~ x10 + x11 + x12
eta5 =~ x13 + x14 + x15

# Structural model
eta4 ~ eta1 + eta2 + eta3 + eta1:eta2 + eta1:eta3 + eta1:eta1 + eta2:eta2
eta5 ~ eta4 + eta1 + eta2 + eta3 + eta1:eta4 + eta2:eta4 + eta1:eta1 + eta3:eta3
"

# full factorial design: 2 × 3 × 3 × 2 × 2 × 2 = 144 conditions
conditions <- expand.grid(
  N              = SAMPLE_SIZES,
  Rel            = RELIABILITIES,
  Distr_Exo      = c("normal", "nonnormal", "uniform"),
  Distr_Epsilon  = c("normal", "exp.rate1"),
  Distr_Zeta     = c("normal", "exp.rate1"),
  Model_Type     = c("full", "linear"),
  stringsAsFactors = FALSE
)

conditions$model_name <- ifelse(
  conditions$Model_Type == "linear",
  paste0("null_model_rel", gsub("\\.", "", as.character(conditions$Rel))),
  paste0("normal_rel", gsub("\\.", "", as.character(conditions$Rel)))
)

# track robust SE setting in conditions
conditions$robust_se <- USE_ROBUST_SE

cat("Total conditions:", nrow(conditions), "\n")

# setup parallel processing
n_cores <- max(1, detectCores() - 4)
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# cluster export
clusterExport(cl, c("GenerateData", "method_lsam", "method_analytic", "method_upi",
                    "analysis.model", "all_models", "USE_ROBUST_SE"))

# main simulation study 2 loop
start_time <- Sys.time()

for (cond in 1:nrow(conditions)) {
  cat("\n==============================================================")
  cat("\nCondition", cond, "of", nrow(conditions))
  cat("\n==============================================================")
  cat("\n- Sample size:", conditions$N[cond])
  cat("\n- Reliability:", conditions$Rel[cond])
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
    nonnormal.shape = if (conditions$Distr_Exo[cond] == "nonnormal") c(1, 1, 1) else NULL,
    nonnormal.rate  = if (conditions$Distr_Exo[cond] == "nonnormal") c(1, 1, 1) else NULL
  )
  
  # condition-specific variables to cluster
  clusterExport(
    cl,
    c("population_model", "dist_params", "conditions", "cond"),
    envir = environment()
  )
  
  # initialize results (warnings only, no error storage) 
  res <- list(
    lsam_tables = vector("list", N_REPLICATIONS),
    qml_tables  = vector("list", N_REPLICATIONS),
    upi_tables  = vector("list", N_REPLICATIONS),
    
    timing = data.frame(
      lsam = numeric(N_REPLICATIONS),
      qml  = numeric(N_REPLICATIONS),
      upi  = numeric(N_REPLICATIONS)
    ),
    
    # warning tracking - lists of character vectors (one per replication)
    warnings = list(
      lsam = vector("list", N_REPLICATIONS),
      qml  = vector("list", N_REPLICATIONS),
      upi  = vector("list", N_REPLICATIONS)
    ),
    
    observed_r2  = matrix(NA, nrow = N_REPLICATIONS, ncol = 2),  # eta4, eta5 
    observed_rel = matrix(NA, nrow = N_REPLICATIONS, ncol = 5),  # 5 latent variables 
    robust_se_used = USE_ROBUST_SE
  )
  
  # replications in parallel
  parallel_results <- foreach(
    i = 1:N_REPLICATIONS,
    .packages = c("lavaan", "modsem", "covsim", "rvinecopulib"),  
    .errorhandling = "pass",
    .options.RNG = SEED_START + cond * 10000
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
      add.eta         = FALSE,
      return.info     = TRUE
    ), silent = TRUE)
    
    if (inherits(Data, "try-error")) return(NULL)
    
    # observed metrics 
    observed_metrics <- list(
      r2  = c(attr(Data, "observed_R2")$eta4,
              attr(Data, "observed_R2")$eta5),
      rel = sapply(attr(Data, "observed_reliabilities")[1:5], mean)  # only 5 factors 
    )
    
    # clean data for methods
    data_clean <- as.data.frame(Data)
    attributes(data_clean) <- attributes(data_clean)[c("names", "row.names", "class")]
    
    results <- list(observed_metrics = observed_metrics)
    
    # run a method and capture warnings only 
    run_with_warnings <- function(expr, filter_pattern = NULL) {
      warns <- NULL
      t0 <- Sys.time()
      out <- withCallingHandlers(
        try(expr, silent = TRUE),
        warning = function(w) {
          msg <- conditionMessage(w)
          # add warning if it doesn't match the filter pattern (or no filter specified)
          if (is.null(filter_pattern) || !grepl(filter_pattern, msg, ignore.case = TRUE)) {
            warns <<- c(warns, msg)
          }
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
    
    # LSAM 
    lsam_res <- run_with_warnings(
      method_lsam(Data = data_clean, model.fit = analysis.model, joint = TRUE)
    )
    results$lsam_table    <- lsam_res$table
    results$lsam_timing   <- lsam_res$timing
    results$lsam_warnings <- lsam_res$warnings
    
    # QML 
    # filter out the specific QML bias warning about exogenous/endogenous interactions
    qml_res <- run_with_warnings(
      method_analytic(Data = data_clean, model.fit = analysis.model, 
                      method = "qml", robust.se = USE_ROBUST_SE),
      filter_pattern = "Interactions between exogenous and (endogenous|enodgenous).*QML.*approach.*biased"
      # enodgenous because there was a typo in the message from modsem
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
    
    results
  }
  
  # RNG states (for reproducibility if needed)
  rng_states_for_condition <- attr(parallel_results, "rng")
  
  # collect results (skip NULLs and foreach error objects)
  for (i in 1:N_REPLICATIONS) {
    iter <- parallel_results[[i]]
    if (is.null(iter) || inherits(iter, "error")) next
    
    # observed metrics
    if (!is.null(iter$observed_metrics)) {
      res$observed_r2[i, ]  <- iter$observed_metrics$r2
      res$observed_rel[i, ] <- iter$observed_metrics$rel
    }
    
    # LSAM
    if (!is.null(iter$lsam_table))  res$lsam_tables[[i]]   <- iter$lsam_table
    if (!is.null(iter$lsam_timing)) res$timing$lsam[i]     <- iter$lsam_timing
    if (!is.null(iter$lsam_warnings)) res$warnings$lsam[[i]] <- iter$lsam_warnings
    
    # QML
    if (!is.null(iter$qml_table))  res$qml_tables[[i]]   <- iter$qml_table
    if (!is.null(iter$qml_timing)) res$timing$qml[i]     <- iter$qml_timing
    if (!is.null(iter$qml_warnings)) res$warnings$qml[[i]] <- iter$qml_warnings
    
    # UPI
    if (!is.null(iter$upi_table))  res$upi_tables[[i]] <- iter$upi_table
    if (!is.null(iter$upi_timing)) res$timing$upi[i]   <- iter$upi_timing
    if (!is.null(iter$upi_warnings)) res$warnings$upi[[i]] <- iter$upi_warnings
  }
  
  res$rng_states <- rng_states_for_condition
  
  # summary (warnings only)
  cat("\n\nObserved metrics across replications:")
  cat("\nMean R² (eta4, eta5):", round(colMeans(res$observed_r2, na.rm = TRUE), 3))
  cat("\nMean reliabilities:", round(colMeans(res$observed_rel, na.rm = TRUE), 3))
  
  cat("\n\nWarnings encountered:")
  cat("\nLSAM:", sum(lengths(res$warnings$lsam) > 0), "iterations with warnings")
  cat("\nQML:",  sum(lengths(res$warnings$qml) > 0),  "iterations with warnings")
  cat("\nUPI:",  sum(lengths(res$warnings$upi) > 0),  "iterations with warnings")
  cat("\n")
  
  # store condition results with updated true parameters
  condition_result <- list(
    condition   = conditions[cond, ],
    results     = res,
    true_parameters = if (conditions$Model_Type[cond] == "linear") {
      list(
        # intercepts
        eta4_intercept = 0.1,
        eta5_intercept = 0.1,
        # main effects only
        eta4_eta1 = 0.20, eta4_eta2 = 0.20, eta4_eta3 = 0.20,
        eta5_eta4 = 0.16, eta5_eta1 = 0.16, eta5_eta2 = 0.16, eta5_eta3 = 0.16,
        # interactions and quadratics set to 0 for linear model
        eta4_eta1eta2 = 0, eta4_eta1eta3 = 0, eta4_eta1eta1 = 0, eta4_eta2eta2 = 0,
        eta5_eta1eta4 = 0, eta5_eta2eta4 = 0, eta5_eta1eta1 = 0, eta5_eta3eta3 = 0
      )
    } else {
      list(
        # intercepts  
        eta4_intercept = 0.1,
        eta5_intercept = 0.1,
        # main effects
        eta4_eta1 = 0.20, eta4_eta2 = 0.20, eta4_eta3 = 0.20,
        eta5_eta4 = 0.16, eta5_eta1 = 0.16, eta5_eta2 = 0.16, eta5_eta3 = 0.16,
        # interactions and quadratics (non-zero for full)
        eta4_eta1eta2 = 0.11, eta4_eta1eta3 = 0.11, 
        eta4_eta1eta1 = 0.08, eta4_eta2eta2 = 0.08,
        eta5_eta1eta4 = 0.08, eta5_eta2eta4 = 0.08, 
        eta5_eta1eta1 = 0.06, eta5_eta3eta3 = 0.06
      )
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

# reload all conditions into all_results for final combined save
all_results <- vector("list", nrow(conditions))
for (cond in 1:nrow(conditions)) {
  load(sprintf("%s_condition_%03d.RData", results_base, cond))
  all_results[[cond]] <- condition_result
}

# save all results from all conditions combined
save(all_results, conditions, file = paste0(results_base, "_final.RData"))