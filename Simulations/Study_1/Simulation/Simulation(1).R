############################ 1. General Information ############################

# See README file for more information concerning this file. 

# This file contains the code necessary to run the simulation study. 
# It is dependent on the file "Models.RData" where we store the lavaan-based 
# syntax models for generating the data. It is also dependent on the file
# "Methods.R" where we specify the functions for estimating the different
# approaches

# Relevant to re-start after running this script once to default MT
# RNGkind("Mersenne-Twister", "Inversion", "Rejection")

############################### 2. Simulation ##################################

library(lavaan); library(modsem); library(doParallel); library(doRNG)
library(covsim); library(copula); library(stringr)

# required files
load("Models.RData")  # all_models
source("Methods(1).R")  # methods file

# SIMULATION PARAMETERS

N_REPLICATIONS <- 2
SAMPLE_SIZES   <- c(400, 1000)
SEED_START     <- 123

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

# --- output dirs ---
dir.create("sim_results", showWarnings = FALSE)
timestamp    <- format(Sys.time(), "%Y%m%d_%H%M")
results_dir  <- paste0("sim_results/run_", timestamp)
dir.create(results_dir, showWarnings = FALSE)

# --- conditions grid ---
conditions <- expand.grid(
  N             = SAMPLE_SIZES,
  Rel           = c(0.4, 0.6, 0.8),
  Distribution  = names(distributions),
  Model_Type    = c("alternative", "null"),  # null = linear
  stringsAsFactors = FALSE
)

# model names based on model type and reliability
conditions$model_name <- ifelse(
  conditions$Model_Type == "null",
  paste0("null_normal_rel", gsub("\\.", "", as.character(conditions$Rel))),
  paste0("normal_rel", gsub("\\.", "", as.character(conditions$Rel)))
)

# track actual distribution used for generation
conditions$generation_distribution <- conditions$Distribution

# --- parallel setup ---
n_cores <- max(1, detectCores() - 6)
cl <- makeCluster(n_cores)
registerDoParallel(cl)
clusterExport(cl, c(
  "GenerateData", "method_analytic", "method_dblcent", "method_sam",
  "analysis.model", "all_models", "distributions"
))

# --- main loop ---
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
  
  # --- initialize results (warnings only, no error storage) ---
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
    
    # warning tracking - lists of character vectors (one per replication)
    warnings = list(
      lms     = vector("list", N_REPLICATIONS),
      qml     = vector("list", N_REPLICATIONS),
      dblcent = vector("list", N_REPLICATIONS),
      sam     = vector("list", N_REPLICATIONS)
    ),
    
    # observed metrics
    observed_r2  = numeric(N_REPLICATIONS),
    observed_rel = matrix(NA, nrow = N_REPLICATIONS, ncol = 9)
  )
  
  # --- run replications in parallel ---
  parallel_results <- foreach(
    i = 1:N_REPLICATIONS,
    .packages = c("lavaan", "modsem", "covsim", "copula"),
    .errorhandling = "pass",
    .options.RNG = SEED_START + cond * 1000
  ) %dorng% {
    
    #----- data generation (skip iteration on hard error) -----
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
    
    # observed metrics
    observed_metrics <- list(
      r2  = attr(Data, "observed_R2")$eta3,
      rel = unlist(attr(Data, "observed_reliabilities"))
    )
    
    # strip attributes for methods
    data_clean <- as.data.frame(Data)
    attributes(data_clean) <- attributes(data_clean)[c("names", "row.names", "class")]
    
    results <- list(observed_metrics = observed_metrics)
    
    #----- helper: run a method & capture warnings only -----
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
    
    #=============== LMS ===============#
    lms_res <- run_with_warnings(
      method_analytic(Data = data_clean, model.fit = analysis.model, method = "lms")
    )
    results$lms_table    <- lms_res$table
    results$lms_timing   <- lms_res$timing
    results$lms_warnings <- lms_res$warnings
    
    #=============== QML ===============#
    qml_res <- run_with_warnings(
      method_analytic(Data = data_clean, model.fit = analysis.model, method = "qml")
    )
    results$qml_table    <- qml_res$table
    results$qml_timing   <- qml_res$timing
    results$qml_warnings <- qml_res$warnings
    
    #=============== UCA/DBLCENT ===============#
    dblcent_res <- run_with_warnings(
      method_dblcent(Data = data_clean, model.fit = analysis.model)
    )
    results$dblcent_table    <- dblcent_res$table
    results$dblcent_timing   <- dblcent_res$timing
    results$dblcent_warnings <- dblcent_res$warnings
    
    #=============== SAM ===============#
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
  
  # --- collect results (skip NULLs and foreach error objects) ---
  for (i in 1:N_REPLICATIONS) {
    iter <- parallel_results[[i]]
    if (is.null(iter) || inherits(iter, "error")) next
    
    # observed metrics
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
    
    # DBLCENT
    if (!is.null(iter$dblcent_table))  res$dblcent_tables[[i]] <- iter$dblcent_table
    if (!is.null(iter$dblcent_timing)) res$timing$dblcent[i]   <- iter$dblcent_timing
    if (!is.null(iter$dblcent_warnings)) res$warnings$dblcent[[i]] <- iter$dblcent_warnings
    
    # SAM
    if (!is.null(iter$sam_table))  res$sam_tables[[i]]   <- iter$sam_table
    if (!is.null(iter$sam_timing)) res$timing$sam[i]     <- iter$sam_timing
    if (!is.null(iter$sam_warnings)) res$warnings$sam[[i]] <- iter$sam_warnings
  }
  
  res$rng_states <- rng_states_for_condition
  
  # --- summary (warnings only) ---
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
    true_parameters = if (conditions$Model_Type[cond] == "null") {
      c(0.316, 0.316, 0, 0, 0)   # null (linear): interaction & quadratics 0
    } else {
      c(0.316, 0.316, 0.139, 0.101, 0.101)  # full model
    }
  )
  
  # checkpoint every 5 conditions (and final)
  if (cond %% 5 == 0 || cond == nrow(conditions)) {
    save(all_results, conditions, file = sprintf("%s/checkpoint_%d.RData", results_dir, cond))
  }
  
  gc()
}

stopCluster(cl)

#save(all_results, conditions, file = paste0(results_dir, "/final_results.RData"))

total_time <- difftime(Sys.time(), start_time, units = "hours")
cat(sprintf("\n\nSimulation completed in %.2f hours\n", total_time))