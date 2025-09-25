############################ 1. General Information ############################

# See README file for more information concerning this file. 

# This file contains the code necessary to run the simulation study 2. 
# It is dependent on the file "Models.RData" where we store the lavaan-based 
# syntax models for generating the data. It is also dependent on the file
# "Methods.R" where we specify the functions for estimating the different
# approaches

# Relevant to re-start after running this script once to default Mer-Twi.
# RNGkind("Mersenne-Twister", "Inversion", "Rejection")

############################### 2. Simulation ##################################

library(lavaan); library(modsem); library(covsim)
library(doParallel); library(doRNG); library(copula)

load("Simulations/Study_2/Simulation/Models(2).RData")  # all_models
source("Simulations/Methods.R")  # methods file
source("Simulations/GenerateData.R") # generate data function

# output directory
base_dir <- "Simulations/Study_2/Data"
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
results_base <- file.path(base_dir, "Results_Study_2")

# analysis model
analysis.model <- "
# Measurement model
eta1 =~ x1 + x2 + x3
eta2 =~ x4 + x5 + x6
eta3 =~ x7 + x8 + x9
eta4 =~ x10 + x11 + x12
eta5 =~ x13 + x14 + x15
eta6 =~ x16 + x17 + x18

# Structural model
eta4 ~ eta1 + eta2 + eta3 + eta1:eta2 + eta1:eta1 
eta5 ~ eta4 + eta1 + eta2 + eta3 + eta2:eta4 + eta2:eta2
eta6 ~ eta5 + eta1 + eta2 + eta3 + eta3:eta5 + eta3:eta3
"

# SIMULATION PARAMETERS

N_REPLICATIONS <- 50
SAMPLE_SIZES <- c(400, 1000)
RELIABILITIES <- c(0.4, 0.6, 0.8)
SEED_START <- 123

# distribution parameters
distributions <- list(
  normal = list(
    skewness = c(0, 0, 0), 
    excesskurtosis = c(0, 0, 0), 
    distr.exo = "normal.rIG"
  ),
  nonnormal = list(
    skewness = c(2, 2, 2), 
    excesskurtosis = c(7, 7, 7), 
    distr.exo = "nonnormal.rIG"
  ),
  uniform = list(
    skewness = c(0, 0, 0), 
    excesskurtosis = c(0, 0, 0), 
    distr.exo = "unif"
  )
)

# conditions

conditions <- expand.grid(
  N = SAMPLE_SIZES,
  Rel = RELIABILITIES,
  Distribution = names(distributions),
  Model_Type = c("full", "linear"),  
  stringsAsFactors = FALSE
)

conditions$model_name <- ifelse(
  conditions$Model_Type == "linear",
  paste0("null_model_rel", gsub("\\.", "", as.character(conditions$Rel))),
  paste0("normal_rel", gsub("\\.", "", as.character(conditions$Rel)))
)

# SETUP PARALLEL PROCESSING
n_cores <- detectCores() - 2
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# cluster export to include all models
clusterExport(cl, c("GenerateData", "method_sam", "method_analytic", "method_dblcent",
                    "analysis.model", "all_models", "distributions"))

# MAIN SIMULATION 

all_results <- list()
start_time <- Sys.time()

for (cond in 1:nrow(conditions)) {
  cat("\n", paste(rep("=", 60), collapse = ""))
  cat("\nCondition", cond, "of", nrow(conditions))
  cat("\n- Sample size:", conditions$N[cond])
  cat("\n- Reliability:", conditions$Rel[cond])
  cat("\n- Distribution:", conditions$Distribution[cond])
  cat("\n- Model type:", conditions$Model_Type[cond])  
  cat("\n- Using model:", conditions$model_name[cond]) 
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  
  # population model from all_models based on model_name
  population_model <- all_models[[conditions$model_name[cond]]]
  
  if (is.null(population_model)) {
    cat("WARNING: Model", conditions$model_name[cond], "not found. Skipping...\n")
    next
  }
  
  dist_params <- distributions[[conditions$Distribution[cond]]]
  
  # condition-specific variables to cluster
  clusterExport(
    cl,
    c("population_model", "dist_params", "conditions", "cond"),
    envir = environment()
  )
  
  # --- initialize results (warnings only, no error storage) ---
  res <- list(
    sam_tables     = vector("list", N_REPLICATIONS),
    qml_tables     = vector("list", N_REPLICATIONS),
    dblcent_tables = vector("list", N_REPLICATIONS),
    
    timing = data.frame(
      sam     = numeric(N_REPLICATIONS),
      qml     = numeric(N_REPLICATIONS),
      dblcent = numeric(N_REPLICATIONS)
    ),
    
    # warning tracking - lists of character vectors (one per replication)
    warnings = list(
      sam     = vector("list", N_REPLICATIONS),
      qml     = vector("list", N_REPLICATIONS),
      dblcent = vector("list", N_REPLICATIONS)
    ),
    
    observed_r2  = matrix(NA, nrow = N_REPLICATIONS, ncol = 3),  # eta4, eta5, eta6
    observed_rel = matrix(NA, nrow = N_REPLICATIONS, ncol = 6)   # 6 latent variables
  )
  
  # replications in parallel
  parallel_results <- foreach(
    i = 1:N_REPLICATIONS,
    .packages = c("lavaan", "modsem", "covsim", "copula"),
    .errorhandling = "pass",
    .options.RNG = SEED_START + cond * 10000
  ) %dorng% {
    
    #----- data generation (skip iteration on hard error) -----
    Data <- try(GenerateData(
      model          = population_model,
      N              = conditions$N[cond],
      skewness       = dist_params$skewness,
      excesskurtosis = dist_params$excesskurtosis,
      distr.exo      = dist_params$distr.exo,
      distr.epsilon  = "normal",
      distr.zeta     = "normal",
      add.eta        = FALSE,
      return.info    = TRUE
    ), silent = TRUE)
    
    if (inherits(Data, "try-error")) return(NULL)
    
    # observed metrics
    observed_metrics <- list(
      r2  = c(attr(Data, "observed_R2")$eta4,
              attr(Data, "observed_R2")$eta5,
              attr(Data, "observed_R2")$eta6),
      rel = sapply(attr(Data, "observed_reliabilities"), mean)
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
    
    #=============== SAM ===============#
    sam_res <- run_with_warnings(
      method_sam(Data = data_clean, model.fit = analysis.model)
    )
    results$sam_table    <- sam_res$table
    results$sam_timing   <- sam_res$timing
    results$sam_warnings <- sam_res$warnings
    
    #=============== QML ===============#
    # filter out the specific QML bias warning about exogenous/endogenous interactions
    qml_res <- run_with_warnings(
      method_analytic(Data = data_clean, model.fit = analysis.model, method = "qml"),
      filter_pattern = "Interactions between exogenous and (endogenous|enodgenous).*QML.*approach.*biased"
      # enodgenous because there was a typo in the message from modsem
    )
    results$qml_table    <- qml_res$table
    results$qml_timing   <- qml_res$timing
    results$qml_warnings <- qml_res$warnings
    
    #=============== DBLCENT ===============#
    dblcent_res <- run_with_warnings(
      method_dblcent(Data = data_clean, model.fit = analysis.model)
    )
    results$dblcent_table    <- dblcent_res$table
    results$dblcent_timing   <- dblcent_res$timing
    results$dblcent_warnings <- dblcent_res$warnings
    
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
      res$observed_r2[i, ]  <- iter$observed_metrics$r2
      res$observed_rel[i, ] <- iter$observed_metrics$rel
    }
    
    # SAM
    if (!is.null(iter$sam_table))  res$sam_tables[[i]]   <- iter$sam_table
    if (!is.null(iter$sam_timing)) res$timing$sam[i]     <- iter$sam_timing
    if (!is.null(iter$sam_warnings)) res$warnings$sam[[i]] <- iter$sam_warnings
    
    # QML
    if (!is.null(iter$qml_table))  res$qml_tables[[i]]   <- iter$qml_table
    if (!is.null(iter$qml_timing)) res$timing$qml[i]     <- iter$qml_timing
    if (!is.null(iter$qml_warnings)) res$warnings$qml[[i]] <- iter$qml_warnings
    
    # DBLCENT
    if (!is.null(iter$dblcent_table))  res$dblcent_tables[[i]] <- iter$dblcent_table
    if (!is.null(iter$dblcent_timing)) res$timing$dblcent[i]   <- iter$dblcent_timing
    if (!is.null(iter$dblcent_warnings)) res$warnings$dblcent[[i]] <- iter$dblcent_warnings
  }
  
  res$rng_states <- rng_states_for_condition
  
  # --- summary (warnings only) ---
  cat("\nObserved metrics across replications:")
  cat("\n- Mean R² (eta4, eta5, eta6):", round(colMeans(res$observed_r2, na.rm = TRUE), 3))
  cat("\n- Mean reliabilities:", round(colMeans(res$observed_rel, na.rm = TRUE), 3))
  
  cat("\n\nWarnings encountered:")
  cat("\n- SAM:",     sum(lengths(res$warnings$sam) > 0),     "iterations with warnings")
  cat("\n- QML:",     sum(lengths(res$warnings$qml) > 0),     "iterations with warnings")
  cat("\n- DBLCENT:", sum(lengths(res$warnings$dblcent) > 0), "iterations with warnings")
  cat("\n")
  
  # store condition results
  all_results[[cond]] <- list(
    condition   = conditions[cond, ],
    results     = res,
    true_parameters = if (conditions$Model_Type[cond] == "linear") {
      list(
        # main effects only (same as full)
        eta4_eta1 = 0.21, eta4_eta2 = 0.21, eta4_eta3 = 0.21,
        eta5_eta4 = 0.18, eta5_eta1 = 0.18, eta5_eta2 = 0.18, eta5_eta3 = 0.18,
        eta6_eta5 = 0.15, eta6_eta1 = 0.15, eta6_eta2 = 0.15, eta6_eta3 = 0.15,
        # interactions and quadratics set to 0 for linear model
        eta4_eta1eta2 = 0, eta4_eta1eta1 = 0, 
        eta5_eta2eta4 = 0, eta5_eta2eta2 = 0,
        eta6_eta3eta5 = 0, eta6_eta3eta3 = 0
      )
    } else {
      list(
        # main effects
        eta4_eta1 = 0.21, eta4_eta2 = 0.21, eta4_eta3 = 0.21,
        eta5_eta4 = 0.18, eta5_eta1 = 0.18, eta5_eta2 = 0.18, eta5_eta3 = 0.18,
        eta6_eta5 = 0.15, eta6_eta1 = 0.15, eta6_eta2 = 0.15, eta6_eta3 = 0.15,
        # interactions and quadratics (non-zero for full)
        eta4_eta1eta2 = 0.13, eta4_eta1eta1 = 0.09, 
        eta5_eta2eta4 = 0.11, eta5_eta2eta2 = 0.09,
        eta6_eta3eta5 = 0.10, eta6_eta3eta3 = 0.09
      )
    }
  )
  
  # checkpoint every 5 conditions (and final)
  if (cond %% 5 == 0 || cond == nrow(conditions)) {
    save(all_results, conditions, file = sprintf("%s_checkpoint_%d.RData", results_base, cond))
  }
  
  gc()  
}

stopCluster(cl)

save(all_results, conditions, file = paste0(results_base, "_final.RData"))
total_time <- difftime(Sys.time(), start_time, units = "hours")
cat(sprintf("\n\nSimulation completed in %.2f hours\n", total_time))

# POST-SIMULATION ANALYSIS OF ERRORS AND WARNINGS

# categorize warnings
categorize_warnings <- function(warnings) {
  list(
    quadrature = sum(grepl("adaptive quadrature", warnings, ignore.case = TRUE)),
    null_values = sum(grepl("is NULL", warnings)),
    std_errors = sum(grepl("Standard errors.*could not", warnings)),
    large_variances = sum(grepl("factor 1000", warnings)),
    vcov_not_pd = sum(grepl("not.*positive definite", warnings)),
    identification = sum(grepl("not identified", warnings))
  )
}

for(cond_idx in seq_along(all_results)) {
  cond_info <- all_results[[cond_idx]]$condition
  res <- all_results[[cond_idx]]$results
  
  # total warnings by type for each method
  sam_all_warnings <- unlist(res$warnings$sam)
  qml_all_warnings <- unlist(res$warnings$qml)
  dblcent_all_warnings <- unlist(res$warnings$dblcent)
  
  if(length(c(sam_all_warnings, qml_all_warnings, dblcent_all_warnings)) > 0) {
    cat(sprintf("\nCondition %d (N=%d, Rel=%.1f, Dist=%s, Model=%s):\n", 
                cond_idx, cond_info$N, cond_info$Rel, cond_info$Distribution, 
                cond_info$Model_Type)) 
    
    if(length(sam_all_warnings) > 0) {
      sam_cats <- categorize_warnings(sam_all_warnings)
      cat("  SAM warnings:", paste(names(sam_cats)[sam_cats > 0], 
                                   "=", sam_cats[sam_cats > 0], collapse=", "), "\n")
    }
    
    if(length(qml_all_warnings) > 0) {
      qml_cats <- categorize_warnings(qml_all_warnings)
      cat("  QML warnings:", paste(names(qml_cats)[qml_cats > 0], 
                                   "=", qml_cats[qml_cats > 0], collapse=", "), "\n")
    }
    
    if(length(dblcent_all_warnings) > 0) {
      dblcent_cats <- categorize_warnings(dblcent_all_warnings)
      cat("  DBLCENT warnings:", paste(names(dblcent_cats)[dblcent_cats > 0], 
                                       "=", dblcent_cats[dblcent_cats > 0], collapse=", "), "\n")
    }
  }
}