# PARALLEL SIMULATION WITH SAM, QML, DBLCENT

library(lavaan)
library(modsem)
library(parallel)
library(doParallel)
library(doRNG)
library(copula)


# SIMULATION PARAMETERS

N_REPLICATIONS <- 1500
SAMPLE_SIZES <- c(400, 1000)
RELIABILITIES <- c(0.8, 0.6, 0.4)
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

# CONDITIONS - Now includes Model_Type

conditions <- expand.grid(
  N = SAMPLE_SIZES,
  Rel = RELIABILITIES,
  Distribution = names(distributions),
  Model_Type = c("alternative", "null"),  
  stringsAsFactors = FALSE
)

# Updated model naming to account for null vs alternative
conditions$model_name <- ifelse(
  conditions$Model_Type == "null",
  paste0("null_model_rel", gsub("\\.", "", as.character(conditions$Rel))),
  paste0("normal_rel", gsub("\\.", "", as.character(conditions$Rel)))
)

# SETUP PARALLEL PROCESSING
n_cores <- detectCores() - 2
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Updated cluster export to include all models
clusterExport(cl, c("GenerateData", "method_sam", "method_qml", "method_dblcent",
                    "analysis_model", "all_models", "distributions"))

# MAIN SIMULATION — store messages

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
  
  res <- list(
    sam_tables     = vector("list", N_REPLICATIONS),
    qml_tables     = vector("list", N_REPLICATIONS),
    dblcent_tables = vector("list", N_REPLICATIONS),
    
    timing = data.frame(
      sam     = numeric(N_REPLICATIONS),
      qml     = numeric(N_REPLICATIONS),
      dblcent = numeric(N_REPLICATIONS)
    ),
    
    # ERROR STORAGE: character vectors ("" if none)
    errors = data.frame(
      sam     = character(N_REPLICATIONS),
      qml     = character(N_REPLICATIONS),
      dblcent = character(N_REPLICATIONS),
      stringsAsFactors = FALSE
    ),
    
    # WARNING STORAGE: list of character vectors per replication
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
    
    # generate data
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
    
    # SAM 
    t0 <- Sys.time()
    sam_warnings <- NULL
    
    sam_table <- withCallingHandlers(
      try(method_sam(Data = data_clean, model.fit = analysis_model), silent = TRUE),
      warning = function(w) {
        sam_warnings <<- c(sam_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    
    if (!inherits(sam_table, "try-error")) {
      results$sam_table  <- sam_table
      results$sam_timing <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    } else {
      results$sam_error <- as.character(sam_table)
    }
    results$sam_warnings <- sam_warnings
    
    # QML 
    t0 <- Sys.time()
    qml_warnings <- NULL
    
    qml_table <- withCallingHandlers(
      try(method_qml(Data = data_clean, model.fit = analysis_model), silent = TRUE),
      warning = function(w) {
        qml_warnings <<- c(qml_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    
    if (!inherits(qml_table, "try-error")) {
      results$qml_table  <- qml_table
      results$qml_timing <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    } else {
      results$qml_error <- as.character(qml_table)
    }
    results$qml_warnings <- qml_warnings
    
    # DBLCENT 
    t0 <- Sys.time()
    dblcent_warnings <- NULL
    
    dblcent_table <- withCallingHandlers(
      try(method_dblcent(Data = data_clean, model.fit = analysis_model), silent = TRUE),
      warning = function(w) {
        dblcent_warnings <<- c(dblcent_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    
    if (!inherits(dblcent_table, "try-error")) {
      results$dblcent_table  <- dblcent_table
      results$dblcent_timing <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    } else {
      results$dblcent_error <- as.character(dblcent_table)
    }
    results$dblcent_warnings <- dblcent_warnings
    
    results
  }
  
  # store RNG states
  rng_states_for_condition <- attr(parallel_results, "rng")
  
  # process parallel results
  convergence_count <- list(sam = 0, qml = 0, dblcent = 0)
  
  for (i in 1:N_REPLICATIONS) {
    if (!is.null(parallel_results[[i]]) && !inherits(parallel_results[[i]], "error")) {
      
      if (!is.null(parallel_results[[i]]$observed_metrics)) {
        res$observed_r2[i, ]  <- parallel_results[[i]]$observed_metrics$r2
        res$observed_rel[i, ] <- parallel_results[[i]]$observed_metrics$rel
      }
      
      # SAM results, errors, warnings
      if (!is.null(parallel_results[[i]]$sam_table)) {
        res$sam_tables[[i]] <- parallel_results[[i]]$sam_table
        res$timing$sam[i]   <- parallel_results[[i]]$sam_timing
        convergence_count$sam <- convergence_count$sam + 1
      }
      if (!is.null(parallel_results[[i]]$sam_error)) {
        res$errors$sam[i] <- parallel_results[[i]]$sam_error
      }
      if (!is.null(parallel_results[[i]]$sam_warnings)) {
        res$warnings$sam[[i]] <- parallel_results[[i]]$sam_warnings
      }
      
      # QML results, errors, warnings
      if (!is.null(parallel_results[[i]]$qml_table)) {
        res$qml_tables[[i]] <- parallel_results[[i]]$qml_table
        res$timing$qml[i]   <- parallel_results[[i]]$qml_timing
        convergence_count$qml <- convergence_count$qml + 1
      }
      if (!is.null(parallel_results[[i]]$qml_error)) {
        res$errors$qml[i] <- parallel_results[[i]]$qml_error
      }
      if (!is.null(parallel_results[[i]]$qml_warnings)) {
        res$warnings$qml[[i]] <- parallel_results[[i]]$qml_warnings
      }
      
      # DBLCENT results, errors, warnings
      if (!is.null(parallel_results[[i]]$dblcent_table)) {
        res$dblcent_tables[[i]] <- parallel_results[[i]]$dblcent_table
        res$timing$dblcent[i]   <- parallel_results[[i]]$dblcent_timing
        convergence_count$dblcent <- convergence_count$dblcent + 1
      }
      if (!is.null(parallel_results[[i]]$dblcent_error)) {
        res$errors$dblcent[i] <- parallel_results[[i]]$dblcent_error
      }
      if (!is.null(parallel_results[[i]]$dblcent_warnings)) {
        res$warnings$dblcent[[i]] <- parallel_results[[i]]$dblcent_warnings
      }
    }
  }
  
  res$rng_states <- rng_states_for_condition
  
  # summary during sim
  cat("\nConvergence rates:")
  cat("\n- SAM:", convergence_count$sam, "/", N_REPLICATIONS,
      sprintf("(%.1f%%)", 100 * convergence_count$sam / N_REPLICATIONS))
  cat("\n- QML:", convergence_count$qml, "/", N_REPLICATIONS,
      sprintf("(%.1f%%)", 100 * convergence_count$qml / N_REPLICATIONS))
  cat("\n- DBLCENT:", convergence_count$dblcent, "/", N_REPLICATIONS,
      sprintf("(%.1f%%)", 100 * convergence_count$dblcent / N_REPLICATIONS))
  
  cat("\n\nObserved metrics (means):")
  cat("\n- R² (eta4, eta5, eta6):", round(colMeans(res$observed_r2, na.rm = TRUE), 3))
  cat("\n- Reliabilities:", round(colMeans(res$observed_rel, na.rm = TRUE), 3))
  
  cat("\n\nMean computation time (seconds):")
  cat("\n- SAM:", round(mean(res$timing$sam,     na.rm = TRUE), 2))
  cat("\n- QML:", round(mean(res$timing$qml,     na.rm = TRUE), 2))
  cat("\n- DBLCENT:", round(mean(res$timing$dblcent, na.rm = TRUE), 2))
  cat("\n")
  
  # true parameters based on model type
  all_results[[cond]] <- list(
    condition   = conditions[cond, ],
    results     = res,
    convergence = convergence_count,
    true_parameters = if (conditions$Model_Type[cond] == "null") {
      list(
        # main effects only (same as alternative)
        eta4_eta1 = 0.21, eta4_eta2 = 0.21, eta4_eta3 = 0.21,
        eta5_eta4 = 0.18, eta5_eta1 = 0.18, eta5_eta2 = 0.18, eta5_eta3 = 0.18,
        eta6_eta5 = 0.15, eta6_eta1 = 0.15, eta6_eta2 = 0.15, eta6_eta3 = 0.15,
        # interactions and quadratics set to 0 for null model
        eta4_eta1eta2 = 0, eta4_eta1eta1 = 0, 
        eta5_eta2eta4 = 0, eta5_eta3eta3 = 0,
        eta6_eta1eta5 = 0, eta6_eta2eta3 = 0
      )
    } else {
      list(
        # main effects
        eta4_eta1 = 0.21, eta4_eta2 = 0.21, eta4_eta3 = 0.21,
        eta5_eta4 = 0.18, eta5_eta1 = 0.18, eta5_eta2 = 0.18, eta5_eta3 = 0.18,
        eta6_eta5 = 0.15, eta6_eta1 = 0.15, eta6_eta2 = 0.15, eta6_eta3 = 0.15,
        # interactions and quadratics (non-zero for alternative)
        eta4_eta1eta2 = 0.13, eta4_eta1eta1 = 0.09, 
        eta5_eta2eta4 = 0.11, eta5_eta2eta2 = 0.09,
        eta6_eta3eta5 = 0.10, eta6_eta3eta3 = 0.09
      )
    }
  )
  
  if (cond %% 3 == 0 || cond == nrow(conditions)) {
    saveRDS(all_results, file = sprintf("checkpoint_%d.rds", cond))
  }
  
  gc()  
}

stopCluster(cl)

#saveRDS(all_results, file = "final_simulation_results.rds")

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

total_time <- difftime(Sys.time(), start_time, units = "hours")
cat(sprintf("\n\nSimulation completed in %.2f hours\n", total_time))