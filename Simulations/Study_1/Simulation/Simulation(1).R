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

N_REPLICATIONS <- 1
SAMPLE_SIZES <- c(100, 1000)
SEED_START <- 123

# analysis model (for fitting - no fixed values)
analysis.model <- "
eta1 =~ x1 + x2 + x3
eta2 =~ x4 + x5 + x6
eta3 =~ x7 + x8 + x9

eta3 ~ eta1 + eta2 + eta1:eta2 + eta1:eta1 + eta2:eta2
"

# distribution parameters
distributions <- list(
  normal = list(skewness = c(0, 0), excesskurtosis = c(0, 0), distr.exo = "normal.rIG"),
  nonnormal = list(skewness = c(2, 2), excesskurtosis = c(7, 7), distr.exo = "nonnormal.rIG"),
  uniform = list(skewness = c(0, 0), excesskurtosis = c(0, 0), distr.exo = "unif")
)

# MAIN SIMULATION

dir.create("sim_results", showWarnings = FALSE)
timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
results_dir <- paste0("sim_results/run_", timestamp)
dir.create(results_dir, showWarnings = FALSE)

# conditions
conditions <- expand.grid(
  N = SAMPLE_SIZES,
  Rel = c(0.4, 0.6, 0.8),
  Distribution = names(distributions),
  Model_Type = c("alternative", "null"),  # model type; null = linear
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

# parallel setup
n_cores <- detectCores() - 6
cl <- makeCluster(n_cores)
registerDoParallel(cl)
clusterExport(cl, c("GenerateData", "method_analytic", "method_dblcent", "method_sam", 
                    "analysis.model", "all_models", "distributions"))

# run simulation
all_results <- list()
start_time <- Sys.time()

for(cond in 1:nrow(conditions)) {
  cat("\n========================================")
  cat("\nCondition", cond, "of", nrow(conditions))
  cat("\n- Sample size:", conditions$N[cond])
  cat("\n- Target reliability:", conditions$Rel[cond])
  cat("\n- Model type:", conditions$Model_Type[cond])
  cat("\n- Actual distribution:", conditions$Distribution[cond])
  cat("\n- Using model:", conditions$model_name[cond])
  cat("\n========================================\n")
  
  # get the appropriate model (normal or null/linear)
  population_model <- all_models[[conditions$model_name[cond]]]
  
  if(is.null(population_model)) {
    cat("WARNING: Model", conditions$model_name[cond], "not found. Skipping...\n")
    next
  }
  
  # distribution parameters for data generation
  dist_params <- distributions[[conditions$Distribution[cond]]]
  
  # condition-specific variables to cluster
  clusterExport(cl, c("population_model", "dist_params", "conditions", "cond"), 
                envir = environment())
  
  # Initialize results storage with error and warning tracking
  res <- list(
    lms_tables = vector("list", N_REPLICATIONS),
    qml_tables = vector("list", N_REPLICATIONS),
    dblcent_tables = vector("list", N_REPLICATIONS),
    sam_tables = vector("list", N_REPLICATIONS),
    timing = data.frame(lms = numeric(N_REPLICATIONS), 
                        qml = numeric(N_REPLICATIONS),
                        dblcent = numeric(N_REPLICATIONS), 
                        sam = numeric(N_REPLICATIONS)),
    # error tracking - stores error messages as strings
    errors = data.frame(lms = character(N_REPLICATIONS),
                        qml = character(N_REPLICATIONS),
                        dblcent = character(N_REPLICATIONS),
                        sam = character(N_REPLICATIONS),
                        stringsAsFactors = FALSE),
    # warning tracking - stores lists of warning messages
    warnings = list(lms = vector("list", N_REPLICATIONS),
                    qml = vector("list", N_REPLICATIONS),
                    dblcent = vector("list", N_REPLICATIONS),
                    sam = vector("list", N_REPLICATIONS)),
    # observed R^2 and reliabilities
    observed_r2 = numeric(N_REPLICATIONS),
    observed_rel = matrix(NA, nrow = N_REPLICATIONS, ncol = 9)
  )
  
  # run replications in parallel
  parallel_results <- foreach(i = 1:N_REPLICATIONS,
                              .packages = c("lavaan", "modsem", "covsim", "copula"),
                              .errorhandling = "pass",
                              .options.RNG = SEED_START + cond * 1000) %dorng% {
                                
                                # generate data using normal-based model with actual distribution
                                Data <- try(GenerateData(
                                  model = population_model,
                                  N = conditions$N[cond],
                                  skewness = dist_params$skewness,
                                  excesskurtosis = dist_params$excesskurtosis,
                                  distr.exo = dist_params$distr.exo,
                                  distr.epsilon = "normal",
                                  distr.zeta = "normal",
                                  add.eta = FALSE,
                                  return.info = TRUE
                                ), silent = TRUE)
                                
                                if(inherits(Data, "try-error")) return(NULL)
                                
                                # observed metrics
                                observed_metrics <- list(
                                  r2 = attr(Data, "observed_R2")$eta3,
                                  rel = unlist(attr(Data, "observed_reliabilities"))
                                )
                                
                                # remove attributes for methods
                                data_clean <- as.data.frame(Data)
                                attributes(data_clean) <- attributes(data_clean)[c("names", "row.names", "class")]
                                
                                # run all methods
                                results <- list(observed_metrics = observed_metrics)
                                
                                #=============== LMS METHOD ===============#
                                t0 <- Sys.time()
                                lms_warnings <- NULL 
                                
                                # withCallingHandlers wraps the try() to catch warnings
                                lms_table <- withCallingHandlers(
                                  # try() catches errors and returns them as objects instead of stopping
                                  try(method_analytic(Data = data_clean, 
                                                      model.fit = analysis.model, 
                                                      method = "lms"), silent = TRUE),
                                  # this function executes when a warning occurs
                                  warning = function(w) {
                                    # <<- assigns to parent environment (lms_warnings outside this function)
                                    # conditionMessage(w) extracts the warning text
                                    # c() concatenates to existing warnings (multiple warnings possible)
                                    lms_warnings <<- c(lms_warnings, conditionMessage(w))
                                    invokeRestart("muffleWarning")
                                  }
                                )
                                
                                # check if try() caught an error
                                if(!inherits(lms_table, "try-error")) {
                                  # no error: store successful results
                                  results$lms_table <- lms_table
                                  results$lms_timing <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
                                } else {
                                  # error occurred: store error message
                                  results$lms_error <- as.character(lms_table)
                                }
                                # store any warnings that were caught
                                results$lms_warnings <- lms_warnings
                                
                                #=============== QML METHOD ===============#
                                t0 <- Sys.time()
                                qml_warnings <- NULL
                                
                                qml_table <- withCallingHandlers(
                                  try(method_analytic(Data = data_clean, 
                                                      model.fit = analysis.model, 
                                                      method = "qml"), silent = TRUE),
                                  warning = function(w) {
                                    # same pattern: capture warning text and add to list
                                    qml_warnings <<- c(qml_warnings, conditionMessage(w))
                                    invokeRestart("muffleWarning")
                                  }
                                )
                                
                                if(!inherits(qml_table, "try-error")) {
                                  results$qml_table <- qml_table
                                  results$qml_timing <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
                                } else {
                                  results$qml_error <- as.character(qml_table)
                                }
                                results$qml_warnings <- qml_warnings
                                
                                #=============== UCA/DBLCENT METHOD ===============#
                                t0 <- Sys.time()
                                dblcent_warnings <- NULL  
                                
                                uca_table <- withCallingHandlers(
                                  try(method_dblcent(Data = data_clean, 
                                                     model.fit = analysis.model), silent = TRUE),
                                  warning = function(w) {
                                    dblcent_warnings <<- c(dblcent_warnings, conditionMessage(w))
                                    invokeRestart("muffleWarning")
                                  }
                                )
                                
                                if(!inherits(uca_table, "try-error")) {
                                  results$dblcent_table <- uca_table
                                  results$dblcent_timing <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
                                } else {
                                  results$dblcent_error <- as.character(uca_table)
                                }
                                results$dblcent_warnings <- dblcent_warnings
                                
                                #=============== SAM METHOD ===============#
                                t0 <- Sys.time()
                                sam_warnings <- NULL  
                                
                                sam_table <- withCallingHandlers(
                                  try(method_sam(Data = data_clean, 
                                                 model.fit = analysis.model), silent = TRUE),
                                  warning = function(w) {
                                    sam_warnings <<- c(sam_warnings, conditionMessage(w))
                                    invokeRestart("muffleWarning")
                                  }
                                )
                                
                                if(!inherits(sam_table, "try-error")) {
                                  results$sam_table <- sam_table
                                  results$sam_timing <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
                                } else {
                                  results$sam_error <- as.character(sam_table)
                                }
                                results$sam_warnings <- sam_warnings
                                
                                results
                              }
  
  # RNG states 
  rng_states_for_condition <- attr(parallel_results, "rng")
  
  # store results from parallel execution
  for(i in 1:N_REPLICATIONS) {
    if(!is.null(parallel_results[[i]]) && !inherits(parallel_results[[i]], "error")) {
      
      # observed metrics
      if(!is.null(parallel_results[[i]]$observed_metrics)) {
        res$observed_r2[i] <- parallel_results[[i]]$observed_metrics$r2
        res$observed_rel[i,] <- parallel_results[[i]]$observed_metrics$rel[1:9]
      }
      
      # LMS results, errors, and warnings
      if(!is.null(parallel_results[[i]]$lms_table)) {
        res$lms_tables[[i]] <- parallel_results[[i]]$lms_table
        res$timing$lms[i] <- parallel_results[[i]]$lms_timing
      }
      if(!is.null(parallel_results[[i]]$lms_error)) {
        res$errors$lms[i] <- parallel_results[[i]]$lms_error
      }
      if(!is.null(parallel_results[[i]]$lms_warnings)) {
        res$warnings$lms[[i]] <- parallel_results[[i]]$lms_warnings
      }
      
      # QML results, errors, and warnings
      if(!is.null(parallel_results[[i]]$qml_table)) {
        res$qml_tables[[i]] <- parallel_results[[i]]$qml_table
        res$timing$qml[i] <- parallel_results[[i]]$qml_timing
      }
      if(!is.null(parallel_results[[i]]$qml_error)) {
        res$errors$qml[i] <- parallel_results[[i]]$qml_error
      }
      if(!is.null(parallel_results[[i]]$qml_warnings)) {
        res$warnings$qml[[i]] <- parallel_results[[i]]$qml_warnings
      }
      
      # UCA results, errors, and warnings
      if(!is.null(parallel_results[[i]]$dblcent_table)) {
        res$dblcent_tables[[i]] <- parallel_results[[i]]$dblcent_table
        res$timing$dblcent[i] <- parallel_results[[i]]$dblcent_timing
      }
      if(!is.null(parallel_results[[i]]$dblcent_error)) {
        res$errors$dblcent[i] <- parallel_results[[i]]$dblcent_error
      }
      if(!is.null(parallel_results[[i]]$dblcent_warnings)) {
        res$warnings$dblcent[[i]] <- parallel_results[[i]]$dblcent_warnings
      }
      
      # SAM results, errors, and warnings
      if(!is.null(parallel_results[[i]]$sam_table)) {
        res$sam_tables[[i]] <- parallel_results[[i]]$sam_table
        res$timing$sam[i] <- parallel_results[[i]]$sam_timing
      }
      if(!is.null(parallel_results[[i]]$sam_error)) {
        res$errors$sam[i] <- parallel_results[[i]]$sam_error
      }
      if(!is.null(parallel_results[[i]]$sam_warnings)) {
        res$warnings$sam[[i]] <- parallel_results[[i]]$sam_warnings
      }
    }
  }
  
  res$rng_states <- rng_states_for_condition
  
  # summary including error/warning counts
  cat("\nObserved metrics across replications:")
  cat("\n- Mean R²:", mean(res$observed_r2, na.rm = TRUE))
  cat("\n- Mean reliabilities:", round(colMeans(res$observed_rel, na.rm = TRUE), 3))
  
  # error and warning summary
  cat("\n\nErrors encountered:")
  cat("\n- LMS:", sum(res$errors$lms != ""), "errors")
  cat("\n- QML:", sum(res$errors$qml != ""), "errors")
  cat("\n- DBLCENT:", sum(res$errors$dblcent != ""), "errors")
  cat("\n- SAM:", sum(res$errors$sam != ""), "errors")
  
  cat("\n\nWarnings encountered:")
  cat("\n- LMS:", sum(lengths(res$warnings$lms) > 0), "iterations with warnings")
  cat("\n- QML:", sum(lengths(res$warnings$qml) > 0), "iterations with warnings")
  cat("\n- DBLCENT:", sum(lengths(res$warnings$dblcent) > 0), "iterations with warnings")
  cat("\n- SAM:", sum(lengths(res$warnings$sam) > 0), "iterations with warnings")
  cat("\n")
  
  all_results[[cond]] <- list(
    condition = conditions[cond, ], 
    results = res,
    true_parameters = if(conditions$Model_Type[cond] == "null") {
      c(0.316, 0.316, 0, 0, 0)  # null (linear) model: interaction and quadratic terms are 0
    } else {
      c(0.316, 0.316, 0.139, 0.101, 0.101)  # full model with all effects
    }
  )
  
  # checkpoint
  if(cond %% 5 == 0 || cond == nrow(conditions)) {
    save(all_results, conditions, file = sprintf("%s/checkpoint_%d.RData", results_dir, cond))
  }
  
  gc()
}

stopCluster(cl)

#save(all_results, conditions, file = paste0(results_dir, "/final_results.RData"))

total_time <- difftime(Sys.time(), start_time, units = "hours")
cat(sprintf("\n\nSimulation completed in %.2f hours\n", total_time))