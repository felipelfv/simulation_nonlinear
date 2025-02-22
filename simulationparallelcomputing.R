#library(parallel)
library(foreach)
library(doParallel)

# Other scripts needed
source("GenerateData_yr.R")
source("Method.analytic.R")
source("Method.uca.R")
source("Method.sam.R")
source("Models.R")

# Conditions as in Brand et al. (2014) - not fully factorial
conditions_interaction <- expand.grid(
  pop_model      = c("population.linear.model", "population.interaction.model"),
  analysis_model = "fit.interaction.model",
  distribution   = c("normal", "nonnormal"),
  stringsAsFactors = FALSE
)
conditions_full <- expand.grid(
  pop_model      = c("population.linear.model", "population.full.model"),
  analysis_model = "fit.full.model", 
  distribution   = c("normal", "nonnormal"),
  stringsAsFactors = FALSE
)
conditions <- rbind(conditions_interaction, conditions_full)
conditions$epsilon <- ifelse(conditions$distribution == "normal", "normal", "exp.rate1")


# Setup parallel backend
n_cores <- detectCores() - 2 # 8 in total, we use 6; change this when using the semlab pc
cl <- makeCluster(n_cores)
registerDoParallel(cl)

clusterExport(cl, c("GenerateData", "method_analytic", "method_uca", "method_sam", "population.interaction.model",
                    "population.linear.model", "population.full.model", "fit.interaction.model", "fit.full.model"))

# Basic parameters 
N = 400L # Later this should vary: 200L, 400L, and 800L 
exo.mean <- rep(0,2)
rel = 0.80 # Later this should vary: 0.2, 0.6, 0.8 
target.var <- list("eta3" = 1.0)
R2 <- list("eta3" = 0.20)
rep <- 25

start_time <- Sys.time()
all_results <- list()

for(cond in 1:nrow(conditions)) { #
  condition_start <- Sys.time()
  cat("\nRunning condition", cond, "of", nrow(conditions), "\n")
  
  if(conditions$distribution[cond] == "normal") {
    skewness <- rep(0, 2)
    excesskurtosis <- rep(0, 2)
  } else {
    skewness <- rep(2, 2)
    excesskurtosis <- rep(7, 2)
  }
  
  n_params <- if(conditions$analysis_model[cond] == "fit.full.model") 5 else 3
  
  # Initialize results structure
  res <- list(
    lms = array(NA, dim = c(rep, n_params, 3),
                dimnames = list(NULL, NULL, c("beta", "se", "pval"))),
    qml = array(NA, dim = c(rep, n_params, 3),
                dimnames = list(NULL, NULL, c("beta", "se", "pval"))),
    uca = array(NA, dim = c(rep, n_params, 3),
                dimnames = list(NULL, NULL, c("beta", "se", "pval"))),
    sam = matrix(NA, rep, n_params)
  )
  
  # Parallel processing of replications
  results <- foreach(i = seq_len(rep), 
                     .packages = c("modsem", "lavaan", "covsim"), 
                     .errorhandling = "pass") %dopar% {
                       set.seed(1234 + i)
                       
                       local_res <- list(
                         lms = array(NA, dim = c(1, n_params, 3)),
                         qml = array(NA, dim = c(1, n_params, 3)),
                         uca = array(NA, dim = c(1, n_params, 3)),
                         sam = rep(NA, n_params)
                       )
                       
                       Data <- try(GenerateData(
                         model = get(conditions$pop_model[cond]),
                         N = N,
                         skewness = skewness,
                         excesskurtosis = excesskurtosis,
                         exo.mean = exo.mean,
                         distr.zeta = "normal",
                         distr.epsilon = conditions$epsilon[cond],
                         rel = rel,
                         target.var = target.var,
                         R2 = R2,
                         add.eta = FALSE), silent = TRUE)
                       
                       if(!inherits(Data, "try-error")) {
                         analysis_model <- get(conditions$analysis_model[cond])
                         
                         # LMS and QML methods
                         for(m in c("lms", "qml")) {
                           result <- try(method_analytic(Data = Data, model.fit = analysis_model, method = m))
                           if(!inherits(result, "try-error")) {
                             local_res[[m]][1, , 1] <- result$Estimates
                             local_res[[m]][1, , 2] <- result$`Standard Errors`
                             local_res[[m]][1, , 3] <- result$`P-values`
                           }
                         }
                         
                         # UCA method
                         result <- try(method_uca(Data = Data, model.fit = analysis_model))
                         if(!inherits(result, "try-error")) {
                           local_res$uca[1, , 1] <- result$Estimates
                           local_res$uca[1, , 2] <- result$`Standard Errors`
                           local_res$uca[1, , 3] <- result$`P-values`
                         }
                         
                         # SAM method
                         result <- try(method_sam(Data = Data, model.fit = analysis_model))
                         if(!inherits(result, "try-error")) {
                           local_res$sam <- result$Estimates
                         }
                       }
                       
                       local_res
                     }
  
  # Combine parallel results
  # As of now, we have a list of 6(?) local_res for each replication 
  for(i in seq_len(rep)) {
    if(!inherits(results[[i]], "try-error")) {
      res$lms[i,,] <- results[[i]]$lms[1,,]
      res$qml[i,,] <- results[[i]]$qml[1,,]
      res$uca[i,,] <- results[[i]]$uca[1,,]
      res$sam[i,] <- results[[i]]$sam # dimensions as this for now because we dont use SE and p-values
    }
  }
  
  all_results[[cond]] <- list(
    condition = conditions[cond, ],
    results = res
  )
  
  save(all_results, file = paste0("simulation_results_condition_", cond, ".RData"))
  
  condition_time <- difftime(Sys.time(), condition_start, units = "mins")
  cat(sprintf("\nCondition %d completed in %.2f minutes\n", cond, condition_time))
}

# Check also `stopImplicitCluster()`
stopCluster(cl)
total_time <- difftime(Sys.time(), start_time, units = "hours")
cat(sprintf("\nTotal simulation completed in %.2f hours\n", total_time))

# We need to combine results from all conditions (!)
all_results <- list()
n_conditions <- nrow(conditions)

for(cond in 1:n_conditions) {
  filename <- paste0("simulation_results_condition_", cond, ".RData")
  load(filename)
  all_results[[cond]] <- get("all_results")[[cond]]
}

# Save combined results
save(all_results, file = "complete_simulation_results.RData")


# standardized solution to see the standard regression coefficients


