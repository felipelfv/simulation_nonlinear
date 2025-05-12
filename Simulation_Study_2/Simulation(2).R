#### 1. General Information ####

# This code is not processing multiple conditions in parallel. 
# It processes the conditions one after another. 
# The parallelization happens only within each condition, 
# where multiple replications are run simultaneously across X cores.

# Other scripts needed
source("Simulation/GenerateData.R")
source("Simulation/Methods.R")
source("Simulation/Models.R")
source("Simulation/Design.R")


# Run comparison with all 4 methods
run_comparison <- function(N = 500, seed = 123) {
  # Generate data
  data <- GenerateData(
    model = population.large.model,
    N = N,
    skewness = skewness,
    excesskurtosis = excesskurtosis,
    exo.mean = exo.mean,
    center.exogenous.latent = TRUE,
    center.exogenous.manifest = FALSE,
    center.lv.dependent = FALSE,
    center.lv.prod = FALSE,
    center.indicators = FALSE,
    distr.exo = "rIG", # Remember to pick skewness and kurtosis accordingly
    distr.zeta = "normal",
    distr.epsilon = "normal",
    rel = Rel,
    #target.var = target.var,
    R2 = R2,
    add.eta = FALSE)
  
  # Apply all methods with timing
  results <- list()
  timing <- list()
  
  # UCA method
  timing$UCA <- system.time({
    results$UCA <- method_uca(Data = data, model.fit = fit.large.model)
  })
  
  # SAM method
  timing$SAM <- system.time({
    results$SAM <- method_sam(Data = data, model.fit = fit.large.model)
  })
  
  # LMS method
  timing$LMS <- system.time({
    results$LMS <- method_analytic(Data = data, model.fit = fit.large.model, method = "lms")
  })
  
  # QML method
  timing$QML <- system.time({
    results$QML <- method_analytic(Data = data, model.fit = fit.large.model, method = "qml")
  })
  
  # Create comparison table
  params <- c("eta1", "eta2", "eta3", "eta2:eta3", "eta1:eta3", "eta6:eta2", 
              "eta5:eta6", "eta5:eta5", "eta1:eta2")
  
  # Extract estimates
  estimates <- data.frame(Parameter = params)
  for (method in names(results)) {
    estimates[[method]] <- results[[method]]$Estimates
  }
  
  # Extract standard errors (not available for SAM)
  se_methods <- c("UCA", "LMS", "QML")
  se <- data.frame(Parameter = params)
  for (method in se_methods) {
    se[[method]] <- results[[method]]$`Standard Errors`
  }
  
  # Extract p-values (not available for SAM)
  pvalues <- data.frame(Parameter = params)
  for (method in se_methods) {
    pvalues[[method]] <- results[[method]]$`P-values`
  }
  
  # Create timing table
  timing_df <- data.frame(
    Method = names(timing),
    User_Time = sapply(timing, function(x) x["user.self"]),
    System_Time = sapply(timing, function(x) x["sys.self"]),
    Elapsed_Time = sapply(timing, function(x) x["elapsed"])
  )
  
  return(list(
    estimates = estimates,
    standard_errors = se,
    p_values = pvalues,
    timing = timing_df,
    detailed_results = results
  ))
}


results <- run_comparison(N = 500)

(results$estimates)
(results$standard_errors)
(results$p_values)
(results$timing)

# Save results to file if needed
# saveRDS(results, "method_comparison_results.rds")
