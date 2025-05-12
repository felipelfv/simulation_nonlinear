library(dplyr)
library(ggplot2)

# Load data function
loadResultData <- function() {
  # Load the main data file
  load("run_20250407_1017/final_results.RData")
  
  # Create all result variables in an environment
  results_env <- new.env()
  
  # (a) Subset for latent distributions varying, epsilon normal, 200, 0.4
  results_env$results_12_lms <- lapply(combined_results[1:12], function(x) x$results$lms[, , "beta"])
  results_env$results_12_qml <- lapply(combined_results[1:12], function(x) x$results$qml[, , "beta"])
  results_env$results_12_uca <- lapply(combined_results[1:12], function(x) x$results$uca[, , "beta"])
  results_env$results_12_sam <- lapply(combined_results[1:12], function(x) x$results$sam)
  
  # (b) Subset for latent distributions varying, epsilon exp.rate1, 200, 0.4
  results_env$results_24_lms <- lapply(combined_results[13:24], function(x) x$results$lms[, , "beta"])
  results_env$results_24_qml <- lapply(combined_results[13:24], function(x) x$results$qml[, , "beta"])
  results_env$results_24_uca <- lapply(combined_results[13:24], function(x) x$results$uca[, , "beta"])
  results_env$results_24_sam <- lapply(combined_results[13:24], function(x) x$results$sam)
  
  # (c) Subset for latent distributions varying, epsilon normal, 500, 0.4
  results_env$results_36_lms <- lapply(combined_results[25:36], function(x) x$results$lms[, , "beta"])
  results_env$results_36_qml <- lapply(combined_results[25:36], function(x) x$results$qml[, , "beta"])
  results_env$results_36_uca <- lapply(combined_results[25:36], function(x) x$results$uca[, , "beta"])
  results_env$results_36_sam <- lapply(combined_results[25:36], function(x) x$results$sam)
  
  # (d) Subset for latent distributions varying, epsilon exp.rate1, 500, 0.4
  results_env$results_48_lms <- lapply(combined_results[37:48], function(x) x$results$lms[, , "beta"])
  results_env$results_48_qml <- lapply(combined_results[37:48], function(x) x$results$qml[, , "beta"])
  results_env$results_48_uca <- lapply(combined_results[37:48], function(x) x$results$uca[, , "beta"])
  results_env$results_48_sam <- lapply(combined_results[37:48], function(x) x$results$sam)
  
  # (e) Subset for latent distributions varying, epsilon normal, 800, 0.4
  results_env$results_60_lms <- lapply(combined_results[49:60], function(x) x$results$lms[, , "beta"])
  results_env$results_60_qml <- lapply(combined_results[49:60], function(x) x$results$qml[, , "beta"])
  results_env$results_60_uca <- lapply(combined_results[49:60], function(x) x$results$uca[, , "beta"])
  results_env$results_60_sam <- lapply(combined_results[49:60], function(x) x$results$sam)
  
  # (f) Subset for latent distributions varying, epsilon exp.rate1, 800, 0.4
  results_env$results_72_lms <- lapply(combined_results[61:72], function(x) x$results$lms[, , "beta"])
  results_env$results_72_qml <- lapply(combined_results[61:72], function(x) x$results$qml[, , "beta"])
  results_env$results_72_uca <- lapply(combined_results[61:72], function(x) x$results$uca[, , "beta"])
  results_env$results_72_sam <- lapply(combined_results[61:72], function(x) x$results$sam)
  
  # (g) Subset for latent distributions varying, epsilon normal, 200, 0.6
  results_env$results_84_lms <- lapply(combined_results[73:84], function(x) x$results$lms[, , "beta"])
  results_env$results_84_qml <- lapply(combined_results[73:84], function(x) x$results$qml[, , "beta"])
  results_env$results_84_uca <- lapply(combined_results[73:84], function(x) x$results$uca[, , "beta"])
  results_env$results_84_sam <- lapply(combined_results[73:84], function(x) x$results$sam)
  
  # (h) Subset for latent distributions varying, epsilon exp.rate1, 200, 0.6
  results_env$results_96_lms <- lapply(combined_results[85:96], function(x) x$results$lms[, , "beta"])
  results_env$results_96_qml <- lapply(combined_results[85:96], function(x) x$results$qml[, , "beta"])
  results_env$results_96_uca <- lapply(combined_results[85:96], function(x) x$results$uca[, , "beta"])
  results_env$results_96_sam <- lapply(combined_results[85:96], function(x) x$results$sam)
  
  # (i) Subset for latent distributions varying, epsilon normal, 500, 0.6
  results_env$results_108_lms <- lapply(combined_results[97:108], function(x) x$results$lms[, , "beta"])
  results_env$results_108_qml <- lapply(combined_results[97:108], function(x) x$results$qml[, , "beta"])
  results_env$results_108_uca <- lapply(combined_results[97:108], function(x) x$results$uca[, , "beta"])
  results_env$results_108_sam <- lapply(combined_results[97:108], function(x) x$results$sam)
  
  # (j) Subset for latent distributions varying, epsilon exp.rate1, 500, 0.6
  results_env$results_120_lms <- lapply(combined_results[109:120], function(x) x$results$lms[, , "beta"])
  results_env$results_120_qml <- lapply(combined_results[109:120], function(x) x$results$qml[, , "beta"])
  results_env$results_120_uca <- lapply(combined_results[109:120], function(x) x$results$uca[, , "beta"])
  results_env$results_120_sam <- lapply(combined_results[109:120], function(x) x$results$sam)
  
  # (k) Subset for latent distributions varying, epsilon normal, 800, 0.6
  results_env$results_132_lms <- lapply(combined_results[121:132], function(x) x$results$lms[, , "beta"])
  results_env$results_132_qml <- lapply(combined_results[121:132], function(x) x$results$qml[, , "beta"])
  results_env$results_132_uca <- lapply(combined_results[121:132], function(x) x$results$uca[, , "beta"])
  results_env$results_132_sam <- lapply(combined_results[121:132], function(x) x$results$sam)
  
  # (l) Subset for latent distributions varying, epsilon exp.rate1, 800, 0.6
  results_env$results_144_lms <- lapply(combined_results[133:144], function(x) x$results$lms[, , "beta"])
  results_env$results_144_qml <- lapply(combined_results[133:144], function(x) x$results$qml[, , "beta"])
  results_env$results_144_uca <- lapply(combined_results[133:144], function(x) x$results$uca[, , "beta"])
  results_env$results_144_sam <- lapply(combined_results[133:144], function(x) x$results$sam)
  
  # (m) Subset for latent distributions varying, epsilon normal, 200, 0.8
  results_env$results_156_lms <- lapply(combined_results[145:156], function(x) x$results$lms[, , "beta"])
  results_env$results_156_qml <- lapply(combined_results[145:156], function(x) x$results$qml[, , "beta"])
  results_env$results_156_uca <- lapply(combined_results[145:156], function(x) x$results$uca[, , "beta"])
  results_env$results_156_sam <- lapply(combined_results[145:156], function(x) x$results$sam)
  
  # (n) Subset for latent distributions varying, epsilon exp.rate1, 200, 0.8
  results_env$results_168_lms <- lapply(combined_results[157:168], function(x) x$results$lms[, , "beta"])
  results_env$results_168_qml <- lapply(combined_results[157:168], function(x) x$results$qml[, , "beta"])
  results_env$results_168_uca <- lapply(combined_results[157:168], function(x) x$results$uca[, , "beta"])
  results_env$results_168_sam <- lapply(combined_results[157:168], function(x) x$results$sam)
  
  # (o) Subset for latent distributions varying, epsilon normal, 500, 0.8
  results_env$results_180_lms <- lapply(combined_results[169:180], function(x) x$results$lms[, , "beta"])
  results_env$results_180_qml <- lapply(combined_results[169:180], function(x) x$results$qml[, , "beta"])
  results_env$results_180_uca <- lapply(combined_results[169:180], function(x) x$results$uca[, , "beta"])
  results_env$results_180_sam <- lapply(combined_results[169:180], function(x) x$results$sam)
  
  # (p) Subset for latent distributions varying, epsilon exp.rate1, 500, 0.8
  results_env$results_192_lms <- lapply(combined_results[181:192], function(x) x$results$lms[, , "beta"])
  results_env$results_192_qml <- lapply(combined_results[181:192], function(x) x$results$qml[, , "beta"])
  results_env$results_192_uca <- lapply(combined_results[181:192], function(x) x$results$uca[, , "beta"])
  results_env$results_192_sam <- lapply(combined_results[181:192], function(x) x$results$sam)
  
  # (q) Subset for latent distributions varying, epsilon normal, 800, 0.8
  results_env$results_204_lms <- lapply(combined_results[193:204], function(x) x$results$lms[, , "beta"])
  results_env$results_204_qml <- lapply(combined_results[193:204], function(x) x$results$qml[, , "beta"])
  results_env$results_204_uca <- lapply(combined_results[193:204], function(x) x$results$uca[, , "beta"])
  results_env$results_204_sam <- lapply(combined_results[193:204], function(x) x$results$sam)
  
  # (r) Subset for latent distributions varying, epsilon exp.rate1, 800, 0.8
  results_env$results_216_lms <- lapply(combined_results[205:216], function(x) x$results$lms[, , "beta"])
  results_env$results_216_qml <- lapply(combined_results[205:216], function(x) x$results$qml[, , "beta"])
  results_env$results_216_uca <- lapply(combined_results[205:216], function(x) x$results$uca[, , "beta"])
  results_env$results_216_sam <- lapply(combined_results[205:216], function(x) x$results$sam)
  
  results_env
}

# analyze and visualize results based on input parameters
analyzeResults <- function(comparison_type, condition_num, parameter_name, remove_outliers = TRUE) {

  results_data <- loadResultData()
  # all results available in the function environment
  list2env(as.list(results_data), envir = environment())
  
  # true beta values for each model
  true_betas <- list(
    "population.linear.model" = c(0.316, 0.316, 0, 0, 0),
    "population.interaction.model" = c(0.316, 0.316, 0.139, 0, 0),
    "population.full.model" = c(0.316, 0.316, 0.139, 0.101, 0.101)
  )
  
  # Map from population model names
  population_model_map <- c(
    "Linear-Interaction" = "population.linear.model",
    "Linear-Full" = "population.linear.model",
    "Interaction-Interaction" = "population.interaction.model",
    "Full-Full" = "population.full.model"
  )
  
  # Map condition number to letter for clarity
  condition_letter <- letters[as.numeric(condition_num)]
  
  # Step 1: Set up the data for the specified condition
  if (comparison_type == "sample_size") {
    if (condition_num == 1) {
      # Condition 1a: Sample size comparison, reliability 0.4, normal residuals
      result_sets <- list(
        list(lms = results_12_lms, qml = results_12_qml, uca = results_12_uca, sam = results_12_sam),
        list(lms = results_36_lms, qml = results_36_qml, uca = results_36_uca, sam = results_36_sam),
        list(lms = results_60_lms, qml = results_60_qml, uca = results_60_uca, sam = results_60_sam)
      )
      group_values <- c(200, 500, 800)
      title <- "Sample Size Effect: Normal Residuals, 0.4 Reliability"
      factor_name <- "Sample Size"
    } else if (condition_num == 2) {
      # Condition 1b: Sample size comparison, reliability 0.4, exponential residuals
      result_sets <- list(
        list(lms = results_24_lms, qml = results_24_qml, uca = results_24_uca, sam = results_24_sam),
        list(lms = results_48_lms, qml = results_48_qml, uca = results_48_uca, sam = results_48_sam),
        list(lms = results_72_lms, qml = results_72_qml, uca = results_72_uca, sam = results_72_sam)
      )
      group_values <- c(200, 500, 800)
      title <- "Sample Size Effect: Exponential Residuals, 0.4 Reliability"
      factor_name <- "Sample Size"
    } else if (condition_num == 3) {
      # Condition 1c: Sample size comparison, reliability 0.6, normal residuals
      result_sets <- list(
        list(lms = results_84_lms, qml = results_84_qml, uca = results_84_uca, sam = results_84_sam),
        list(lms = results_108_lms, qml = results_108_qml, uca = results_108_uca, sam = results_108_sam),
        list(lms = results_132_lms, qml = results_132_qml, uca = results_132_uca, sam = results_132_sam)
      )
      group_values <- c(200, 500, 800)
      title <- "Sample Size Effect: Normal Residuals, 0.6 Reliability"
      factor_name <- "Sample Size"
    } else if (condition_num == 4) {
      # Condition 1d: Sample size comparison, reliability 0.6, exponential residuals
      result_sets <- list(
        list(lms = results_96_lms, qml = results_96_qml, uca = results_96_uca, sam = results_96_sam),
        list(lms = results_120_lms, qml = results_120_qml, uca = results_120_uca, sam = results_120_sam),
        list(lms = results_144_lms, qml = results_144_qml, uca = results_144_uca, sam = results_144_sam)
      )
      group_values <- c(200, 500, 800)
      title <- "Sample Size Effect: Exponential Residuals, 0.6 Reliability"
      factor_name <- "Sample Size"
    } else if (condition_num == 5) {
      # Condition 1e: Sample size comparison, reliability 0.8, normal residuals
      result_sets <- list(
        list(lms = results_156_lms, qml = results_156_qml, uca = results_156_uca, sam = results_156_sam),
        list(lms = results_180_lms, qml = results_180_qml, uca = results_180_uca, sam = results_180_sam),
        list(lms = results_204_lms, qml = results_204_qml, uca = results_204_uca, sam = results_204_sam)
      )
      group_values <- c(200, 500, 800)
      title <- "Sample Size Effect: Normal Residuals, 0.8 Reliability"
      factor_name <- "Sample Size"
    } else if (condition_num == 6) {
      # Condition 1f: Sample size comparison, reliability 0.8, exponential residuals
      result_sets <- list(
        list(lms = results_168_lms, qml = results_168_qml, uca = results_168_uca, sam = results_168_sam),
        list(lms = results_192_lms, qml = results_192_qml, uca = results_192_uca, sam = results_192_sam),
        list(lms = results_216_lms, qml = results_216_qml, uca = results_216_uca, sam = results_216_sam)
      )
      group_values <- c(200, 500, 800)
      title <- "Sample Size Effect: Exponential Residuals, 0.8 Reliability"
      factor_name <- "Sample Size"
    } else {
      stop("Invalid condition number (must be 1-6)")
    }
  } else if (comparison_type == "reliability") {
    if (condition_num == 1) {
      # Condition 2a: Reliability comparison, sample size 200, normal residuals
      result_sets <- list(
        list(lms = results_12_lms, qml = results_12_qml, uca = results_12_uca, sam = results_12_sam),
        list(lms = results_84_lms, qml = results_84_qml, uca = results_84_uca, sam = results_84_sam),
        list(lms = results_156_lms, qml = results_156_qml, uca = results_156_uca, sam = results_156_sam)
      )
      group_values <- c(0.4, 0.6, 0.8)
      title <- "Reliability Effect: Normal Residuals, N=200"
      factor_name <- "Reliability"
    } else if (condition_num == 2) {
      # Condition 2b: Reliability comparison, sample size 200, exponential residuals
      result_sets <- list(
        list(lms = results_24_lms, qml = results_24_qml, uca = results_24_uca, sam = results_24_sam),
        list(lms = results_96_lms, qml = results_96_qml, uca = results_96_uca, sam = results_96_sam),
        list(lms = results_168_lms, qml = results_168_qml, uca = results_168_uca, sam = results_168_sam)
      )
      group_values <- c(0.4, 0.6, 0.8)
      title <- "Reliability Effect: Exponential Residuals, N=200"
      factor_name <- "Reliability"
    } else if (condition_num == 3) {
      # Condition 2c: Reliability comparison, sample size 500, normal residuals
      result_sets <- list(
        list(lms = results_36_lms, qml = results_36_qml, uca = results_36_uca, sam = results_36_sam),
        list(lms = results_108_lms, qml = results_108_qml, uca = results_108_uca, sam = results_108_sam),
        list(lms = results_180_lms, qml = results_180_qml, uca = results_180_uca, sam = results_180_sam)
      )
      group_values <- c(0.4, 0.6, 0.8)
      title <- "Reliability Effect: Normal Residuals, N=500"
      factor_name <- "Reliability"
    } else if (condition_num == 4) {
      # Condition 2d: Reliability comparison, sample size 500, exponential residuals
      result_sets <- list(
        list(lms = results_48_lms, qml = results_48_qml, uca = results_48_uca, sam = results_48_sam),
        list(lms = results_120_lms, qml = results_120_qml, uca = results_120_uca, sam = results_120_sam),
        list(lms = results_192_lms, qml = results_192_qml, uca = results_192_uca, sam = results_192_sam)
      )
      group_values <- c(0.4, 0.6, 0.8)
      title <- "Reliability Effect: Exponential Residuals, N=500"
      factor_name <- "Reliability"
    } else if (condition_num == 5) {
      # Condition 2e: Reliability comparison, sample size 800, normal residuals
      result_sets <- list(
        list(lms = results_60_lms, qml = results_60_qml, uca = results_60_uca, sam = results_60_sam),
        list(lms = results_132_lms, qml = results_132_qml, uca = results_132_uca, sam = results_132_sam),
        list(lms = results_204_lms, qml = results_204_qml, uca = results_204_uca, sam = results_204_sam)
      )
      group_values <- c(0.4, 0.6, 0.8)
      title <- "Reliability Effect: Normal Residuals, N=800"
      factor_name <- "Reliability"
    } else if (condition_num == 6) {
      # Condition 2f: Reliability comparison, sample size 800, exponential residuals
      result_sets <- list(
        list(lms = results_72_lms, qml = results_72_qml, uca = results_72_uca, sam = results_72_sam),
        list(lms = results_144_lms, qml = results_144_qml, uca = results_144_uca, sam = results_144_sam),
        list(lms = results_216_lms, qml = results_216_qml, uca = results_216_uca, sam = results_216_sam)
      )
      group_values <- c(0.4, 0.6, 0.8)
      title <- "Reliability Effect: Exponential Residuals, N=800"
      factor_name <- "Reliability"
    } else {
      stop("Invalid condition number (must be 1-6)")
    }
  } else if (comparison_type == "residual_distribution") {
    # Create residual comparisons configurations
    residual_comparisons <- list(
      # N=200, Reliability=0.4
      "1" = list(
        result_sets = list(
          list(lms = results_12_lms, qml = results_12_qml, uca = results_12_uca, sam = results_12_sam),
          list(lms = results_24_lms, qml = results_24_qml, uca = results_24_uca, sam = results_24_sam)
        ),
        group_values = c("Normal", "Exponential"),
        title = "Residual Distribution Effect: N=200, 0.4 Reliability",
        factor_name = "Residual Distribution"
      ),
      # N=200, Reliability=0.6
      "2" = list(
        result_sets = list(
          list(lms = results_84_lms, qml = results_84_qml, uca = results_84_uca, sam = results_84_sam),
          list(lms = results_96_lms, qml = results_96_qml, uca = results_96_uca, sam = results_96_sam)
        ),
        group_values = c("Normal", "Exponential"),
        title = "Residual Distribution Effect: N=200, 0.6 Reliability",
        factor_name = "Residual Distribution"
      ),
      # N=200, Reliability=0.8
      "3" = list(
        result_sets = list(
          list(lms = results_156_lms, qml = results_156_qml, uca = results_156_uca, sam = results_156_sam),
          list(lms = results_168_lms, qml = results_168_qml, uca = results_168_uca, sam = results_168_sam)
        ),
        group_values = c("Normal", "Exponential"),
        title = "Residual Distribution Effect: N=200, 0.8 Reliability",
        factor_name = "Residual Distribution"
      ),
      # N=500, Reliability=0.4
      "4" = list(
        result_sets = list(
          list(lms = results_36_lms, qml = results_36_qml, uca = results_36_uca, sam = results_36_sam),
          list(lms = results_48_lms, qml = results_48_qml, uca = results_48_uca, sam = results_48_sam)
        ),
        group_values = c("Normal", "Exponential"),
        title = "Residual Distribution Effect: N=500, 0.4 Reliability",
        factor_name = "Residual Distribution"
      ),
      # N=500, Reliability=0.6
      "5" = list(
        result_sets = list(
          list(lms = results_108_lms, qml = results_108_qml, uca = results_108_uca, sam = results_108_sam),
          list(lms = results_120_lms, qml = results_120_qml, uca = results_120_uca, sam = results_120_sam)
        ),
        group_values = c("Normal", "Exponential"),
        title = "Residual Distribution Effect: N=500, 0.6 Reliability",
        factor_name = "Residual Distribution"
      ),
      # N=500, Reliability=0.8
      "6" = list(
        result_sets = list(
          list(lms = results_180_lms, qml = results_180_qml, uca = results_180_uca, sam = results_180_sam),
          list(lms = results_192_lms, qml = results_192_qml, uca = results_192_uca, sam = results_192_sam)
        ),
        group_values = c("Normal", "Exponential"),
        title = "Residual Distribution Effect: N=500, 0.8 Reliability",
        factor_name = "Residual Distribution"
      ),
      # N=800, Reliability=0.4
      "7" = list(
        result_sets = list(
          list(lms = results_60_lms, qml = results_60_qml, uca = results_60_uca, sam = results_60_sam),
          list(lms = results_72_lms, qml = results_72_qml, uca = results_72_uca, sam = results_72_sam)
        ),
        group_values = c("Normal", "Exponential"),
        title = "Residual Distribution Effect: N=800, 0.4 Reliability",
        factor_name = "Residual Distribution"
      ),
      # N=800, Reliability=0.6
      "8" = list(
        result_sets = list(
          list(lms = results_132_lms, qml = results_132_qml, uca = results_132_uca, sam = results_132_sam),
          list(lms = results_144_lms, qml = results_144_qml, uca = results_144_uca, sam = results_144_sam)
        ),
        group_values = c("Normal", "Exponential"),
        title = "Residual Distribution Effect: N=800, 0.6 Reliability",
        factor_name = "Residual Distribution"
      ),
      # N=800, Reliability=0.8
      "9" = list(
        result_sets = list(
          list(lms = results_204_lms, qml = results_204_qml, uca = results_204_uca, sam = results_204_sam),
          list(lms = results_216_lms, qml = results_216_qml, uca = results_216_uca, sam = results_216_sam)
        ),
        group_values = c("Normal", "Exponential"),
        title = "Residual Distribution Effect: N=800, 0.8 Reliability",
        factor_name = "Residual Distribution"
      )
    )
    
    if (condition_num < 1 || condition_num > 9) {
      stop("Invalid condition number for residual distribution comparison (must be 1-9)")
    }
    
    # select comparison data
    comparison_data <- residual_comparisons[[as.character(condition_num)]]
    
    # components 
    result_sets <- comparison_data$result_sets
    group_values <- comparison_data$group_values
    title <- comparison_data$title
    factor_name <- comparison_data$factor_name
  }
  else {
    stop("Invalid comparison type. Use 'sample_size', 'reliability', or 'residual_distribution'")
  }
  
  # parameter index
  param_index <- as.numeric(gsub("beta", "", parameter_name))
  
  # Step 2: Prepare the combined data from all result sets
  all_data <- data.frame()
  
  # Process each result set
  for (i in seq_along(result_sets)) {
    # Extract the results
    lms_results <- result_sets[[i]]$lms
    qml_results <- result_sets[[i]]$qml 
    uca_results <- result_sets[[i]]$uca
    sam_results <- result_sets[[i]]$sam
    
    # Process each condition
    for (j in seq_along(lms_results)) {
      # Determine model and distribution from indices
      model_index <- (j - 1) %% 4 + 1
      dist_index <- ceiling(j / 4)
      
      # Define model and distribution names
      model_names <- c("Linear-Interaction", "Linear-Full", "Interaction-Interaction", "Full-Full")
      dist_names <- c("Normal-rIG", "Nonnormal-Unif", "Nonnormal-rIG")
      
      model <- model_names[model_index]
      distribution <- dist_names[dist_index]
      
      # Skip if any result is missing
      if (is.null(dim(lms_results[[j]])) || is.null(dim(qml_results[[j]])) || 
          is.null(dim(uca_results[[j]])) || is.null(dim(sam_results[[j]]))) {
        next
      }
      
      # Process LMS method
      if (!is.null(dim(lms_results[[j]]))) {
        n_params <- ncol(lms_results[[j]])
        for (p in 1:n_params) {
          param <- paste0("beta", p)
          if (param == parameter_name) {
            all_data <- rbind(all_data, data.frame(
              Beta = lms_results[[j]][, p],
              Parameter = param,
              Method = "LMS",
              Model = model,
              Distribution = distribution,
              Group = group_values[i]
            ))
          }
        }
      }
      
      # Process QML method
      if (!is.null(dim(qml_results[[j]]))) {
        n_params <- ncol(qml_results[[j]])
        for (p in 1:n_params) {
          param <- paste0("beta", p)
          if (param == parameter_name) {
            all_data <- rbind(all_data, data.frame(
              Beta = qml_results[[j]][, p],
              Parameter = param,
              Method = "QML",
              Model = model,
              Distribution = distribution,
              Group = group_values[i]
            ))
          }
        }
      }
      
      # Process UCA method
      if (!is.null(dim(uca_results[[j]]))) {
        n_params <- ncol(uca_results[[j]])
        for (p in 1:n_params) {
          param <- paste0("beta", p)
          if (param == parameter_name) {
            all_data <- rbind(all_data, data.frame(
              Beta = uca_results[[j]][, p],
              Parameter = param,
              Method = "UCA",
              Model = model,
              Distribution = distribution,
              Group = group_values[i]
            ))
          }
        }
      }
      
      # Process SAM method
      if (!is.null(dim(sam_results[[j]]))) {
        n_params <- ncol(sam_results[[j]])
        for (p in 1:n_params) {
          param <- paste0("beta", p)
          if (param == parameter_name) {
            all_data <- rbind(all_data, data.frame(
              Beta = sam_results[[j]][, p],
              Parameter = param,
              Method = "SAM",
              Model = model,
              Distribution = distribution,
              Group = group_values[i]
            ))
          }
        }
      }
    }
  }
  
  # check if we found any data for the parameter
  if (nrow(all_data) == 0) {
    stop(paste("No data found for parameter", parameter_name))
  }
  
  # Step 3: Process outliers
  
  # outlier thresholds - now grouped by Model, Distribution, Method, AND Group
  outlier_info <- all_data %>%
    group_by(Model, Distribution, Method, Group) %>%
    summarize(
      q1 = quantile(Beta, 0.25, na.rm = TRUE),
      q3 = quantile(Beta, 0.75, na.rm = TRUE),
      iqr = q3 - q1,
      lower_bound = q1 - 5 * iqr,  # 5 IQRs for extreme outliers
      upper_bound = q3 + 5 * iqr,
      .groups = "drop"
    )
  
  # flg outliers - now joined on all four grouping variables
  all_data <- all_data %>%
    left_join(outlier_info, by = c("Model", "Distribution", "Method", "Group")) %>%
    mutate(
      is_extreme_outlier = Beta < lower_bound | Beta > upper_bound
    )
  
  # Count outliers
  outlier_count <- sum(all_data$is_extreme_outlier, na.rm = TRUE)
  
  # Apply outlier handling
  if (remove_outliers) {
    all_data <- all_data %>%
      mutate(Beta_plot = ifelse(is_extreme_outlier, NA, Beta))
    outlier_note <- ifelse(outlier_count > 0,
                           paste0("\nNote: ", outlier_count, " extreme outliers removed"),
                           "")
  } else {
    all_data <- all_data %>%
      mutate(Beta_plot = Beta)
    outlier_note <- "\nNote: All outliers included"
  }
  
  # Step 4: Prepare summary statistics
  
  # Overall NA counts by method
  na_summary <- all_data %>%
    group_by(Method) %>%
    summarize(
      Total = n(),
      NAs = sum(is.na(Beta)),
      NA_percent = round(sum(is.na(Beta)) / n() * 100, 2)
    )
  
  # Overall Outlier counts by method
  outlier_summary <- all_data %>%
    filter(!is.na(Beta)) %>%
    group_by(Method) %>%
    summarize(
      Total_non_NA = sum(!is.na(Beta)),
      Outliers = sum(is_extreme_outlier, na.rm = TRUE),
      Outlier_percent = round(sum(is_extreme_outlier, na.rm = TRUE) / sum(!is.na(Beta)) * 100, 2)
    )
  
  # Detailed NA counts by the full grouping
  na_detailed <- all_data %>%
    group_by(Method, Model, Distribution, Group) %>%
    summarize(
      Total = n(),
      NAs = sum(is.na(Beta)),
      NA_percent = round(sum(is.na(Beta)) / n() * 100, 2),
      .groups = "drop"
    ) %>%
    arrange(Method, Model, Distribution, Group)
  
  # Detailed Outlier counts by the full grouping
  outlier_detailed <- all_data %>%
    filter(!is.na(Beta)) %>%
    group_by(Method, Model, Distribution, Group) %>%
    summarize(
      Total_non_NA = sum(!is.na(Beta)),
      Outliers = sum(is_extreme_outlier, na.rm = TRUE),
      Outlier_percent = round(sum(is_extreme_outlier, na.rm = TRUE) / sum(!is.na(Beta)) * 100, 2),
      .groups = "drop"
    ) %>%
    arrange(Method, Model, Distribution, Group)
  
  # Step 5: Prepare true values for plotting
  
  # Create a data frame to store true values for each model
  true_values <- data.frame()
  
  # Get unique model and group combinations
  unique_models <- unique(all_data$Model)
  unique_groups <- unique(all_data$Group)
  unique_distributions <- unique(all_data$Distribution)
  
  # For each model, add the true parameter value
  for (model_name in unique_models) {
    # Get the corresponding population model
    pop_model <- population_model_map[model_name]
    
    # Get the true value for this parameter
    if (param_index <= length(true_betas[[pop_model]])) {
      true_value <- true_betas[[pop_model]][param_index]
      
      # Add a row for each group and distribution combination
      for (group_val in unique_groups) {
        for (dist_val in unique_distributions) {
          true_values <- rbind(true_values, data.frame(
            Model = model_name,
            Distribution = dist_val,
            Group = group_val,
            TrueValue = true_value
          ))
        }
      }
    }
  }
  
  # Filter data for plotting (remove NAs after outlier handling)
  all_data_plot <- all_data %>%
    filter(!is.na(Beta_plot))
  
  # list of results
  list(
    all_data = all_data,
    all_data_plot = all_data_plot,
    true_values = true_values,
    title = title,
    factor_name = factor_name,
    group_values = group_values,
    outlier_note = outlier_note,
    na_summary = na_summary,
    outlier_summary = outlier_summary,
    na_detailed = na_detailed,
    outlier_detailed = outlier_detailed,
    parameter_name = parameter_name,
    comparison_type = comparison_type,
    condition_num = condition_num
  )
}

# Function to create and save the parameter boxplot
createParameterPlot <- function(results, filename = NULL, width = 10, height = 7) {
  all_data_plot <- results$all_data_plot
  true_values <- results$true_values
  
  if (results$comparison_type == "residual_distribution") {
    all_data_plot$Group <- factor(all_data_plot$Group, levels = results$group_values)
  } else {
    all_data_plot$Group <- factor(all_data_plot$Group)
    true_values$Group <- factor(true_values$Group)
  }
  
  plot <- ggplot() +
    # Add boxplots
    geom_boxplot(data = all_data_plot, 
                 aes(x = Group, y = Beta_plot, fill = Method),
                 alpha = 0.7, position = position_dodge(width = 0.85), outlier.size = 1) +
    # Add true value lines
    geom_hline(data = true_values %>% 
                 select(Model, Distribution, TrueValue) %>%
                 distinct(), # Keep only unique Model/Distribution/TrueValue combinations
               aes(yintercept = TrueValue), 
               color = "red", linetype = "solid", linewidth = 0.5) +
    # Facet and themeing
    facet_grid(. ~ Model + Distribution, scales = "free_y") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 14, face = "bold"),
      strip.text = element_text(size = 8),
      legend.position = "bottom"
    ) +
    labs(
      title = paste0(results$title, " - ", results$parameter_name),
      subtitle = paste0("Comparison factor: ", results$factor_name, results$outlier_note),
      x = results$factor_name,
      y = paste0(results$parameter_name, " Coefficient Value")
    ) +
    scale_fill_brewer(palette = "Set1")
  
  # Save plot if filename is provided
  if (!is.null(filename)) {
    ggsave(filename, plot = plot, width = width, height = height, dpi = 300)
  }
  
  plot
}

# Function to create and save the missing values (NA) plot
createNAPlot <- function(results, filename = NULL, width = 10, height = 7) {
  # Create data summary for NA plot
  na_summary <- results$all_data %>%
    group_by(Method, Model, Distribution, Group) %>%
    summarize(
      Total = n(),
      NAs = sum(is.na(Beta)),
      NA_percent = round(sum(is.na(Beta)) / n() * 100, 2),
      .groups = "drop"
    )
  
  # Create facet labels for the plot
  na_summary$facet_label <- paste(na_summary$Model, na_summary$Distribution, sep = "\n")
  
  # Convert Group to factor for proper ordering
  na_summary$Group <- factor(na_summary$Group)
  
  plot <- ggplot(na_summary, aes(x = Group, y = NA_percent, fill = Method)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    facet_wrap(~ facet_label, ncol = 4) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 14, face = "bold"),
      strip.text = element_text(size = 8),
      legend.position = "bottom"
    ) +
    labs(
      title = paste0("Missing Values ", results$parameter_name),
      subtitle = paste0("Comparison factor: ", results$factor_name),
      x = results$factor_name,
      y = "Percentage of Missing Values (%)"
    ) +
    scale_fill_brewer(palette = "Set1") +
    geom_text(aes(label = NA_percent), position = position_dodge(width = 0.9), 
              vjust = -0.5, size = 2.5)
  
  # Save plot if filename is provided
  if (!is.null(filename)) {
    ggsave(filename, plot = plot, width = width, height = height, dpi = 300)
  }
  
  plot
}

# Function to create and save the outliers plot
createOutlierPlot <- function(results, filename = NULL, width = 10, height = 7) {
  # Create data summary for outlier plot
  outlier_summary <- results$all_data %>%
    filter(!is.na(Beta)) %>%
    group_by(Method, Model, Distribution, Group) %>%
    summarize(
      Total_non_NA = sum(!is.na(Beta)),
      Outliers = sum(is_extreme_outlier, na.rm = TRUE),
      Outlier_percent = round(sum(is_extreme_outlier, na.rm = TRUE) / sum(!is.na(Beta)) * 100, 2),
      .groups = "drop"
    )
  
  # Create facet labels for the plot
  outlier_summary$facet_label <- paste(outlier_summary$Model, outlier_summary$Distribution, sep = "\n")
  
  # Convert Group to factor for proper ordering
  outlier_summary$Group <- factor(outlier_summary$Group)
  
  plot <- ggplot(outlier_summary, aes(x = Group, y = Outlier_percent, fill = Method)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    facet_wrap(~ facet_label, ncol = 4) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 14, face = "bold"),
      strip.text = element_text(size = 8),
      legend.position = "bottom"
    ) +
    labs(
      title = paste0("Extreme Outliers ", results$parameter_name),
      subtitle = paste0("Comparison factor: ", results$factor_name, 
                        "\nNote: Outliers defined as values >5 IQR from Q1/Q3"),
      x = results$factor_name,
      y = "Percentage of Outliers (%)"
    ) +
    scale_fill_brewer(palette = "Set1") +
    geom_text(aes(label = Outlier_percent), position = position_dodge(width = 0.9), 
              vjust = -0.5, size = 2.5)
  
  # Save plot if filename is provided
  if (!is.null(filename)) {
    ggsave(filename, plot = plot, width = width, height = height, dpi = 300)
  }
  
  plot
}

load("run_20250407_1017/final_results.RData")
# combined_results <- all_results

# sample size effect with condition 1 (Normal Residuals, 0.4 Reliability) for beta1:
results <- analyzeResults(comparison_type = "sample_size", condition_num = 1, parameter_name = "beta1", remove_outliers = TRUE)

parameter_plot <- createParameterPlot(results)
na_plot <- createNAPlot(results)
outlier_plot <- createOutlierPlot(results)

# summary
# (results$na_summary)
# (results$outlier_summary)
# (results$na_detailed)
# (results$outlier_detailed)