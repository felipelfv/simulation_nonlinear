library(simhelpers)

# Integrated function for bias, SE/SD ratio, CI coverage, power, and Type I error calculations
CalculateBiasAndCoverage <- function(all_results, confidence_level = 0.95, 
                                     remove_outliers = FALSE, outlier_threshold = 5,
                                     alpha = 0.05) {
  
  options(scipen = 999)
  
  if (is.null(all_results) || length(all_results) == 0) {
    stop("all_results is NULL or empty. Please load your data first.")
  }
  
  # true beta values 
  true_betas <- list(
    "population.linear.model" = c(0.316, 0.316, 0, 0, 0),
    "population.full.model" = c(0.316, 0.316, 0.139, 0.101, 0.101)
  )
  
  results_summary <- data.frame()
  
  # each condition
  for (i in seq_along(all_results)) {
    condition <- all_results[[i]]$condition
    sam_array <- all_results[[i]]$results$sam
    if (is.null(sam_array) || !is.array(sam_array)) next
    
    # model type and get true values
    pop_model <- as.character(condition$Population)
    true_values <- true_betas[[pop_model]]
    
    # model and distribution names for plotting
    model_name <- ifelse(pop_model == "population.linear.model", "Linear", "Full")
    
    # map latent distribution
    latent_dist <- as.character(condition$Distribution)
    latent_type <- as.character(condition$Exo_method)
    if (latent_dist == "normal") {
      distribution <- "Normal-rIG"
    } else if (latent_dist == "nonnormal" && latent_type == "rIG") {
      distribution <- "Nonnormal-rIG"
    } else if (latent_dist == "nonnormal" && latent_type == "unif") {
      distribution <- "Nonnormal-Unif"
    }
    
    # sample size and reliability
    sample_size <- condition$N
    reliability <- condition$Rel
    
    # each parameter
    for (p in 1:5) {
      # estimates, standard errors, p-values, and CI bounds for each parameter
      all_estimates <- sam_array[, p, 1]  
      all_std_errors <- sam_array[, p, 2]
      all_pvalues <- sam_array[, p, 3] 
      all_ci_lower <- sam_array[, p, 4]
      all_ci_upper <- sam_array[, p, 5]
      
      # count NAs
      n_total <- length(all_estimates)
      n_na_estimates <- sum(is.na(all_estimates))
      n_na_std_errors <- sum(is.na(all_std_errors))
      n_na_pvalues <- sum(is.na(all_pvalues))
      n_na_ci <- sum(is.na(all_ci_lower) | is.na(all_ci_upper))
      n_na_either <- sum(is.na(all_estimates) | is.na(all_std_errors) | 
                           is.na(all_pvalues) | is.na(all_ci_lower) | is.na(all_ci_upper))
      
      # remove NA values for calculations - need all components to be valid
      valid_idx <- !is.na(all_estimates) & !is.na(all_std_errors) & 
        !is.na(all_pvalues) & !is.na(all_ci_lower) & !is.na(all_ci_upper)
      estimates <- all_estimates[valid_idx]
      std_errors <- all_std_errors[valid_idx]
      pvalues <- all_pvalues[valid_idx]
      ci_lower <- all_ci_lower[valid_idx]
      ci_upper <- all_ci_upper[valid_idx]
      
      if (length(estimates) < 10) next # skip if too few valid estimates
      
      n_original <- length(estimates)
      n_outliers <- 0
      
      # remove outliers if requested
      if (remove_outliers && length(estimates) > 0) {
        q1 <- quantile(estimates, 0.25, na.rm = TRUE)
        q3 <- quantile(estimates, 0.75, na.rm = TRUE)
        iqr <- q3 - q1
        lower_bound <- q1 - outlier_threshold * iqr
        upper_bound <- q3 + outlier_threshold * iqr
        is_outlier <- estimates < lower_bound | estimates > upper_bound
        n_outliers <- sum(is_outlier, na.rm = TRUE)
        
        # Remove outliers from all arrays
        estimates <- estimates[!is_outlier]
        std_errors <- std_errors[!is_outlier]
        pvalues <- pvalues[!is_outlier]
        ci_lower <- ci_lower[!is_outlier]
        ci_upper <- ci_upper[!is_outlier]
      }
      
      if (length(estimates) < 10) next # skip if too few estimates remain
      
      # (i). Mean-based statistics
      mean_estimate <- mean(estimates)
      mean_bias <- mean_estimate - true_values[p]
      squared_bias <- mean_bias^2
      
      if (true_values[p] == 0) {
        relative_bias <- NA
        percent_bias_mean <- NA
      } else {
        relative_bias <- mean_bias / true_values[p]  # relative bias
        percent_bias_mean <- relative_bias * 100     # relative bias * 100
      }
      
      # (ii). Variance, SD, and SE-related statistics
      sd_estimates <- sd(estimates) # empirical standard error
      variance_estimates <- sd_estimates^2
      mse <- squared_bias + variance_estimates 
      mean_se <- mean(std_errors) # average of the SEs within each simulation
      se_sd_ratio <- mean_se / sd_estimates
      
      # (iii). Median-based (robust) statistics
      median_estimate <- median(estimates, na.rm = TRUE)
      mad_estimate <- mad(estimates, na.rm = TRUE) 
      
      median_bias <- median_estimate - true_values[p] 
      
      if (true_values[p] == 0) {
        percent_bias_median <- NA 
      } else {
        percent_bias_median <- (median_bias / true_values[p]) * 100 
      }
      
      abs_percent_bias_median <- abs(percent_bias_median) 
      is_small_bias <- ifelse(!is.na(abs_percent_bias_median), abs_percent_bias_median < 10, NA) 
      rmse_robust <- sqrt(median_bias^2 + mad_estimate^2)
      
      # (iv). CI Coverage calculations
      true_value <- true_values[p]
      covered <- (ci_lower <= true_value) & (ci_upper >= true_value)
      coverage_rate <- mean(covered) * 100
      
      # (v). Power and Type I Error Rate calculations
      # Test H0: beta = 0
      significant <- pvalues < alpha
      rejection_rate <- mean(significant) * 100
      
      # Determine if this is a Type I error or Power calculation
      if (true_values[p] == 0) {
        # True value is 0, so H0 is true - this is Type I error rate
        type_i_error <- rejection_rate
        power <- NA
      } else {
        # True value is not 0, so H0 is false - this is Power
        type_i_error <- NA
        power <- rejection_rate
      }
      
      # Additional p-value statistics
      #mean_pvalue <- mean(pvalues)
      #median_pvalue <- median(pvalues)
      
      results_summary <- rbind(results_summary, data.frame(
        Condition = i,
        Model = model_name,
        Distribution = distribution,
        SampleSize = sample_size,
        Reliability = reliability,
        Parameter = paste0("beta", p),
        TrueValue = true_values[p],
        
        # Mean-based 
        MeanEstimate = mean_estimate,
        MeanBias = mean_bias,
        RelativeBias = relative_bias,       
        PercentBiasMean = percent_bias_mean,   # this is relative bias * 100
        Variance = variance_estimates,
        MSE = mse,
        
        # Median-based (robust)
        MedianEstimate = median_estimate,
        MedianBias = median_bias,
        PercentBiasMedian = percent_bias_median,
        AbsPercentBiasMedian = abs_percent_bias_median,
        RMSE_Robust = rmse_robust,
        
        # Standard Error 
        MeanSE = mean_se,
        SD_Estimates = sd_estimates,
        SE_SD_Ratio = se_sd_ratio,
        
        # CI Coverage 
        CoverageRate = coverage_rate,
        
        # Power and Type I Error
        RejectionRate = rejection_rate,
        TypeI_Error = type_i_error,
        Power = power,
        #MeanPValue = mean_pvalue,
        #MedianPValue = median_pvalue,
        
        # Data quality 
        N_total = n_total,
        N_NA_either = n_na_either,
        N_valid = length(estimates),
        N_outliers = n_outliers,
        Outlier_percent = round(n_outliers / n_original * 100, 2),
        
        stringsAsFactors = FALSE
      ))
    }
  }
  
  results_summary 
}

# Run your analysis first
results_summary <- CalculateBiasAndCoverage(all_results)

# Add MCSE columns to results_summary
mcse_results <- data.frame()

# Loop through all conditions and parameters
for (i in seq_along(all_results)) {
  sam_array <- all_results[[i]]$results$sam
  if (is.null(sam_array) || !is.array(sam_array)) next
  
  # Get true values for this condition
  pop_model <- as.character(all_results[[i]]$condition$Population)
  true_betas <- list(
    "population.linear.model" = c(0.316, 0.316, 0, 0, 0),
    "population.full.model" = c(0.316, 0.316, 0.139, 0.101, 0.101)
  )
  
  for (p in 1:5) {
    # Extract data
    estimates <- sam_array[, p, 1]
    pvalues <- sam_array[, p, 3]
    ci_lower <- sam_array[, p, 4]
    ci_upper <- sam_array[, p, 5]
    
    # Remove NAs
    valid_idx <- !is.na(estimates) & !is.na(pvalues) & !is.na(ci_lower) & !is.na(ci_upper)
    estimates <- estimates[valid_idx]
    pvalues <- pvalues[valid_idx]
    ci_lower <- ci_lower[valid_idx]
    ci_upper <- ci_upper[valid_idx]
    
    if (length(estimates) < 10) next
    
    true_value <- true_betas[[pop_model]][p]
    
    # Calculate all MCSEs at once
    abs_mcse <- calc_absolute(
      data = data.frame(est = estimates, true_param = true_value),
      estimates = est, true_param = true_param,
      criteria = c("bias", "variance", "stddev", "mse")
    )
    
    cov_mcse <- calc_coverage(
      data = data.frame(lower_bound = ci_lower, upper_bound = ci_upper, true_param = true_value),
      lower_bound = lower_bound, upper_bound = upper_bound, true_param = true_param
    )
    
    rej_mcse <- calc_rejection(
      data = data.frame(p_val = pvalues),
      p_values = p_val, alpha = 0.05
    )
    
    mcse_results <- rbind(mcse_results, data.frame(
      Condition = i,
      Parameter = paste0("beta", p),
      Bias_MCSE = abs_mcse$bias_mcse,
      Variance_MCSE = abs_mcse$var_mcse,
      EmpSE_MCSE = abs_mcse$stddev_mcse,
      MSE_MCSE = abs_mcse$mse_mcse,
      Coverage_MCSE = cov_mcse$coverage_mcse,
      RejRate_MCSE = rej_mcse$rej_rate_mcse
    ))
  }
}

# Merge with results_summary
results_with_mcse <- merge(results_summary, mcse_results, 
                           by = c("Condition", "Parameter"))

# View first few rows
head(results_with_mcse)