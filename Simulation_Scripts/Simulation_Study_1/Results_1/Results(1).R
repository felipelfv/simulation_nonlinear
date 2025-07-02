# ============================================================================
# ANALYSIS FUNCTIONS
# ============================================================================

library(simhelpers)

# Note: the function uses complete cases only for now

# Note: for clarity, when calculating the relative bias, when the true_value = 0, 
# we divided by zero. That is why we get multiple NAs when the pop model is linear.

# Type I error rate and power:
# As in Brandt et al. (2014), the type I error and power for the different 
# conditions is calculated. We follow their approach where "the percentage of
# significant t-values for the parameter estimates was calculated for each 
# linear and nonlinear effect of the structural model. These percentages are 
# interpreted as power for detecting the corresponding parameter when the 
# population parameter is different from zero, or as type I error when the 
# population parameter is equal to zero" (p. 186).

# Bias and squared bias:
# Following Rosseel et al. (2025), it is the difference between the estimate and
# the true population value (bias). If squared, then we have the squared bias. 

#### 2. Calculations ####

# Note: 
# When the true parameter value (gamma) is 0, the percentage bias formula 
# (Mdn(gamma_hat)) - gamma)/gamma · 100% would involve division by zero
# That is why the code returns NA for these cases.

# Hence, for parameters with true values of 0, 
# we use the absolute bias (median_bias) (?). To be discussed.

# Integrated function for bias, SE/SD ratio, CI coverage, power, Type I error, and MCSE calculations
CalculateBiasAndCoverage <- function(all_results, confidence_level = 0.95, 
                                     remove_outliers = FALSE, outlier_threshold = 5,
                                     alpha = 0.05, calculate_mcse = TRUE) {
  
  options(scipen = 999)
  
  if (is.null(all_results) || length(all_results) == 0) {
    stop("all_results is NULL or empty. Please load your data first.")
  }
  
  # true beta values 
  true_betas <- list(
    "population.linear.model" = c(0.316, 0.316, 0, 0, 0),
    "population.full.model" = c(0.316, 0.316, 0.139, 0.101, 0.101)
  )
  
  # Define the methods to process
  methods <- c("sam", "lms", "qml", "dblcent")
  
  results_summary <- data.frame()
  
  # each condition
  for (i in seq_along(all_results)) {
    condition <- all_results[[i]]$condition
    
    # process each method
    for (method in methods) {
      method_array <- all_results[[i]]$results[[method]]
      if (is.null(method_array) || !is.array(method_array)) next
      
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
        # get all estimates, standard errors, p-values, and CI bounds for this parameter
        all_estimates <- method_array[, p, 1]  
        all_std_errors <- method_array[, p, 2]
        all_pvalues <- method_array[, p, 3] 
        all_ci_lower <- method_array[, p, 4]
        all_ci_upper <- method_array[, p, 5]
        
        # ============================================================================
        # DATA QUALITY TRACKING - PART 1: NA COUNTS
        # ============================================================================
        n_total <- length(all_estimates)
        
        # Count NAs for each component
        n_na_estimates <- sum(is.na(all_estimates))
        n_na_std_errors <- sum(is.na(all_std_errors))
        n_na_pvalues <- sum(is.na(all_pvalues))
        n_na_ci_lower <- sum(is.na(all_ci_lower))
        n_na_ci_upper <- sum(is.na(all_ci_upper))
        
        # Count cases where ANY component has NA
        n_na_any_component <- sum(is.na(all_estimates) | is.na(all_std_errors) | 
                                    is.na(all_pvalues) | is.na(all_ci_lower) | 
                                    is.na(all_ci_upper))
        
        # Number of complete cases (no NAs in any component)
        n_complete_cases <- n_total - n_na_any_component
        
        # ============================================================================
        # REMOVE NAs - Keep only complete cases for analysis
        # ============================================================================
        valid_idx <- !is.na(all_estimates) & !is.na(all_std_errors) & 
          !is.na(all_pvalues) & !is.na(all_ci_lower) & !is.na(all_ci_upper)
        
        estimates <- all_estimates[valid_idx]
        std_errors <- all_std_errors[valid_idx]
        pvalues <- all_pvalues[valid_idx]
        ci_lower <- all_ci_lower[valid_idx]
        ci_upper <- all_ci_upper[valid_idx]
        
        if (length(estimates) < 10) next # skip if too few valid estimates
        
        # ============================================================================
        # DATA QUALITY TRACKING - PART 2: OUTLIER COUNTS
        # ============================================================================
        n_before_outlier_removal <- length(estimates)
        n_outliers_removed <- 0
        
        # remove outliers if requested
        if (remove_outliers && length(estimates) > 0) {
          q1 <- quantile(estimates, 0.25, na.rm = TRUE)
          q3 <- quantile(estimates, 0.75, na.rm = TRUE)
          iqr <- q3 - q1
          lower_bound <- q1 - outlier_threshold * iqr
          upper_bound <- q3 + outlier_threshold * iqr
          is_outlier <- estimates < lower_bound | estimates > upper_bound
          n_outliers_removed <- sum(is_outlier, na.rm = TRUE)
          
          # remove outliers from all arrays
          estimates <- estimates[!is_outlier]
          std_errors <- std_errors[!is_outlier]
          pvalues <- pvalues[!is_outlier]
          ci_lower <- ci_lower[!is_outlier]
          ci_upper <- ci_upper[!is_outlier]
        }
        
        if (length(estimates) < 10) next # skip if too few estimates remain
        
        # ============================================================================
        # DATA QUALITY TRACKING - PART 3: FINAL COUNTS
        # ============================================================================
        n_final_valid <- length(estimates)
        n_total_excluded <- n_total - n_final_valid
        
        # ============================================================================
        # CALCULATIONS
        # ============================================================================
        
        # (i). Mean-based statistics
        mean_estimate <- mean(estimates)
        mean_bias <- mean_estimate - true_values[p]
        squared_bias <- mean_bias^2
        
        if (true_values[p] == 0) {
          relative_bias_mean <- NA
          percent_relative_bias_mean <- NA
        } else {
          relative_bias_mean <- mean_bias / true_values[p]  
          percent_relative_bias_mean <- relative_bias_mean * 100     
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
          relative_bias_median <- NA
          percent_relative_bias_median <- NA 
        } else {
          relative_bias_median <- median_bias / true_values[p]
          percent_relative_bias_median <- relative_bias_median * 100 
        }
        
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
        
        # ============================================================================
        # MCSE CALCULATIONS (if requested)
        # ============================================================================
        if (calculate_mcse) {
          # Calculate MCSEs for key statistics
          # Note: Each MCSE corresponds to the statistic with the same name
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
            p_values = p_val, alpha = alpha
          )
          
          mean_bias_mcse <- abs_mcse$bias_mcse
          variance_mcse <- abs_mcse$var_mcse
          sd_estimates_mcse <- abs_mcse$stddev_mcse
          mse_mcse <- abs_mcse$mse_mcse
          coverage_rate_mcse <- cov_mcse$coverage_mcse
          rejection_rate_mcse <- rej_mcse$rej_rate_mcse
          
          # Assign MCSE to either Power or Type I Error
          if (true_values[p] == 0) {
            type_i_error_mcse <- rejection_rate_mcse
            power_mcse <- NA
          } else {
            type_i_error_mcse <- NA
            power_mcse <- rejection_rate_mcse
          }
        } else {
          # Set MCSEs to NA if not calculated
          mean_bias_mcse <- NA
          variance_mcse <- NA
          sd_estimates_mcse <- NA
          mse_mcse <- NA
          coverage_rate_mcse <- NA
          rejection_rate_mcse <- NA
          type_i_error_mcse <- NA
          power_mcse <- NA
        }
        
        results_summary <- rbind(results_summary, data.frame(
          Condition = i,
          Method = toupper(method),  
          Model = model_name,
          Distribution = distribution,
          SampleSize = sample_size,
          Reliability = reliability,
          Parameter = paste0("beta", p),
          TrueValue = true_values[p],
          
          # Mean-based statistics
          MeanEstimate = mean_estimate,
          MeanBias = mean_bias,
          RelativeBiasMean = relative_bias_mean,       
          PercentRelativeBiasMean = percent_relative_bias_mean,   
          Variance = variance_estimates,
          MSE = mse,
          
          # Median-based (robust) statistics
          MedianEstimate = median_estimate,
          MedianBias = median_bias,
          RelativeBiasMedian = relative_bias_median,
          PercentRelativeBiasMedian = percent_relative_bias_median,
          RMSE_Robust = rmse_robust,
          
          # Standard Error statistics
          MeanSE = mean_se,
          SD_Estimates = sd_estimates,
          SE_SD_Ratio = se_sd_ratio,
          
          # CI Coverage 
          CoverageRate = coverage_rate,
          
          # Power and Type I Error
          RejectionRate = rejection_rate,
          TypeI_Error = type_i_error,
          Power = power,
          
          # MCSE values
          MeanBias_MCSE = mean_bias_mcse,
          Variance_MCSE = variance_mcse,
          SD_Estimates_MCSE = sd_estimates_mcse,
          MSE_MCSE = mse_mcse,
          CoverageRate_MCSE = coverage_rate_mcse,
          RejectionRate_MCSE = rejection_rate_mcse,
          TypeI_Error_MCSE = type_i_error_mcse,
          Power_MCSE = power_mcse,
          
          # Data quality - NA counts
          N_Total = n_total,
          N_NA_Estimates = n_na_estimates,
          N_NA_StdErrors = n_na_std_errors,
          N_NA_PValues = n_na_pvalues,
          N_NA_CI_Lower = n_na_ci_lower,
          N_NA_CI_Upper = n_na_ci_upper,
          N_NA_Any = n_na_any_component,
          N_Complete_Cases = n_complete_cases,
          
          # Data quality - Outlier counts
          N_Outliers_Removed = n_outliers_removed,
          Percent_Outliers = ifelse(n_before_outlier_removal > 0, 
                                    round(n_outliers_removed / n_before_outlier_removal * 100, 2), 
                                    0),
          
          # Data quality - Final counts
          N_Final_Valid = n_final_valid,
          N_Total_Excluded = n_total_excluded,
          Percent_Total_Excluded = round(n_total_excluded / n_total * 100, 2),
          
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  results_summary 
}

# Usage:
# With MCSE calculations (default)
results_with_mcse <- CalculateBiasAndCoverage(all_results)

# Without MCSE calculations (faster if not needed)
results_no_mcse <- CalculateBiasAndCoverage(all_results, calculate_mcse = FALSE)

