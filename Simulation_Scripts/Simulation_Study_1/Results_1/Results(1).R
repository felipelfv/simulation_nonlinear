#### 1. General Session Information ####

# This script contains the code for calculating the performance measures

##### 1.1 Information for each measure #####

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

options(scipen=999) # Avoid scientific notation

process_beta_data <- function(all_results) {
  # True beta values for each population model
  true_betas <- list(
    # B1(eta1), B2(eta2), B3(eta1:eta2), B4(eta1:eta1), and B5(eta2:eta2)
    "population.linear.model" = c(0.316, 0.316, 0, 0, 0), 
    "population.interaction.model" = c(0.316, 0.316, 0.139, 0, 0),
    "population.full.model" = c(0.316, 0.316, 0.139, 0.101, 0.101)
  )
  
  beta_data <- data.frame()
  
  for (i in seq_along(all_results)) {
    data <- all_results[[i]]
    
    condition <- data$condition
    pop_model <- condition$Population
    true_beta <- true_betas[[pop_model]]
    methods <- c("lms", "qml", "uca", "sam")
    
    for (method in methods) {
      if (!is.null(data$results[[method]])) {
        beta_matrix <- data$results[[method]][,,"beta"]
        se_matrix <- data$results[[method]][,,"se"]
        
        na_counts <- colSums(is.na(beta_matrix)) # NAs for each beta coefficient
        # Outliers (values with absolute value > 1)
        outlier_matrix <- !is.na(beta_matrix) & abs(beta_matrix) > 1
        outlier_counts <- colSums(outlier_matrix)
        working_matrix <- beta_matrix # working copy of the matrix
        # Outliers to NA in the working matrix
        working_matrix[outlier_matrix] <- NA
        
        # Mean-based statistics
        mean_betas <- colMeans(working_matrix, na.rm = TRUE)
        variance <- apply(working_matrix, 2, var, na.rm = TRUE)
        
        # Calculate Monte-Carlo SD (standard deviation of parameter estimates)
        mc_sd <- sqrt(variance)
        
        # Calculate average SE for each coefficient
        # Remove outliers from SE matrix too (using same outlier criteria from beta matrix)
        se_working_matrix <- se_matrix
        se_working_matrix[outlier_matrix] <- NA
        avg_se <- colMeans(se_working_matrix, na.rm = TRUE)
        
        # Calculate SE/SD ratio
        se_sd_ratio <- avg_se / mc_sd
        
        # Median-based statistics (as described in Brandt et al., 2020) 
        median_betas <- apply(working_matrix, 2, median, na.rm = TRUE)
        
        # Calculate MAD for each beta (robust measure of variability)
        # in Brandt et al., 2020 (p. 332)
        mad_betas <- numeric(length(median_betas))
        for (j in seq_along(median_betas)) {
          abs_deviations <- abs(working_matrix[,j] - median_betas[j])
          mad_betas[j] <- median(abs_deviations, na.rm = TRUE)
        }
        
        # Here we have the data for each beta coefficient
        for (j in seq_along(mean_betas)) {
          true_val <- if(j <= length(true_beta)) true_beta[j] else 0
          
          # Bias as per formula 13 in Brandt et al., 2020 (p. 332)
          percent_bias_mdn <- ifelse(true_val != 0,
                                     ((median_betas[j] - true_val) / true_val) * 100,
                                     NA)  # avoid division by zero
          
          # Absolute percent bias 
          abs_percent_bias <- abs(percent_bias_mdn)
          
          # RMSE as per formula 14 in Brandt et al., 2020 (p. 332)
          rmse_mdn <- sqrt((median_betas[j] - true_val)^2 + mad_betas[j]^2)
          
          beta_data <- rbind(beta_data, data.frame(
            condition_id = i,
            population = pop_model,
            distribution = condition$Distribution,
            exo_method = condition$Exo_method,
            epsilon = condition$Epsilon,
            n = condition$N,
            reliability = condition$Rel,
            analysis_model = condition$Analysis_model,
            method = method,
            beta_index = j,
            beta_name = paste0("B", j),
            
            # Mean-based statistics
            mean_estimate = mean_betas[j],
            mean_bias = mean_betas[j] - true_val,
            squared_bias = (mean_betas[j] - true_val)^2,
            variance = variance[j], 
            mse = (mean_betas[j] - true_val)^2 + variance[j],
            
            # SE/SD ratio calculation
            avg_se = avg_se[j],
            mc_sd = mc_sd[j],
            se_sd_ratio = se_sd_ratio[j],
            
            # Median-based statistics (from Brandt et al., 2020)
            median_estimate = median_betas[j],
            median_bias = median_betas[j] - true_val,
            percent_bias_mdn = percent_bias_mdn,
            abs_percent_bias = abs_percent_bias,
            small_bias = ifelse(!is.na(abs_percent_bias), abs_percent_bias < 10, NA),
            mad = mad_betas[j],
            rmse_mdn = rmse_mdn,
            
            # Counts
            na_count = na_counts[j],
            outlier_count = outlier_counts[j]
          ))
        }
      }
    }
  }
  
  beta_data
}

results <- process_beta_data(all_results)
