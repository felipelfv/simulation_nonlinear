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

options(scipen=999) # To avoid scientific notation

##### 2.1 Average Betas (Bias) #####
all_beta_averages_with_bias <- function(all_results) {
    avgs <- list()
    diff_bias <- list()
    squared_bias <- list()
  
  # True beta values based on model type
  true_betas <- list(
    "population.linear.model" = c(0.316, 0.316, 0, 0, 0),
    "population.interaction.model" = c(0.316, 0.316, 0.139, 0, 0),
    "population.full.model" = c(0.316, 0.316, 0.139, 0.101, 0.101)
  )
  
  # each element in all_results
  for(i in seq_along(all_results)) {
    data <- all_results[[i]]
    result <- list()
    bias_result <- list()
    squared_bias_results <- list()
    
    # True betas for this model type
    model_type <- data$condition$Population
    true_beta <- true_betas[[model_type]]
    
    # Process standard methods (lms, qml, uca)
    for(method in c("lms", "qml", "uca")) {
      if(!is.null(data$results[[method]])) {
        # Calculate average betas
        avg_beta <- colMeans(data$results[[method]][,,"beta"])
        result[[method]] <- avg_beta
        
        # Calculate bias (estimated - true)
        # Make sure vectors are the same length
        len <- min(length(avg_beta), length(true_beta))
        bias <- avg_beta[1:len] - true_beta[1:len]
        bias_result[[method]] <- bias
        squared_bias_results[[method]] <- bias^2
      }
    }
    
    # Handle SAM method
    if(!is.null(data$results$sam)) {
      avg_beta <- colMeans(data$results$sam)
      result$sam <- avg_beta
      
      # Calculate bias for SAM
      len <- min(length(avg_beta), length(true_beta))
      bias <- avg_beta[1:len] - true_beta[1:len]
      bias_result$sam <- bias
      squared_bias_results$sam <- bias^2
    }
    
    # Add condition information to bias
    bias_result$condition <- data$condition
    squared_bias_results$condition <- data$condition
    
    avgs[[i]] <- result
    diff_bias[[i]] <- bias_result
    squared_bias[[i]] <- squared_bias_results
    
  }
  
  list(averages = avgs, bias = diff_bias, squared_bias = squared_bias)
}

results <- all_beta_averages_with_bias(all_results)




