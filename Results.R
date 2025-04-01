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

options(scipen=999) # To avoid scientific notation, 999 high threshold

# with na and outlier
process_beta_data <- function(all_results) {
  # True beta values for each population model
  true_betas <- list(
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
        if (method == "sam") {
          beta_matrix <- data$results[[method]]
        } else {
          beta_matrix <- data$results[[method]][,,"beta"]
        }

        na_counts <- colSums(is.na(beta_matrix)) # NAs for each beta coefficient
        # Oiutliers (values with absolute value > 1), so above or below
        outlier_matrix <- !is.na(beta_matrix) & abs(beta_matrix) > 1
        outlier_counts <- colSums(outlier_matrix)
        working_matrix <- beta_matrix # working copy of the matrix
        # Outliers to NA in the working matrix
        working_matrix[outlier_matrix] <- NA
        mean_betas <- colMeans(working_matrix, na.rm = TRUE)
        
        # Here we have the data for each beta coefficient
        for (j in seq_along(mean_betas)) {
          true_val <- if(j <= length(true_beta)) true_beta[j] else 0
          
          variance <- var(working_matrix[,j], na.rm = TRUE)
          
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
            estimate = mean_betas[j],
            true_value = true_val,
            bias = mean_betas[j] - true_val,
            squared_bias = (mean_betas[j] - true_val)^2,
            variance = variance, 
            mse = (mean_betas[j] - true_val)^2 + variance,
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
results <- process_beta_data(combined_results)






