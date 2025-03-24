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

all_beta_averages <- function(all_results) {

  all_avgs <- list()
  
  # each element in all_results
  for(i in seq_along(all_results)) {
    data <- all_results[[i]]
    result <- list()
    
    # Process standard methods (lms, qml, uca)
    for(method in c("lms", "qml", "uca")) {
      if(!is.null(data$results[[method]])) {
        result[[method]] <- colMeans(data$results[[method]][,,"beta"])
      }
    }
    
    # Eventually all together when SE work
    if(!is.null(data$results$sam)) {
      result$sam <- colMeans(data$results$sam)
    }
    
    all_avgs[[i]] <- result
  }
  
  all_avgs
}

beta_avgs <- all_beta_averages(all_results)

