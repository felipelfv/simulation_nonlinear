#### 1. General Session Information ####

# ADD NUMBER OF PARAMETERS HERE FOR THE LMS AND QML. 
# This will make it easier in the simulation loop for the correct results 

#### 2. Design ####

##### 2.1 Libraries used for  #####

library(covsim) # in GenerateData.R for rIG()
library(lavaan) # in GenerateData.R and Methods.R for SAM approach
library(modsem) # in Methods.R for LMS, QML, and UCA approach
# In Simulation.R for parallel computing:
library(foreach)
library(doParallel)
library(doRNG)

##### 2.2 Fixed parameters in study 2 #####

##### 2.3 Design for model comparison #####
create_conditions <- function() {
  
  model_numbers <- 1:7  # 1 through 7
  
  # fixed parameters for all conditions
  N <- 1000L
  Rel <- 0.8
  
  result <- data.frame()
  
  for (model_num in model_numbers) {
    # number of exogenous variables based on model
    num_exo_vars <- ifelse(model_num <= 2, 3,  # Models 1-2 have 3 exo variables
                           ifelse(model_num == 3, 4,  # Model 3 has 4 exo variables
                                  5))  # Models 4-7 have 5 exo variables
    
    new_row <- data.frame(
      Population = paste0("population.large.model.", model_num),
      Distribution = "normal",  # (skewness=0, kurtosis=0)
      Exo_method = "rIG",
      Epsilon = "normal",
      N = N,
      Rel = Rel,
      R2 = 0.30,
      Analysis_model = paste0("fit.large.model.", model_num),
      Num_exo_vars = num_exo_vars,  # number of exogenous variables
      stringsAsFactors = FALSE
    )
    result <- rbind(result, new_row)
  }
  
  row.names(result) <- NULL
  result
}

conditions_2 <- create_conditions()

#target.var <- list("eta4" = 1.0) # target variance for eta 
#R2 <- list("eta4" = 0.30)
#skewness <- rep(0, 5)
#excesskurtosis <- rep(0, 5)
#If nonnormal: skewness <- rep(2, 4) and excesskurtosis <- rep(7, 4)
#exo.mean <- rep(0, 5)



