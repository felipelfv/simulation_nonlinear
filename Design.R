#### 1. General Session Information ####

# In this function, I have all the parameters (and variations) that are called
# in the Simulation.R script. This file needs to be sourced.

# Importantly, this script creates more conditions than necessary (or needed).
# The objective here was simply to have a general version (at the expense of 
# having something parsimonious)

# Possible source of confusion:
# The "latent_exo_distribution" refers to generating the exogenous latent variable
# as a normally distributed random variable or nonnormal. However, we have 
# two ways of generating nonnormal latent variables: (i) with skewness and kurtosis 
# passed to the rIG() function, and (ii) uniformly distributed (as in GAPI). 
# This is the reason for having "exo_methods", to differentiate both nonnormal 
# scenarios.Importantly, is that rIG() is also used to generate normal variable 
# as well, given the values for skewness and kurtosis. 

#### 2. Design ####

##### 2.1 Libraries used for  #####

library(covsim) # in GenerateData.R for rIG()
library(lavaan) # in GenerateData.R and Methods.R for SAM approach
library(modsem) # in Methods.R for LMS, QML, and UCA approach
# In Simulation.R for parallel computing:
library(foreach)
library(doParallel)
library(doRNG)


##### 2.2 Parameters that we vary in the simulation study #####

sample_sizes <- c(200L, 500L, 800L) 
reliability_values <- c(0.4, 0.6, 0.8)
population_models <- c("population.linear.model", "population.interaction.model", "population.full.model")
latent_exo_distribution <- c("normal", "nonnormal")
exo_methods <- c("rIG", "unif")
epsilon_distributions <- c("normal", "exp.rate1")


create_conditions <- function(sample_sizes, reliability_values,
                              population_models, latent_exo_distribution, 
                              exo_methods, epsilon_distributions) {
  # Base conditions
  base_conditions <- expand.grid(
    Population = population_models,
    Distribution = latent_exo_distribution,
    Exo_method = exo_methods,
    Epsilon = epsilon_distributions,
    N = sample_sizes,
    Rel = reliability_values,
    stringsAsFactors = FALSE
  )
  
  # Remove invalid combinations
  base_conditions <- base_conditions[!(base_conditions$Distribution == "normal" & 
                                         base_conditions$Exo_method == "unif"), ]
  
  result <- data.frame()
  
  # For each row in base_conditions
  for (i in 1:nrow(base_conditions)) {
    row <- base_conditions[i, ]
    
    # Which analysis models to use
    if (row$Population == "population.linear.model") {
      analysis_models <- c("fit.interaction.model", "fit.full.model")
    } else if (row$Population == "population.full.model") {
      analysis_models <- c("fit.full.model")
    } else {
      analysis_models <- c("fit.interaction.model")
    }
    
    # Create a row for each analysis model
    for (model in analysis_models) {
      new_row <- row
      new_row$Analysis_model <- model
      result <- rbind(result, new_row)
    }
  }
  row.names(result) <- NULL # Reset row indices
  result
}

conditions <- create_conditions(
  sample_sizes, 
  reliability_values, 
  population_models, 
  latent_exo_distribution, 
  exo_methods,
  epsilon_distributions
)

# For multiple clusters in HPC
#conditions$Seed <- 1:nrow(conditions)
# or this may be better: 
conditions$Seed <- (sample(1:1e9, size = nrow(conditions), replace = FALSE)) %% 1000

#RNGkind("L'Ecuyer-CMRG")
#Seeds <- list(.Random.seed)
#for(i in 2:216) {
#  Seeds[[i]] <-  nextRNGStream(Seeds[[i-1]])
#}



##### 2.3 Fixed parameters #####

# Parameters that remain constant (as of 04/03/25)
exo.mean <- rep(0, 2)
#target.var <- list("eta3" = 1.0) # target variance for eta 
R2 <- list("eta3" = 0.30)
rep <- 1000  # Repetitions (to be increased for the actual study): 1000(?)

