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

##### 2.2 Parameters that we vary in the simulation study #####

sample_sizes <- c(200L, 500L, 800L, 1000L) 
reliability_values <- c(0.8)
population_models <- c("population.large.model", "population.large.model.null")
latent_exo_distribution <- c("normal", "nonnormal")
epsilon_distributions <- c("normal", "exp.rate1")
create_conditions <- function(sample_sizes, reliability_values,
                              population_models, latent_exo_distribution, 
                              epsilon_distributions) {
  # Base conditions
  base_conditions <- expand.grid(
    Population = population_models,
    Distribution = latent_exo_distribution,
    Epsilon = epsilon_distributions,
    N = sample_sizes,
    Rel = reliability_values,
    stringsAsFactors = FALSE
  )

  result <- data.frame()
  
  # each row in base_conditions
  for (i in 1:nrow(base_conditions)) {
    row <- base_conditions[i, ]
    new_row <- row
    new_row$Analysis_model <- "fit.large.model"
    result <- rbind(result, new_row)
  }
  
  row.names(result) <- NULL # reset row indices
  result
}

conditions <- create_conditions(
  sample_sizes, 
  reliability_values, 
  population_models, 
  latent_exo_distribution, 
  epsilon_distributions
)

set.seed(123)
conditions$Seed <- (sample(1:1e9, size = nrow(conditions), replace = FALSE))

##### 2.3 Fixed parameters #####

exo.mean <- rep(0, 4)

R2 = list(
  "norm" = 0.4,
  "eta3" = 0.5,
  "eta4" = 0.2,
  "eta5" = 0.6
)
