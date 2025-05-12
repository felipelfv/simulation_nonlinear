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

##### 2.2 Fixed parameters in study 2 #####

set.seed(123)

#target.var <- list("eta4" = 1.0) # target variance for eta 
R2 <- list("eta4" = 0.30)

skewness <- rep(0, 4)
excesskurtosis <- rep(0, 4)
# If nonnormal:
#skewness <- rep(2, 4)
#excesskurtosis <- rep(7, 4)

exo.mean <- rep(0, 4)

N <- 1000L
Rel <- 0.8


