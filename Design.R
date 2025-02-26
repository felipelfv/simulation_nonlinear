#### 1. General Session Information ####

# In this function, I have all the parameters (and variations) that are called
# in the Simulation.R script. This file needs to be sourced.

#### 2. Design ####

N <- c(400L, 600L, 800L) # For the different sample sizes (i.e., dataset size)
RELIABILITY <- c(0.2, 0.6, 0.8) # For the indicators (i.e., measurement model)

# Conditions as in Brand et al. (2014) - not fully factorial
Conditions <- expand.grid(
  POPULATION      = c("population.linear.model", "population.interaction.model", "population.full.model"),
  DISTRIBUTION   = c("normal", "nonnormal"),
  EXO_METHOD     = c("rIG", "unif"),  # Factor for exogenous generation method
  N              = N,
  REL            = RELIABILITY,
  stringsAsFactors = FALSE
)

Conditions$analysis_model <- ifelse(
  Conditions$POPULATION == "population.full.model",
  "fit.full.model",
  "fit.interaction.model"
)

# Remove invalid combinations:
# If distribution is "normal", exo_method should only be "rIG"
Conditions <- Conditions[!(Conditions$DISTRIBUTION == "normal" & 
                             Conditions$EXO == "unif"), ]

# Add epsilon distribution based on the distribution factor
Conditions$epsilon <- ifelse(Conditions$DISTRIBUTION == "normal", 
                             "normal", "exp.rate1")
