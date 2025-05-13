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

set.seed(123)
conditions$Seed <- (sample(1:1e9, size = nrow(conditions), replace = FALSE))

# The only check performed here was to ensure seeds are not the same
# The iterations per condition make use of streams

# n of unique seeds equals the total number of rows
all_unique <- length(unique(conditions$Seed)) == nrow(conditions)
(all_unique)  # TRUE

# any duplicated values?
any_duplicates <- any(duplicated(conditions$Seed))
print(any_duplicates)  # FALSE
