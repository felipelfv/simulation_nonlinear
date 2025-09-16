# Define population models for each reliability level
population_model_rel08 <- "
# Measurement model
eta1 =~ 1*x1 + 1*x2 + 1*x3
eta2 =~ 1*x4 + 1*x5 + 1*x6
eta3 =~ 1*x7 + 1*x8 + 1*x9
eta4 =~ 1*x10 + 1*x11 + 1*x12
eta5 =~ 1*x13 + 1*x14 + 1*x15
# Structural model
eta3 ~ 0.3*eta1 + 0.3*eta2 + 0.15*eta1:eta2 + 0.1*eta1:eta1
eta4 ~ 0.25*eta1 + 0.25*eta3 + 0.12*eta1:eta3 + 0.08*eta2:eta2
eta5 ~ 0.2*eta2 + 0.2*eta3 + 0.2*eta4 + 0.1*eta3:eta4 + 0.1*eta3:eta3 + 0.08*eta4:eta4
# Variances and covariances
eta1 ~~ 0.4*eta2
eta1 ~~ 1*eta1
eta2 ~~ 1*eta2
eta3 ~~ 0.8*eta3
eta4 ~~ 0.75*eta4
eta5 ~~ 1.39*eta5
# Measurement errors for reliability ≈ 0.8
x1 ~~ 0.25*x1; x2 ~~ 0.25*x2; x3 ~~ 0.25*x3
x4 ~~ 0.25*x4; x5 ~~ 0.25*x5; x6 ~~ 0.25*x6
x7 ~~ 0.25*x7; x8 ~~ 0.25*x8; x9 ~~ 0.25*x9
x10 ~~ 0.24*x10; x11 ~~ 0.24*x11; x12 ~~ 0.24*x12
x13 ~~ 0.40*x13; x14 ~~ 0.40*x14; x15 ~~ 0.40*x15
"

population_model_rel06 <- "
# Measurement model
eta1 =~ 1*x1 + 1*x2 + 1*x3
eta2 =~ 1*x4 + 1*x5 + 1*x6
eta3 =~ 1*x7 + 1*x8 + 1*x9
eta4 =~ 1*x10 + 1*x11 + 1*x12
eta5 =~ 1*x13 + 1*x14 + 1*x15
# Structural model
eta3 ~ 0.3*eta1 + 0.3*eta2 + 0.15*eta1:eta2 + 0.1*eta1:eta1
eta4 ~ 0.25*eta1 + 0.25*eta3 + 0.12*eta1:eta3 + 0.08*eta2:eta2
eta5 ~ 0.2*eta2 + 0.2*eta3 + 0.2*eta4 + 0.1*eta3:eta4 + 0.1*eta3:eta3 + 0.08*eta4:eta4
# Variances and covariances
eta1 ~~ 0.4*eta2
eta1 ~~ 1*eta1
eta2 ~~ 1*eta2
eta3 ~~ 0.8*eta3
eta4 ~~ 0.75*eta4
eta5 ~~ 1.39*eta5
# Measurement errors for reliability ≈ 0.6
x1 ~~ 0.67*x1; x2 ~~ 0.67*x2; x3 ~~ 0.67*x3
x4 ~~ 0.67*x4; x5 ~~ 0.67*x5; x6 ~~ 0.67*x6
x7 ~~ 0.78*x7; x8 ~~ 0.78*x8; x9 ~~ 0.78*x9
x10 ~~ 0.71*x10; x11 ~~ 0.71*x11; x12 ~~ 0.71*x12
x13 ~~ 1.32*x13; x14 ~~ 1.32*x14; x15 ~~ 1.32*x15
"

population_model_rel04 <- "
# Measurement model
eta1 =~ 1*x1 + 1*x2 + 1*x3
eta2 =~ 1*x4 + 1*x5 + 1*x6
eta3 =~ 1*x7 + 1*x8 + 1*x9
eta4 =~ 1*x10 + 1*x11 + 1*x12
eta5 =~ 1*x13 + 1*x14 + 1*x15
# Structural model
eta3 ~ 0.3*eta1 + 0.3*eta2 + 0.15*eta1:eta2 + 0.1*eta1:eta1
eta4 ~ 0.25*eta1 + 0.25*eta3 + 0.12*eta1:eta3 + 0.08*eta2:eta2
eta5 ~ 0.2*eta2 + 0.2*eta3 + 0.2*eta4 + 0.1*eta3:eta4 + 0.1*eta3:eta3 + 0.08*eta4:eta4
# Variances and covariances
eta1 ~~ 0.4*eta2
eta1 ~~ 1*eta1
eta2 ~~ 1*eta2
eta3 ~~ 0.8*eta3
eta4 ~~ 0.75*eta4
eta5 ~~ 1.39*eta5
# Measurement errors for reliability ≈ 0.4
x1 ~~ 1.50*x1; x2 ~~ 1.50*x2; x3 ~~ 1.50*x3
x4 ~~ 1.50*x4; x5 ~~ 1.50*x5; x6 ~~ 1.50*x6
x7 ~~ 1.76*x7; x8 ~~ 1.76*x8; x9 ~~ 1.76*x9
x10 ~~ 1.61*x10; x11 ~~ 1.61*x11; x12 ~~ 1.61*x12
x13 ~~ 2.98*x13; x14 ~~ 2.98*x14; x15 ~~ 2.98*x15
"

# Define analysis model (same for all)
analysis_model <- "
# Measurement model
eta1 =~ x1 + x2 + x3
eta2 =~ x4 + x5 + x6
eta3 =~ x7 + x8 + x9
eta4 =~ x10 + x11 + x12
eta5 =~ x13 + x14 + x15
# Structural model
eta3 ~ eta1 + eta2 + eta1:eta2 + eta1:eta1
eta4 ~ eta1 + eta3 + eta1:eta3 + eta2:eta2
eta5 ~ eta2 + eta3 + eta4 + eta3:eta4 + eta3:eta3 + eta4:eta4
"

# conditions
sample_sizes <- c(200, 400, 1000)
distributions <- c("normal.rIG", "uniform", "nonnormal.rIG")
reliabilities <- c(0.8, 0.6, 0.4)
n_sim <- 2
all_results <- list()

# Create list of population models
population_models <- list(
  "0.8" = population_model_rel08,
  "0.6" = population_model_rel06,
  "0.4" = population_model_rel04
)

# simulations for each condition
for(N in sample_sizes) {
  for(distr in distributions) {
    for(rel in reliabilities) {
      
      condition_name <- paste0("N", N, "_", distr, "_rel", rel)
      cat("\n\nStarting condition:", condition_name, "\n")
      
      # Select appropriate population model
      population_model <- population_models[[as.character(rel)]]
      
      # storage for this condition
      condition_results <- list()
      
      for(i in 1:n_sim) {
        cat("\rIteration", i, "/", n_sim)
        
        # Generate data
        data <- GenerateData(
          model = population_model,
          N = N,
          seed = 123 + i,
          skewness = c(0, 0),
          excesskurtosis = c(0, 0),
          distr.exo = distr,
          distr.zeta = "normal",
          distr.epsilon = "normal",
          add.eta = FALSE,
          return.info = TRUE
        )
        
        fit <- sam(analysis_model, data = data, se = "twostep")
        
        condition_results[[i]] <- parameterEstimates(fit, remove.step1 = FALSE)
      }
      
      # results for this condition
      all_results[[condition_name]] <- condition_results
    }
  }
}