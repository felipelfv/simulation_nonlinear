# R^2 very close to 0.30
population_full_model <- "
# measurement model
eta1 =~ 1*x1 + 1*x2 + 1*x3
eta2 =~ 1*x4 + 1*x5 + 1*x6
eta3 =~ 1*x7 + 1*x8 + 1*x9
eta4 =~ 1*x10 + 1*x11 + 1*x12
eta5 =~ 1*x13 + 1*x14 + 1*x15

# structural model
eta3 ~ 0.3*eta1 + 0.3*eta2 + 0.15*eta1:eta2 + 0.1*eta1:eta1
eta4 ~ 0.25*eta1 + 0.25*eta3 + 0.12*eta1:eta3 + 0.08*eta2:eta2
eta5 ~ 0.2*eta2 + 0.2*eta3 + 0.2*eta4 + 0.1*eta3:eta4 + 0.1*eta3:eta3 + 0.08*eta4:eta4

# variances and covariances
eta1 ~~ 0.4*eta2
eta1 ~~ 1*eta1
eta2 ~~ 1*eta2

# residual variances
eta3 ~~ 0.8*eta3   
eta4 ~~ 0.75*eta4  
eta5 ~~ 1.25*eta5 

# measurement errors for reliability approx. 0.8
x1 ~~ 0.25*x1; x2 ~~ 0.25*x2; x3 ~~ 0.25*x3
x4 ~~ 0.25*x4; x5 ~~ 0.25*x5; x6 ~~ 0.25*x6
x7 ~~ 0.25*x7; x8 ~~ 0.25*x8; x9 ~~ 0.25*x9
x10 ~~ 0.24*x10; x11 ~~ 0.24*x11; x12 ~~ 0.24*x12
x13 ~~ 0.42*x13; x14 ~~ 0.42*x14; x15 ~~ 0.42*x15 
"

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
eta5 ~ eta2 + eta3 + eta4 + eta3:eta4 + eta3:eta3 
"

data_final <- GenerateData(
  model = population_full_model,
  N = 1000,
  skewness = c(0, 0),
  excesskurtosis = c(0, 0),
  distr.exo = "normal.rIG",
  distr.zeta = "normal",
  distr.epsilon = "normal",
  seed = 123,
  add.eta = FALSE,
  return.info = TRUE
)

fit <- sam(model = analysis_model, data = data_final)

library(lavaan)

# conditions
sample_sizes <- c(200, 400, 1000)
distributions <- c("normal.rIG", "uniform", "nonnormal.rIG")
n_sim <- 100

all_results <- list()

# simulations for each condition
for(N in sample_sizes) {
  for(distr in distributions) {
    
    condition_name <- paste0("N", N, "_", distr)
    cat("\n\nStarting condition:", condition_name, "\n")
    
    # storage for this condition
    condition_results <- list()
    
    for(i in 1:n_sim) {
      cat("\rIteration", i, "/", n_sim)
      
      # Generate data
      data <- GenerateData(
        model = population_full_model,
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

#saveRDS(all_results, "simulation_all_conditions.rds")

# SE/SD ratios for each condition
interaction_terms <- c("eta3 ~ eta1:eta2", "eta3 ~ eta1:eta1", 
                       "eta4 ~ eta1:eta3", "eta4 ~ eta2:eta2",
                       "eta5 ~ eta3:eta4", "eta5 ~ eta3:eta3")

summary_all <- data.frame()

for(N in sample_sizes) {
  for(distr in distributions) {
    
    condition_name <- paste0("N", N, "_", distr)
    results <- all_results[[condition_name]]
    
    # estimates and SEs
    n_params <- length(interaction_terms)
    estimates_matrix <- matrix(NA, nrow = n_sim, ncol = n_params)
    se_matrix <- matrix(NA, nrow = n_sim, ncol = n_params)
    
    for(i in 1:n_sim) {
      pe <- results[[i]]
      
      for(j in 1:n_params) {
        parts <- strsplit(interaction_terms[j], " ~ ")[[1]]
        lhs <- parts[1]
        rhs <- parts[2]
        
        row_idx <- which(pe$lhs == lhs & pe$op == "~" & pe$rhs == rhs)
        
        if(length(row_idx) > 0) {
          estimates_matrix[i, j] <- pe$est[row_idx]
          se_matrix[i, j] <- pe$se[row_idx]
        }
      }
    }
    
    empirical_sd <- apply(estimates_matrix, 2, sd, na.rm = TRUE)
    average_se <- apply(se_matrix, 2, mean, na.rm = TRUE)
    se_sd_ratio <- average_se / empirical_sd
    
    for(j in 1:n_params) {
      summary_all <- rbind(summary_all, data.frame(
        N = N,
        Distribution = distr,
        Parameter = interaction_terms[j],
        Avg_Estimate = mean(estimates_matrix[, j], na.rm = TRUE),
        Empirical_SD = empirical_sd[j],
        Average_SE = average_se[j],
        SE_SD_Ratio = se_sd_ratio[j]
      ))
    }
  }
}

(summary_all)
