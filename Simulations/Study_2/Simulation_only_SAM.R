library(copula)
library(covsim)
#install.packages("copula")

GenerateData <- function(model, 
                         N = 1000L,
                         skewness = NULL,
                         excesskurtosis = NULL,
                         exo.mean = NULL,
                         distr.exo = "normal.rIG",
                         distr.zeta = "normal",
                         distr.epsilon = "normal",
                         center.exogenous.latent = FALSE,
                         center.exogenous.manifest = FALSE,
                         center.lv.dependent = FALSE,
                         center.lv.prod = FALSE,
                         center.indicators = FALSE,
                         seed = NULL,
                         add.eta = FALSE,
                         return.info = TRUE) {
  
  if(!is.null(seed)) set.seed(seed)
  
  ##============================================================================
  ## MODEL INFORMATION
  ##============================================================================
  
  # lavaan's parsed model information
  fit <- lavaan::sam(model)
  pt <- lavaan::parTable(fit)
  
  # structural model information
  model_info <- list(structural = list(
    dependent = character(), 
    exogenous = character(),
    equations = list(),
    coefficients = list(),
    interactions = list(),
    generation_order = NULL
  ))
  
  # process structural model
  model_info$structural$dependent <- unique(pt$lhs[pt$op == "~"])
  
  # process equations, coefficients, and interactions 
  for(dv in model_info$structural$dependent) {
    eq_rows <- pt$op == "~" & pt$lhs == dv
    model_info$structural$equations[[dv]] <- pt$rhs[eq_rows]
    model_info$structural$coefficients[[dv]] <- pt$est[eq_rows]
    
    # interactions 
    dv_interactions <- intersect(
      pt$rhs[eq_rows],
      c(fit@pta$vnames$lv.interaction[[1]], fit@pta$vnames$ov.interaction[[1]])
    )
    if (length(dv_interactions) > 0) {
      model_info$structural$interactions[[dv]] <- dv_interactions
    }
  }
  
  # all base variables
  base_vars <- c(
    fit@pta$vnames$lv.regular[[1]],
    setdiff(fit@pta$vnames$ov[[1]], fit@pta$vnames$ov.ind[[1]])
  )
  
  # "B" matrix
  B <- matrix(NA, nrow = length(base_vars), ncol = length(base_vars))
  rownames(B) <- colnames(B) <- base_vars
  
  # process each regression equation to get dependencies
  for(dv in model_info$structural$dependent) {
    rhs_terms <- pt$rhs[pt$op == "~" & pt$lhs == dv]
    
    for(term in rhs_terms) {
      components <- unlist(strsplit(term, ":"))
      for(comp in components) {
        if(comp %in% base_vars) {
          B[dv, comp] <- 1
        }
      }
    }
  }
  
  # ancestors 
  ancestors <- lavaan:::lav_utils_get_ancestors(B)
  names(ancestors) <- rownames(B)
  ancestor_lengths <- sapply(ancestors, length)
  # exogenous variables
  model_info$structural$exogenous <- names(ancestor_lengths)[ancestor_lengths == 0]
  # generation order
  model_info$structural$generation_order <- names(ancestor_lengths)[order(ancestor_lengths)]
  
  ##============================================================================
  ## EXTRACT ALL RESIDUAL VARIANCES FROM MODEL
  ##============================================================================
  
  # residual variances specified in the model
  residual_rows <- pt[pt$op == "~~" & pt$lhs == pt$rhs, ]
  all_residual_vars <- setNames(residual_rows$est, residual_rows$lhs)
  
  # distinguish structural and indicator residual variances
  indicator_vars <- fit@pta$vnames$ov.ind[[1]]
  structural_residual_vars <- all_residual_vars[!names(all_residual_vars) %in% indicator_vars]
  indicator_residual_vars <- all_residual_vars[names(all_residual_vars) %in% indicator_vars]
  
  ##============================================================================
  ## COPULA-BASED DATA GENERATION FUNCTION (FROM ORIGINAL AUTHORS)
  ##============================================================================
  
  CopSEM <- function(copmvdc, Sigma, nw = 100000, np = 1000,
                     shift = rep(0, ncol(Sigma))) {
    ## copmvdc ... joint density from mvdc()
    ## Sigma ... model VC-matrix to be approximated
    ## nw ... sample size for warm-up sample
    ## np ... sample size for production sample
    Xw <- rMvdc(nw, copmvdc) ## draw warm-up sample
    Sw <- cov(Xw) ## warm-up VC matrix
    Sigma.eigen <- eigen(Sigma) ## EV decomposition Sigma
    Sigmaroot <- Sigma.eigen$vectors %*% sqrt(diag(Sigma.eigen$values)) %*%
      t(Sigma.eigen$vectors) ## root Sigma
    Sx.eigen <- eigen(solve(Sw)) ## EV decomposition S
    Sxroot <- Sx.eigen$vectors %*% sqrt(diag(Sx.eigen$values)) %*%
      t(Sx.eigen$vectors) ## root S
    X <- rMvdc(np, copmvdc) ## draw production sample
    X <- sweep(X, 2, shift, FUN = "+")
    Y <- (X %*% (Sxroot) %*% Sigmaroot) ## linear combination for Y
    list(Y = Y, covY = (cov(Y))) ## return Y and cov(Y)
  }
  
  ##============================================================================
  ## DATA GENERATION
  ##============================================================================
  
  # matrix for all variables
  all_vars <- c(fit@pta$vnames$lv.regular[[1]],
                fit@pta$vnames$lv.interaction[[1]],
                setdiff(fit@pta$vnames$ov[[1]], 
                        fit@pta$vnames$ov.ind[[1]]))
  
  Values <- matrix(NA, nrow = N, ncol = length(all_vars))
  colnames(Values) <- all_vars
  
  # GENERATE ALL EXOGENOUS VARIABLES AT ONCE 
  exo_vars <- model_info$structural$exogenous
  psi_matrix <- lavaan::lavInspect(fit, "est")$psi
  exo.vcov <- psi_matrix[exo_vars, exo_vars, drop = FALSE]
  
  for(i in 1:nrow(exo.vcov)) { 
    if(is.na(exo.vcov[i,i])) {
      exo.vcov[i,i] <- 1
    }
  }
  
  # generate exogenous variables
  generate_exo <- function() {
    if(distr.exo == "unif") {
      # copula method similar to lonati et al, 2024
      coppar <- gumbelCopula(2, dim = length(exo_vars))
      
      # mvdc object with uniform marginals
      copjoint <- mvdc(coppar, 
                       margins = rep("unif", length(exo_vars)),
                       paramMargins = rep(list(list(min = -2.5, max = 2.5)), 
                                          length(exo_vars)))
      
      # using CopSEM to get correct covariance structure
      result <- CopSEM(copjoint, exo.vcov, nw = 10000, np = N)
      EXO <- result$Y
      
    } else {
      EXO <- covsim::rIG(N, sigma = exo.vcov, skewness = skewness, 
                         excesskurtosis = excesskurtosis)[[1]]
    }
    EXO
  }
  
  # generate initial EXO
  EXO <- generate_exo()
  
  # ensure it's a matrix
  if(!is.matrix(EXO)) {
    EXO <- as.matrix(EXO)
  }
  
  colnames(EXO) <- exo_vars
  
  # center exogenous variables
  eta_cols <- grep("eta", exo_vars, value = TRUE)
  if(center.exogenous.latent && length(eta_cols) > 0) {
    EXO[, eta_cols] <- t(t(EXO[, eta_cols]) - colMeans(EXO[, eta_cols]))
  } else if(!is.null(exo.mean) && length(eta_cols) > 0) {
    EXO[, eta_cols] <- t(t(EXO[, eta_cols]) + exo.mean[1:length(eta_cols)])
  }
  
  manifest_cols <- setdiff(exo_vars, eta_cols)
  if(center.exogenous.manifest && length(manifest_cols) > 0) {
    EXO[, manifest_cols] <- t(t(EXO[, manifest_cols]) - colMeans(EXO[, manifest_cols]))
  } else if(!is.null(exo.mean) && length(manifest_cols) > 0) {
    offset = length(eta_cols)
    EXO[, manifest_cols] <- t(t(EXO[, manifest_cols]) + 
                                exo.mean[(offset + 1):length(exo.mean)])
  }
  
  Values[, exo_vars] <- EXO
  
  observed_R2 <- list()
  deterministic_vars <- list()
  
  # GENERATE DEPENDENT VARIABLES WITH FIXED RESIDUAL VARIANCES
  intercepts <- lavaan::lavInspect(fit, "est")$alpha
  for(var in model_info$structural$generation_order) { 
    if(var %in% model_info$structural$dependent) {
      
      terms <- model_info$structural$equations[[var]]
      # handle interactions first
      if(!is.null(model_info$structural$interactions[[var]])) {
        for(inter in model_info$structural$interactions[[var]]) {
          if(all(is.na(Values[, inter]))) {
            components <- unlist(strsplit(inter, ":"))
            Values[, inter] <- Values[, components[1]] * Values[, components[2]]
            
            if(center.lv.prod) {
              Values[, inter] <- Values[, inter] - mean(Values[, inter])
            }
          }
        }
      }
      
      # coefficients
      equation_coefs <- model_info$structural$coefficients[[var]]
      
      # deterministic part
      deterministic_part <- Values[, terms, drop = FALSE] %*% equation_coefs
      
      if (!is.null(intercepts) && !is.na(intercepts[var,1])) {
        deterministic_part <- intercepts[var,1] + deterministic_part
      }
      
      # variance of deterministic part
      var_det <- var(as.vector(deterministic_part))
      deterministic_vars[[var]] <- var_det
      
      # fixed residual variance from model
      resid_var <- structural_residual_vars[var]
      if(is.na(resid_var)) {
        resid_var <- 1
        warning(paste("No residual variance found for", var, "- using 1"))
      }
      
      # generate residual with fixed variance
      zeta <- switch(distr.zeta,
                     "normal" = rnorm(N, 0, sqrt(resid_var)),
                     "exp.rate1" = rexp(N, rate = 1/sqrt(resid_var)) - sqrt(resid_var),
                     stop(paste("wrong option for distr.zeta:", distr.zeta)))
      
      Values[, var] <- deterministic_part + zeta
      
      if(center.lv.dependent) {
        Values[, var] <- t(t(Values[, var]) - mean(Values[, var]))
      }
      
      # observed R^2
      total_var <- var(Values[, var])
      observed_R2[[var]] <- var_det / total_var
    }
  }
  
  ##============================================================================
  ## MEASUREMENT PART - WITH FIXED VARIANCES FROM MODEL
  ##============================================================================
  
  lambda <- lavaan::lavInspect(fit, "est")$lambda
  LAMBDA <- lambda[startsWith(rownames(lambda), "x"), fit@pta$vnames$lv.regular[[1]], 
                   drop = FALSE]
  
  # measurement errors with fixed variances from model
  n_indicators <- nrow(LAMBDA)
  THETA <- matrix(NA, nrow = N, ncol = n_indicators)
  colnames(THETA) <- rownames(LAMBDA)
  
  for(i in seq_len(n_indicators)) {
    ind_name <- rownames(LAMBDA)[i]
    # variance specified in the model
    ind_var <- indicator_residual_vars[ind_name]
    if(is.na(ind_var)) {
      ind_var <- 1
      warning(paste("No residual variance found for", ind_name, "- using 1"))
    }
    
    if(distr.epsilon == "normal") {
      THETA[, i] <- rnorm(N, 0, sqrt(ind_var))
    } else {
      THETA[, i] <- rexp(N, rate = 1/sqrt(ind_var)) - sqrt(ind_var)
    }
  }
  
  # indicators: Y = Lambda * Eta + Theta
  Y <- as.matrix(Values[, fit@pta$vnames$lv.regular[[1]]]) %*% t(LAMBDA) + THETA 
  colnames(Y) <- rownames(LAMBDA)
  
  # center indicators 
  if(center.indicators) {
    Y <- t(t(Y) - colMeans(Y))
  }
  
  # RESULTING reliability (not used for generation, just for info)
  observed_reliabilities <- list()
  if(return.info) {
    for(i in seq_along(fit@pta$vnames$lv.regular[[1]])) {
      eta <- fit@pta$vnames$lv.regular[[1]][i]
      eta_var <- var(Values[, eta])
      
      ind_idx <- which(LAMBDA[, i] != 0)
      if(length(ind_idx) > 0) {
        reliabilities <- numeric(length(ind_idx))
        for(j in seq_along(ind_idx)) {
          ind_name <- rownames(LAMBDA)[ind_idx[j]]
          loading <- LAMBDA[ind_idx[j], i]
          error_var <- indicator_residual_vars[ind_name]
          # what reliability we get from fixed variances
          reliabilities[j] <- (loading^2 * eta_var) / (loading^2 * eta_var + error_var)
        }
        names(reliabilities) <- rownames(LAMBDA)[ind_idx]
        observed_reliabilities[[eta]] <- reliabilities
      }
    }
  }
  
  if(add.eta) {
    Results <- cbind(Values, Y)
  } else {
    manifest_vars <- setdiff(fit@pta$vnames$ov[[1]], fit@pta$vnames$ov.ind[[1]])
    manifest_vars <- manifest_vars[!grepl("eta", manifest_vars)]
    Values_indicators <- Values[, manifest_vars, drop = FALSE]
    Results <- cbind(Values_indicators, Y)
  }
  
  result <- as.data.frame(Results)
  
  if(return.info) {
    attr(result, "observed_R2") <- observed_R2
    attr(result, "observed_reliabilities") <- observed_reliabilities
    attr(result, "fixed_residual_variances") <- all_residual_vars
    attr(result, "model_info") <- model_info
  }
  
  result
}

######## Simulation trial for study 2 ###########

complex_model <- "
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

# Optimized residual variances
eta3 ~~ 0.8*eta3   # R² ≈ 0.32
eta4 ~~ 0.75*eta4  # R² ≈ 0.30
eta5 ~~ 1.25*eta5  # Slight increase to bring R² closer to 0.30

# Measurement errors for reliability ≈ 0.8
x1 ~~ 0.25*x1; x2 ~~ 0.25*x2; x3 ~~ 0.25*x3
x4 ~~ 0.25*x4; x5 ~~ 0.25*x5; x6 ~~ 0.25*x6
x7 ~~ 0.25*x7; x8 ~~ 0.25*x8; x9 ~~ 0.25*x9
x10 ~~ 0.24*x10; x11 ~~ 0.24*x11; x12 ~~ 0.24*x12
x13 ~~ 0.42*x13; x14 ~~ 0.42*x14; x15 ~~ 0.42*x15  # Slight increase
"

# Test the model with your function
data <- GenerateData(
  model = complex_model,
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
# Check the R² values for all three dependent variables
print("Observed R² values:")
print(attr(data, "observed_R2"))

# Check reliabilities
print("\nObserved reliabilities:")
print(attr(data, "observed_reliabilities"))

#####################

# Expanded simulation with varying conditions
library(lavaan)

# Define models
complex_model <- "
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

# Optimized residual variances
eta3 ~~ 0.8*eta3   # R² ≈ 0.32
eta4 ~~ 0.75*eta4  # R² ≈ 0.30
eta5 ~~ 1.25*eta5  # Slight increase to bring R² closer to 0.30

# Measurement errors for reliability ≈ 0.8
x1 ~~ 0.25*x1; x2 ~~ 0.25*x2; x3 ~~ 0.25*x3
x4 ~~ 0.25*x4; x5 ~~ 0.25*x5; x6 ~~ 0.25*x6
x7 ~~ 0.25*x7; x8 ~~ 0.25*x8; x9 ~~ 0.25*x9
x10 ~~ 0.24*x10; x11 ~~ 0.24*x11; x12 ~~ 0.24*x12
x13 ~~ 0.42*x13; x14 ~~ 0.42*x14; x15 ~~ 0.42*x15  # Slight increase
"

fit_model <- "
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

# conditions
sample_sizes <- c(200, 400, 1000)
distributions <- c("normal.rIG", "unif", "nonnormal.rIG")
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
      
      # generate data
      data <- GenerateData(
        model = complex_model,
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
      
      # fit model
      fit <- sam(fit_model, data = data, se = "twostep")
      
      # save parameter estimates
      condition_results[[i]] <- parameterEstimates(fit, remove.step1 = FALSE)
    }
    
    # store results for this condition
    all_results[[condition_name]] <- condition_results
  }
}

#saveRDS(all_results, "simulation_all_conditions.rds")

# SE/SD ratios for each condition
interaction_terms <- c("eta3 ~ eta1:eta2", "eta3 ~ eta1:eta1", 
                       "eta4 ~ eta1:eta3", "eta4 ~ eta2:eta2",
                       "eta5 ~ eta3:eta4", "eta5 ~ eta3:eta3")

# table for all conditions
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
    
    # compute statistics
    empirical_sd <- apply(estimates_matrix, 2, sd, na.rm = TRUE)
    average_se <- apply(se_matrix, 2, mean, na.rm = TRUE)
    se_sd_ratio <- average_se / empirical_sd
    
    # add to summary
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



