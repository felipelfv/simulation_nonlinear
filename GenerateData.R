GenerateData <- function(model, 
                         N = 1000L,
                         skewness = NULL,
                         excesskurtosis = NULL,
                         exo.mean = NULL,
                         distr.exo = "rIG", # If unif, then generates naive unif, else uses rIG with excesskurtosis and skewness values
                         distr.zeta = "normal",
                         distr.epsilon = "normal",
                         center.exogenous.latent = TRUE,
                         center.exogenous.manifest = TRUE,
                         center.lv.dependent = FALSE,
                         center.lv.prod = FALSE,
                         center.indicators = FALSE,
                         target.var = NULL,
                         R2 = NULL,
                         rel = 0.64,
                         seed = NULL,
                         add.eta = FALSE) {
  
  if(!is.null(seed)) set.seed(seed)
  
  ###################################
  ######## Model Information ########
  ###################################
  
  # Lavaan's parsed model information
  fit <- lavaan::sam(model) # Fit the model without any data
  pt <- lavaan::parTable(fit) # Parameters table of the fitted model
  
  # Structural model information
  model_info <- list(structural = list(
    dependent = character(), 
    exogenous = character(), # "Pure" independent
    equations = list(),
    coefficients = list(),
    interactions = list(),
    generation_order = NULL # Important for tracing the path of each variable
  ))
  
  # Process structural model
  model_info$structural$dependent <- unique(pt$lhs[pt$op == "~"]) # Store all dependent variables
  # Process equations, coefficients, and interactions 
  for(dv in model_info$structural$dependent) { # For each dependent variable
    eq_rows <- pt$op == "~" & pt$lhs == dv # Select the predictors
    # Store equations and coefficients for each dependent 
    model_info$structural$equations[[dv]] <- pt$rhs[eq_rows]
    model_info$structural$coefficients[[dv]] <- pt$est[eq_rows] # Coefficients associated with each predictor
    # Interactions 
    dv_interactions <- intersect( # Rows in both x and y 
      pt$rhs[eq_rows],
      c(fit@pta$vnames$lv.interaction[[1]], fit@pta$vnames$ov.interaction[[1]])
    )
    if (length(dv_interactions) > 0) { # Just to keep the output clean from the parsed model
      model_info$structural$interactions[[dv]] <- dv_interactions # So, if no interaction, nothing instead of character(0) showing
    }
  }
  # All base variables (main and involved in interactions) - this is relevant for the "B" matrix(!)
  # Reasoning: there could be interaction among observed and latent, there could be interaction among latents that are not single predictors for any dependent
  base_vars <- c(
    fit@pta$vnames$lv.regular[[1]], # Latent variables
    setdiff(fit@pta$vnames$ov[[1]], fit@pta$vnames$ov.ind[[1]]) # Observed variables
  )
  
  # "B" matrix instead of lavaan::lavInspect(fit, "est")$beta
  # Reason: Im not interested in the coefficients. Im just interested about which variables are "connected". The "1" will show me that
  B <- matrix(NA, nrow = length(base_vars), ncol = length(base_vars)) # NA´s better than 0´s
  rownames(B) <- colnames(B) <- base_vars
  
  # Process each regression equation to get dependencies
  for(dv in model_info$structural$dependent) { # obtained from above with unique(pt$lhs[pt$op == "~"])
    rhs_terms <- pt$rhs[pt$op == "~" & pt$lhs == dv] # Predictors for each dependent
    
    for(i in rhs_terms) { # For all predictors for each dependent 
      # Split the term if interaction (very important!)
      components <- unlist(strsplit(rhs_terms, ":"))
      for(comp in components) { # And add the dependency information with 1
        if(comp %in% base_vars) { # So we get main effects and any variables involved in interactions
          B[dv, comp] <- 1 # Those will have 1
        }
      }
    }
  }
  # Ancestors 
  ancestors <- lavaan:::lav_utils_get_ancestors(B) # List for each variable in the "B" matrix
  names(ancestors) <- rownames(B)
  ancestor_lengths <- sapply(ancestors, length) # This is the crucial part that works perfectly with lav_utils_get_ancestors 
  # Exogenous variables
  model_info$structural$exogenous <- names(ancestor_lengths)[ancestor_lengths == 0] # 0 ancestors
  # Generation order
  model_info$structural$generation_order <- names(ancestor_lengths)[order(ancestor_lengths)] # based on ancestors order given by the length
  
  #################################
  ######## Data Generation ########
  #################################
  
  # Matrix for all variables
  # Maybe just $eqs.y[[1]] and $eqs.x[[1]] is enough?
  # Also add indicators already?
  all_vars <- c(fit@pta$vnames$lv.regular[[1]], # Regular latent variables
                fit@pta$vnames$lv.interaction[[1]], # Interaction terms
                setdiff(fit@pta$vnames$ov[[1]], 
                        fit@pta$vnames$ov.ind[[1]])) # Observed variables
  
  # Also, matrix instead of data frame: to be confirmed with Yves
  Values <- matrix(NA, nrow = N, ncol = length(all_vars))
  colnames(Values) <- all_vars
  # Could also just skip the matrix initialization here actually..
  
  # GENERATE ALL EXOGENOUS VARIABLES AT ONCE WITH rIG
  # Important: Add check related to variance as Yves has (comment in 06/12/24)
  exo_vars <- model_info$structural$exogenous
  psi_matrix <- lavaan::lavInspect(fit, "est")$psi # For (co)variances among exogenous
  # Extract the exogenous variables in PSI
  exo.vcov <- psi_matrix[exo_vars, exo_vars, drop = FALSE]
  
  for(i in 1:nrow(exo.vcov)) { 
    if(is.na(exo.vcov[i,i])) { # If variance left unspecified
      exo.vcov[i,i] <- 1 # We set to 1. Otherwise, the value given in the model is retained
    }
  }
  
  generate_exo <- function() {
    if(distr.exo == "unif") {
      # "NAIVE" approach; ignoring copulas
      # For uniform distribution as in GAPI article
      Z <- MASS::mvrnorm(N, mu = rep(0, ncol(exo.vcov)), Sigma = exo.vcov)
      EXO <- pnorm(Z) # Transform into an uniform distribution
      sd_desired <- sqrt(diag(exo.vcov)) # sd (from diag of cov matrix)
      scaling <- sqrt(12) * sd_desired
      EXO <- EXO - 0.5 # Center all variables by subtracting 0.5
      # Scale by multiplying each column by its scaling factor:
      EXO <- sweep(EXO, MARGIN = 2, STATS = scaling, FUN = "*")
    } else {
      EXO <- covsim::rIG(N, sigma = exo.vcov, skewness = skewness, 
                         excesskurtosis = excesskurtosis)[[1]] # Correlations are as given now!
    }
    EXO
  }
  
  EXO <- generate_exo()
  # Check variances: if they are too large (>2times the intended value), generate again
  if(any(apply(EXO, 2, var) > 2 * diag(exo.vcov))) {
    EXO <- generate_exo()
  }
  
  colnames(EXO) <- exo_vars
  
  # Would we ever center one and not the other? I believe rarely but I have both options 
  eta_cols <- grep("eta", exo_vars, value = TRUE)
  if(center.exogenous.latent) {
    EXO[, eta_cols] <- t(t(EXO[, eta_cols]) - colMeans(EXO[, eta_cols]))
  } else {
    EXO[, eta_cols] <- t(t(EXO[, eta_cols]) + exo.mean[1:length(eta_cols)])
  }
  
  manifest_cols <- setdiff(exo_vars, eta_cols) # All rows in x not in y
  if(center.exogenous.manifest) {
    EXO[, manifest_cols] <- t(t(EXO[, manifest_cols]) - colMeans(EXO[, manifest_cols]))
  } else {
    EXO[, manifest_cols] <- t(t(EXO[, manifest_cols]) + # For below, we assume latent and then manifest
                                exo.mean[(length(eta_cols) + 1):length(exo.mean)]) 
  }
  
  Values[, exo_vars] <- EXO
  
  # GENERATE DEPENDENT VARIABLES ACCORDING TO THE DEPENDENCIES 
  intercepts <- lavaan::lavInspect(fit, "est")$alpha # Intercepts for the dependent variables equations
  # Remember [[var]] read as string [["var"]]
  for(var in model_info$structural$generation_order) { 
    if(var %in% model_info$structural$dependent) { # Checks if current variable is a dependent variable
      terms <- model_info$structural$equations[[var]] # Equation terms for current variable
      # Handle interactions first
      if(!is.null(model_info$structural$interactions[[var]])) { # Checks if variable has interaction terms
        for(inter in model_info$structural$interactions[[var]]) {
          if(all(is.na(Values[, inter]))) {  # Checks if interaction term values are missing in the dataset
            components <- unlist(strsplit(inter, ":")) # Split the interaction into their individual variables
            # Multiplication for interactions
            # Check this if correctly applied; maybe should be the transpose * X? I dont think so; hence, just element-wise.
            Values[, inter] <- Values[, components[1]] * Values[, components[2]]
            
            # Check this:
            if(center.lv.prod) {
              Values[, inter] <- Values[, inter] - mean(Values[, inter]) # colMeans wont work with 1 column
            }
          }
        }
      }
      # Check coefficients
      if(!is.null(model_info$structural$coefficients[[var]])) { 
        equation_coefs <- model_info$structural$coefficients[[var]] # Coefficients for current variable
        if(length(equation_coefs) != length(terms)) {
          stop(sprintf("Error for %s: Number of coefficients (%d) does not match number of terms (%d)", 
                       var, length(equation_coefs), length(terms)))
        }
      } else {
        stop(sprintf("No coefficients provided for %s. Coefficients must be specified.", var))
      }
      
      # Check if all terms exist in the matrix
      missing_terms <- terms[!terms %in% colnames(Values)]
      if(length(missing_terms) > 0) {
        stop(sprintf("Terms not found in data for variable '%s': %s", 
                     var, paste(missing_terms, collapse = ", ")))
      }
      
      # Matrix algebra for the deterministic part as opposed to before (element-wise):
      deterministic_part <- Values[, terms, drop = FALSE] %*% equation_coefs
      if (!is.null(intercepts) && !is.na(intercepts[var,1])) {
        deterministic_part <- intercepts[var,1] + deterministic_part 
      }
      # Calculate variance without the residual/error
      var.nozeta <- var(deterministic_part)
      
      # Calculate target variance (this is given by the user)
      if(!is.null(target.var[[var]])) {
        if(var.nozeta < target.var[[var]]) {
          current_target_var <- target.var[[var]] - var.nozeta
        } else {
          warning(sprintf(
            "var.nozeta [%f] is larger than target.var [%f] for %s -- 
            using R2 instead",
            var.nozeta, target.var[[var]], var
          ))
          if(is.null(R2[[var]])) {
            stop(sprintf("No R2 value provided for %s. R2 must be specified.", 
                         var))
          }
          current_target_var <- var.nozeta * ((1-R2[[var]])/R2[[var]])
        }
      } else {
        if(var.nozeta < 0.0001) {
          current_target_var <- 1
        } else {
          if(is.null(R2[[var]])) {
            stop(sprintf("No R2 value provided for %s. R2 must be specified.", 
                         var))
          }
          current_target_var <- var.nozeta * ((1-R2[[var]])/R2[[var]])
        }
      }
      
      # Generate and standardize zeta 
      zeta <- switch(distr.zeta,
                     "normal" = rnorm(N, 0, sqrt(current_target_var)),
                     "exp.rate1" = rexp(N, rate = 1/sqrt(current_target_var)) - 1,
                     stop(paste("wrong option for distr.zeta:", distr.zeta)))
      
      zeta <- zeta - mean(zeta)
      zeta <- (zeta - mean(zeta)) * sqrt(c(current_target_var) / c(var(zeta - mean(zeta))))
      # `c()` otherwise warning about deprecated 
      
      # Entries for the variables with the residual
      Values[, var] <- deterministic_part + zeta # Add on top of the deterministic part
      
      if(center.lv.dependent) {
        Values[, var] <- t(t(Values[, var]) - mean(Values[, var])) # Same remark as with lv.prod:
      }
    }
  }
  
  # GENERATE THE MEASUREMENT PART
  # Get LAMBDA and calculate values for each pure eta
  lambda <- lavaan::lavInspect(fit, "est")$lambda
  #pure_etas <- fit@pta$vnames$lv.regular[[1]]
  LAMBDA <- lambda[startsWith(rownames(lambda), "x"), fit@pta$vnames$lv.regular[[1]], 
                   drop = FALSE] # Limitation: assumes indicators as x only
  # The above to drop the manifest variables from rows
  
  # Calculates target variances for the indicators
  # Could also use values from above for loop but for now just making sure things work as expected
  # We first get the variance of the latent variables
  eta_vars <- sapply(fit@pta$vnames$lv.regular[[1]], 
                     function(eta) var(Values[, eta])) # Correct to use N-1 as default in `var()`
  target_var_indicators <- eta_vars * (1/rel - 1)
  
  # Generate errors in theta 
  n_indicators <- nrow(LAMBDA)
  # For the different factor models below: 
  indicator_groups <- lapply(seq_along(fit@pta$vnames$lv.regular[[1]]),
                             function(i) which(LAMBDA[, i] == 1)) # Position of elements (logical vector)
  
  THETA <- matrix(NA, nrow = N, ncol = n_indicators)
  colnames(THETA) <- rownames(LAMBDA)
  
  # Get THETA matrix by the indicators
  # This part is the only part I am not 100% if mathematically correct (in accordance to SEM)
  # Reasoning: All indicators measuring the same latent variable get errors drawn from the same distribution (with the same variance)
  for(i in seq_along(indicator_groups)) { # For each latent variable,
    inds <- indicator_groups[[i]] # We get the indicators for it
    t_var_ind <- target_var_indicators[i] # We get the target variance from above calculation
    
    THETA[, inds] <- if(distr.epsilon == "normal") {
      matrix(rnorm(N * length(inds), 0, sqrt(t_var_ind)), N, length(inds))
    } else {
      matrix(rexp(N * length(inds), rate = 1/sqrt(t_var_ind)) - sqrt(t_var_ind), N, length(inds))
    }
  }
  
  # To get the actual values for indicators: Y = eta %*% t(LAMBDA) + THETA
  Y <- as.matrix(Values[, fit@pta$vnames$lv.regular[[1]]]) %*% t(LAMBDA) + THETA 
  colnames(Y) <- rownames(LAMBDA)
  
  # Center indicators 
  if(center.indicators) {
    Y <- t(t(Y) - colMeans(Y))
  }
  
  # Keep etas or not in the final dataset; relevant for fitting model in lavaan
  if(add.eta) {
    Results <- cbind(Values, Y)
  } else {
    # Get manifest variables from the model excluding indicators
    manifest_vars <- setdiff(fit@pta$vnames$ov[[1]], fit@pta$vnames$ov.ind[[1]])
    # Remove interaction terms and variables latent variables (as in eta..)
    manifest_vars <- manifest_vars[!grepl("eta", manifest_vars)]
    Values_indicators <- Values[, manifest_vars, drop = FALSE]
    Results <- cbind(Values_indicators, Y)
  }
  as.data.frame(Results)
}
  

