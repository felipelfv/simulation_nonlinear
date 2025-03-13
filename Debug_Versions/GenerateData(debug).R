#### 1. General Information ####

# This file contains the function for generating the datasets based on the 
# lavaan-syntax based models. It is literally the same as the original function
# used in the simulations. The only difference is that here I have some debug 
# and checks throughout the lines to ensure that I did things correctly as the 
# function got incrementally more complex (personal statement) and added much 
# more information that actually needed for understanding the function. This 
# amout of visual information could hinder understanding from others. Hence,
# this is the reason for creating and adding this to a separate script. 


#### 2. Debug version for GenerateData() ####

GenerateData <- function(model, 
                         N = 1000L,
                         skewness = NULL,
                         excesskurtosis = NULL,
                         exo.mean = NULL,
                         distr.exo = "rIG",
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
                         add.eta = FALSE,
                         verbose = FALSE) {  # Added verbose parameter
  
  if(!is.null(seed)) set.seed(seed)
  
  # Print debug function for tracking issues
  print_debug <- function(...) {
    if(verbose) {
      cat(..., "\n")
    }
  }
  
  print_debug("\n=== Starting Data Generation Process ===")
  
  ###################################
  ######## Model Information ########
  ###################################
  
  print_debug("\n=== Step 1: Parsing Model Information ===")
  
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
  print_debug("Dependent variables identified:")
  print_debug(paste("- ", model_info$structural$dependent, collapse = "\n"))
  
  # Process equations, coefficients, and interactions 
  for(dv in model_info$structural$dependent) { # For each dependent variable
    print_debug(sprintf("\nProcessing equations for: %s", dv))
    eq_rows <- pt$op == "~" & pt$lhs == dv # Select the predictors
    # Store equations and coefficients for each dependent 
    model_info$structural$equations[[dv]] <- pt$rhs[eq_rows]
    model_info$structural$coefficients[[dv]] <- pt$est[eq_rows] # Coefficients associated with each predictor
    
    print_debug("Predictors:")
    for(i in seq_along(model_info$structural$equations[[dv]])) {
      print_debug(sprintf("- %s: %.4f", 
                          model_info$structural$equations[[dv]][i], 
                          model_info$structural$coefficients[[dv]][i]))
    }
    
    # Interactions 
    dv_interactions <- intersect( # Rows in both x and y 
      pt$rhs[eq_rows],
      c(fit@pta$vnames$lv.interaction[[1]], fit@pta$vnames$ov.interaction[[1]])
    )
    if (length(dv_interactions) > 0) { # Just to keep the output clean from the parsed model
      model_info$structural$interactions[[dv]] <- dv_interactions # So, if no interaction, nothing instead of character(0) showing
      print_debug("Interactions:")
      print_debug(paste("- ", dv_interactions, collapse = "\n"))
    }
  }
  
  # All base variables (main and involved in interactions) - this is relevant for the "B" matrix(!)
  # Reasoning: there could be interaction among observed and latent, there could be interaction among latents that are not single predictors for any dependent
  base_vars <- c(
    fit@pta$vnames$lv.regular[[1]], # Latent variables
    setdiff(fit@pta$vnames$ov[[1]], fit@pta$vnames$ov.ind[[1]]) # Observed variables
  )
  
  print_debug("\nBase variables identified:")
  print_debug(paste("- ", base_vars, collapse = "\n"))
  
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
  
  print_debug("\nDependency matrix B created")
  if(verbose) print(B)
  
  # Ancestors 
  ancestors <- lavaan:::lav_utils_get_ancestors(B) # List for each variable in the "B" matrix
  names(ancestors) <- rownames(B)
  ancestor_lengths <- sapply(ancestors, length) # This is the crucial part that works perfectly with lav_utils_get_ancestors 
  # Exogenous variables
  model_info$structural$exogenous <- names(ancestor_lengths)[ancestor_lengths == 0] # 0 ancestors
  # Generation order
  model_info$structural$generation_order <- names(ancestor_lengths)[order(ancestor_lengths)] # based on ancestors order given by the length
  
  print_debug("\nExogenous variables identified:")
  print_debug(paste("- ", model_info$structural$exogenous, collapse = "\n"))
  
  print_debug("\nGeneration order determined:")
  print_debug(paste("- ", model_info$structural$generation_order, collapse = "\n"))
  
  #################################
  ######## Data Generation ########
  #################################
  
  print_debug("\n=== Step 2: Preparing Data Generation ===")
  
  # Matrix for all variables
  # Maybe just $eqs.y[[1]] and $eqs.x[[1]] is enough?
  # Also add indicators already?
  all_vars <- c(fit@pta$vnames$lv.regular[[1]], # Regular latent variables
                fit@pta$vnames$lv.interaction[[1]], # Interaction terms
                setdiff(fit@pta$vnames$ov[[1]], 
                        fit@pta$vnames$ov.ind[[1]])) # Observed variables
  
  print_debug("All variables to be generated:")
  print_debug(paste("- ", all_vars, collapse = "\n"))
  
  # Also, matrix instead of data frame: to be confirmed with Yves
  Values <- matrix(NA, nrow = N, ncol = length(all_vars))
  colnames(Values) <- all_vars
  # Could also just skip the matrix initialization here actually..
  
  print_debug("\n=== Step 3: Generating Exogenous Variables ===")
  
  # GENERATE ALL EXOGENOUS VARIABLES AT ONCE WITH rIG
  # Important: Add check related to variance as Yves has (comment in 06/12/24)
  exo_vars <- model_info$structural$exogenous
  psi_matrix <- lavaan::lavInspect(fit, "est")$psi # For (co)variances among exogenous
  # Extract the exogenous variables in PSI
  exo.vcov <- psi_matrix[exo_vars, exo_vars, drop = FALSE]
  
  print_debug("Exogenous variables covariance matrix:")
  if(verbose) print(exo.vcov)
  
  for(i in 1:nrow(exo.vcov)) { 
    if(is.na(exo.vcov[i,i])) { # If variance left unspecified
      print_debug(sprintf("Setting unspecified variance for %s to 1", rownames(exo.vcov)[i]))
      exo.vcov[i,i] <- 1 # We set to 1. Otherwise, the value given in the model is retained
    }
  }
  
  print_debug(sprintf("\nGenerating exogenous variables using '%s' distribution", distr.exo))
  
  if(distr.exo == "unif") {
    # "NAIVE" approach; ignoring copulas
    # For uniform distribution as in GAPI article
    print_debug("Generating with uniform distribution")
    Z <- MASS::mvrnorm(N, mu = rep(0, ncol(exo.vcov)), Sigma = exo.vcov)
    EXO <- pnorm(Z) # Transform into an uniform distribution
    sd_desired <- sqrt(diag(exo.vcov)) # sd (from diag of cov matrix)
    scaling <- sqrt(12) * sd_desired
    EXO <- EXO - 0.5 # Center all variables by subtracting 0.5
    # Scale by multiplying each column by its scaling factor:
    EXO <- sweep(EXO, MARGIN = 2, STATS = scaling, FUN = "*")
  } else {
    print_debug("Generating with rIG")
    EXO <- covsim::rIG(N, sigma = exo.vcov, skewness = skewness, 
                       excesskurtosis = excesskurtosis)[[1]] # Correlations are as given now!
  }
  
  colnames(EXO) <- exo_vars
  
  print_debug("\n=== Step 4: Centering Exogenous Variables ===")
  
  # Would we ever center one and not the other? I believe rarely but I have both options 
  eta_cols <- grep("eta", exo_vars, value = TRUE)
  if(center.exogenous.latent) {
    print_debug("Centering latent exogenous variables:")
    for(col in eta_cols) {
      print_debug(sprintf("\nProcessing %s:", col))
      first_val_before <- EXO[1, col]
      current_mean <- mean(EXO[, col])
      EXO[, col] <- EXO[, col] - current_mean
      first_val_after <- EXO[1, col]
      
      print_debug(sprintf("- First value before centering: %.6f", first_val_before))
      print_debug(sprintf("- Mean used for centering: %.6f", current_mean))
      print_debug(sprintf("- First value after centering: %.6f", first_val_after))
      print_debug(sprintf("- Verification: %.6f - %.6f = %.6f", 
                          first_val_before, current_mean, first_val_after))
      print_debug(sprintf("- New mean (should be ~0): %.10f", mean(EXO[, col])))
    }
  } else if(!is.null(exo.mean)) {
    print_debug("Adding means to latent exogenous variables:")
    for(i in seq_along(eta_cols)) {
      print_debug(sprintf("\nProcessing %s:", eta_cols[i]))
      first_val_before <- EXO[1, eta_cols[i]]
      mean_val <- exo.mean[i]
      EXO[, eta_cols[i]] <- EXO[, eta_cols[i]] + mean_val
      first_val_after <- EXO[1, eta_cols[i]]
      
      print_debug(sprintf("- First value before adding mean: %.6f", first_val_before))
      print_debug(sprintf("- Mean value added: %.6f", mean_val))
      print_debug(sprintf("- First value after adding mean: %.6f", first_val_after))
      print_debug(sprintf("- Verification: %.6f + %.6f = %.6f", 
                          first_val_before, mean_val, first_val_after))
    }
  }
  
  manifest_cols <- setdiff(exo_vars, eta_cols) # All rows in x not in y
  if(center.exogenous.manifest) {
    print_debug("Centering manifest exogenous variables:")
    for(col in manifest_cols) {
      print_debug(sprintf("\nProcessing %s:", col))
      first_val_before <- EXO[1, col]
      current_mean <- mean(EXO[, col])
      EXO[, col] <- EXO[, col] - current_mean
      first_val_after <- EXO[1, col]
      
      print_debug(sprintf("- First value before centering: %.6f", first_val_before))
      print_debug(sprintf("- Mean used for centering: %.6f", current_mean))
      print_debug(sprintf("- First value after centering: %.6f", first_val_after))
      print_debug(sprintf("- Verification: %.6f - %.6f = %.6f", 
                          first_val_before, current_mean, first_val_after))
      print_debug(sprintf("- New mean (should be ~0): %.10f", mean(EXO[, col])))
    }
  } else if(!is.null(exo.mean)) {
    print_debug("Adding means to manifest exogenous variables:")
    offset = length(eta_cols)
    for(i in seq_along(manifest_cols)) {
      print_debug(sprintf("\nProcessing %s:", manifest_cols[i]))
      first_val_before <- EXO[1, manifest_cols[i]]
      mean_val <- exo.mean[i + offset]
      EXO[, manifest_cols[i]] <- EXO[, manifest_cols[i]] + mean_val
      first_val_after <- EXO[1, manifest_cols[i]]
      
      print_debug(sprintf("- First value before adding mean: %.6f", first_val_before))
      print_debug(sprintf("- Mean value added: %.6f", mean_val))
      print_debug(sprintf("- First value after adding mean: %.6f", first_val_after))
      print_debug(sprintf("- Verification: %.6f + %.6f = %.6f", 
                          first_val_before, mean_val, first_val_after))
    }
  }
  
  Values[, exo_vars] <- EXO
  
  print_debug("\n=== Step 5: Generating Dependent Variables ===")
  
  # Store zetas for potential covariance checking
  stored_zetas <- list()
  
  # GENERATE DEPENDENT VARIABLES ACCORDING TO THE DEPENDENCIES 
  intercepts <- lavaan::lavInspect(fit, "est")$alpha # Intercepts for the dependent variables equations
  # Remember [[var]] read as string [["var"]]
  for(var in model_info$structural$generation_order) { 
    if(var %in% model_info$structural$dependent) { # Checks if current variable is a dependent variable
      print_debug(sprintf("\nProcessing dependent variable: %s", var))
      
      terms <- model_info$structural$equations[[var]] # Equation terms for current variable
      # Handle interactions first
      if(!is.null(model_info$structural$interactions[[var]])) { # Checks if variable has interaction terms
        print_debug("Creating interaction terms:")
        for(inter in model_info$structural$interactions[[var]]) {
          if(all(is.na(Values[, inter]))) {  # Checks if interaction term values are missing in the dataset
            components <- unlist(strsplit(inter, ":")) # Split the interaction into their individual variables
            print_debug(sprintf("- Creating %s = %s * %s", 
                                inter, components[1], components[2]))
            
            # Multiplication for interactions
            # Check this if correctly applied; maybe should be the transpose * X? I dont think so; hence, just element-wise.
            Values[, inter] <- Values[, components[1]] * Values[, components[2]]
            
            # Check this:
            if(center.lv.prod) {
              print_debug(sprintf("Centering interaction term: %s", inter))
              first_val_before <- Values[1, inter]
              current_mean <- mean(Values[, inter])
              Values[, inter] <- Values[, inter] - mean(Values[, inter]) # colMeans wont work with 1 column
              first_val_after <- Values[1, inter]
              
              print_debug(sprintf("- First value before centering: %.6f", first_val_before))
              print_debug(sprintf("- Mean used for centering: %.6f", current_mean))
              print_debug(sprintf("- First value after centering: %.6f", first_val_after))
              print_debug(sprintf("- Verification: %.6f - %.6f = %.6f", 
                                  first_val_before, current_mean, first_val_after))
              print_debug(sprintf("- New mean (should be ~0): %.10f", mean(Values[, inter])))
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
        
        print_debug("Using coefficients:")
        for(i in seq_along(terms)) {
          print_debug(sprintf("- %s: %.4f", terms[i], equation_coefs[i]))
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
      print_debug("Calculating deterministic part")
      deterministic_part <- Values[, terms, drop = FALSE] %*% equation_coefs
      if (!is.null(intercepts) && !is.na(intercepts[var,1])) {
        print_debug(sprintf("Adding intercept: %.6f", intercepts[var,1]))
        deterministic_part <- intercepts[var,1] + deterministic_part 
      }
      
      # Calculate variance without the residual/error
      var.nozeta <- var(deterministic_part)
      print_debug(sprintf("Deterministic part variance (var.nozeta): %.6f", var.nozeta))
      
      # Calculate target variance (this is given by the user)
      if(!is.null(target.var) && !is.null(target.var[[var]])) {
        print_debug(sprintf("Specified target variance for %s: %.6f", 
                            var, target.var[[var]]))
        
        if(var.nozeta < target.var[[var]]) {
          current_target_var <- target.var[[var]] - var.nozeta
          print_debug("var.nozeta is smaller than target.var - using difference for error variance")
        } else {
          warning(sprintf(
            "var.nozeta [%f] is larger than target.var [%f] for %s -- using R2 instead",
            var.nozeta, target.var[[var]], var
          ))
          if(is.null(R2) || is.null(R2[[var]])) {
            stop(sprintf("No R2 value provided for %s. R2 must be specified.", var))
          }
          current_target_var <- var.nozeta * ((1-R2[[var]])/R2[[var]])
        }
      } else {
        print_debug("No target variance specified - using R² to determine error variance")
        if(var.nozeta < 0.0001) {
          current_target_var <- 1
          print_debug("Near-zero deterministic variance, setting target.var = 1")
        } else {
          if(is.null(R2) || is.null(R2[[var]])) {
            stop(sprintf("No R2 value provided for %s. R2 must be specified.", var))
          }
          current_target_var <- var.nozeta * ((1-R2[[var]])/R2[[var]])
        }
      }
      print_debug(sprintf("Target variance for error: %.6f", current_target_var))
      
      # Generate and standardize zeta 
      print_debug(sprintf("Generating %s error term", distr.zeta))
      zeta <- switch(distr.zeta,
                     "normal" = rnorm(N, 0, sqrt(current_target_var)),
                     "exp.rate1" = rexp(N, rate = 1/sqrt(current_target_var)) - 1,
                     stop(paste("wrong option for distr.zeta:", distr.zeta)))
      
      zeta <- zeta - mean(zeta)
      zeta <- (zeta - mean(zeta)) * sqrt(c(current_target_var) / c(var(zeta - mean(zeta))))
      
      print_debug("Error term diagnostics:")
      print_debug(sprintf("- Mean: %.6f", mean(zeta)))
      print_debug(sprintf("- Variance: %.6f", var(zeta)))
      
      # Store zeta for checking covariances
      stored_zetas[[var]] <- zeta
      
      # Entries for the variables with the residual
      Values[, var] <- deterministic_part + zeta # Add on top of the deterministic part
      
      if(center.lv.dependent) {
        print_debug(sprintf("Centering dependent variable: %s", var))
        first_val_before <- Values[1, var]
        current_mean <- mean(Values[, var])
        Values[, var] <- t(t(Values[, var]) - mean(Values[, var])) # Same remark as with lv.prod:
        first_val_after <- Values[1, var]
        
        print_debug(sprintf("- First value before centering: %.6f", first_val_before))
        print_debug(sprintf("- Mean used for centering: %.6f", current_mean))
        print_debug(sprintf("- First value after centering: %.6f", first_val_after))
        print_debug(sprintf("- Verification: %.6f - %.6f = %.6f", 
                            first_val_before, current_mean, first_val_after))
        print_debug(sprintf("- New mean (should be ~0): %.10f", mean(Values[, var])))
      }
      
      print_debug(sprintf("Final diagnostics for %s:", var))
      print_debug(sprintf("- Mean: %.6f", mean(Values[, var])))
      print_debug(sprintf("- Total variance: %.6f", var(Values[, var])))
      print_debug(sprintf("- Achieved R²: %.6f", 
                          cor(deterministic_part, Values[, var])^2))
    }
  }
  
  # Check zeta covariances
  if(verbose && length(stored_zetas) > 1) {
    zeta_matrix <- do.call(cbind, stored_zetas)
    colnames(zeta_matrix) <- names(stored_zetas)
    
    print_debug("\nZeta correlations:")
    print_debug(cor(zeta_matrix))
    print_debug("\nZeta covariances:")
    print_debug(cov(zeta_matrix))
  }
  
  print_debug("\n=== Step 6: Generating Measurement Model ===")
  
  # GENERATE THE MEASUREMENT PART
  # Get LAMBDA and calculate values for each pure eta
  lambda <- lavaan::lavInspect(fit, "est")$lambda
  #pure_etas <- fit@pta$vnames$lv.regular[[1]]
  LAMBDA <- lambda[startsWith(rownames(lambda), "x"), fit@pta$vnames$lv.regular[[1]], 
                   drop = FALSE] # Limitation: assumes indicators as x only
  # The above to drop the manifest variables from rows
  
  print_debug("Lambda matrix:")
  if(verbose) print(LAMBDA)
  
  # Calculates target variances for the indicators
  # Could also use values from above for loop but for now just making sure things work as expected
  # We first get the variance of the latent variables
  print_debug("Calculating target variances for indicators")
  eta_vars <- sapply(fit@pta$vnames$lv.regular[[1]], 
                     function(eta) var(Values[, eta])) # Correct to use N-1 as default in `var()`
  target_var_indicators <- eta_vars * (1/rel - 1)
  
  print_debug("Indicator target variance calculations:")
  for(i in seq_along(fit@pta$vnames$lv.regular[[1]])) {
    eta <- fit@pta$vnames$lv.regular[[1]][i]
    print_debug(sprintf("- For %s: var(eta) = %.6f, (1/rel - 1) = %.6f, target_var = %.6f", 
                        eta, eta_vars[i], (1/rel - 1), target_var_indicators[i]))
  }
  
  # Generate errors in theta 
  n_indicators <- nrow(LAMBDA)
  # For the different factor models below: 
  indicator_groups <- lapply(seq_along(fit@pta$vnames$lv.regular[[1]]),
                             function(i) which(LAMBDA[, i] == 1)) # Position of elements (logical vector)
  
  THETA <- matrix(NA, nrow = N, ncol = n_indicators)
  colnames(THETA) <- rownames(LAMBDA)
  
  print_debug("Generating indicator errors")
  
  # Get THETA matrix by the indicators
  # This part is the only part I am not 100% if mathematically correct (in accordance to SEM)
  # Reasoning: All indicators measuring the same latent variable get errors drawn from the same distribution (with the same variance)
  for(i in seq_along(indicator_groups)) { # For each latent variable,
    inds <- indicator_groups[[i]] # We get the indicators for it
    t_var_ind <- target_var_indicators[i] # We get the target variance from above calculation
    
    print_debug(sprintf("\nGenerating errors for indicators of %s:", 
                        fit@pta$vnames$lv.regular[[1]][i]))
    print_debug(sprintf("- Target error variance: %.6f", t_var_ind))
    print_debug(sprintf("- Indicators: %s", paste(rownames(LAMBDA)[inds], collapse=", ")))
    
    if(distr.epsilon == "normal") {
      print_debug("- Using normal distribution for errors")
      THETA[, inds] <- matrix(rnorm(N * length(inds), 0, sqrt(t_var_ind)), N, length(inds))
    } else {
      print_debug("- Using exponential distribution for errors")
      THETA[, inds] <- matrix(rexp(N * length(inds), rate = 1/sqrt(t_var_ind)) - sqrt(t_var_ind), N, length(inds))
    }
    
    # Check error properties for the first indicator in the group
    if(length(inds) > 0) {
      first_ind <- inds[1]
      print_debug("Error term diagnostics (first indicator):")
      print_debug(sprintf("- Mean: %.6f", mean(THETA[, first_ind])))
      print_debug(sprintf("- Variance: %.6f", var(THETA[, first_ind])))
    }
  }
  
  # To get the actual values for indicators: Y = eta %*% t(LAMBDA) + THETA
  print_debug("\nCalculating indicator values: Y = eta %*% t(LAMBDA) + THETA")
  Y <- as.matrix(Values[, fit@pta$vnames$lv.regular[[1]]]) %*% t(LAMBDA) + THETA 
  colnames(Y) <- rownames(LAMBDA)
  
  # Center indicators 
  if(center.indicators) {
    print_debug("Centering indicators")
    for(ind in colnames(Y)) {
      print_debug(sprintf("\nCentering indicator: %s", ind))
      first_val_before <- Y[1, ind]
      current_mean <- mean(Y[, ind])
      Y[, ind] <- Y[, ind] - current_mean
      first_val_after <- Y[1, ind]
      
      print_debug(sprintf("- First value before centering: %.6f", first_val_before))
      print_debug(sprintf("- Mean used for centering: %.6f", current_mean))
      print_debug(sprintf("- First value after centering: %.6f", first_val_after))
      print_debug(sprintf("- Verification: %.6f - %.6f = %.6f", 
                          first_val_before, current_mean, first_val_after))
      print_debug(sprintf("- New mean (should be ~0): %.10f", mean(Y[, ind])))
    }
  }
  
  # Check reliabilities
  print_debug("\nVerifying indicator reliabilities:")
  for(i in seq_along(fit@pta$vnames$lv.regular[[1]])) {
    eta <- fit@pta$vnames$lv.regular[[1]][i]
    inds <- rownames(LAMBDA)[LAMBDA[, i] == 1]
    
    for(ind in inds) {
      cor_with_lv <- cor(Y[, ind], Values[, eta])
      achieved_rel <- cor_with_lv^2
      print_debug(sprintf("- %s with %s: target=%.4f, achieved=%.4f, diff=%.4f", 
                          ind, eta, rel, achieved_rel, achieved_rel - rel))
    }
  }
  
  print_debug("\n=== Step 7: Finalizing Dataset ===")
  
  # Keep etas or not in the final dataset; relevant for fitting model in lavaan
  if(add.eta) {
    print_debug("Including latent variables in final dataset")
    Results <- cbind(Values, Y)
  } else {
    print_debug("Excluding latent variables from final dataset")
    # Get manifest variables from the model excluding indicators
    manifest_vars <- setdiff(fit@pta$vnames$ov[[1]], fit@pta$vnames$ov.ind[[1]])
    # Remove interaction terms and variables latent variables (as in eta..)
    manifest_vars <- manifest_vars[!grepl("eta", manifest_vars)]
    Values_indicators <- Values[, manifest_vars, drop = FALSE]
    Results <- cbind(Values_indicators, Y)
  }
  
  print_debug("\n=== Data Generation Complete ===")
  
  as.data.frame(Results)
}

#### 3. Generate Data ####

# Source required functions
source("Models.R")
source("Design.R") # Check here for most of the parameters 

# Remaining that were used in the simulation loop
skewness <- rep(0, 2)
excesskurtosis <- rep(0, 2)
# If nonnormal:
#skewness <- rep(2, 2)
#excesskurtosis <- rep(7, 2)

# For a single run:
conditions <- conditions[80,] # Randomly chosen condition
N <- 500
Rel <- 0.6
exo.mean <- rep(0, 2)
target.var <- list("eta3" = 1.0) # target variance for eta 
R2 <- list("eta3" = 0.20)
  
Data <- GenerateData(
  model = population.interaction.model,
  N = N,
  skewness = skewness,
  excesskurtosis = excesskurtosis,
  exo.mean = exo.mean,
  distr.exo = "unif",
  distr.zeta = "normal",
  distr.epsilon = "normal",
  rel = Rel,
  target.var = target.var,
  R2 = R2,
  add.eta = FALSE, 
  verbose = TRUE)

