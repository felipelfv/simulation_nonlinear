##==============================================================================
## General Information 
##
## Script name: GenerateData.R
##
## Purpose of script: Generate population model according to specified parameter
##
#' @param model A lavaan-based syntax model. 
#' @param N Sample size for the generated dataset. 
#' @param skewness Skewness parameter for exogenous variables when using rIG. 
#' @param excesskurtosis Excess kurtosis parameter for  exogenous variables when using rIG.
#' @param exo.mean Vector of means for exogenous variables when not centering. Important: latent variables first, then manifest variables.
#' @param distr.exo Distribution for exogenous variables. Either "rIG" - uses covsim::rIG with skewness/kurtosis - or "unif" - uniform via normal transformation.
#' @param distr.zeta Distribution for structural residuals: "normal" or "exp.rate1".
#' @param distr.epsilon Distribution for measurement errors. Either "normal" or exponential.
#' @param center.exogenous.latent Logical. Whether to center exogenous latent variables.
#' @param center.exogenous.manifest Logical. Whether to mean center exogenous manifest. 
#' @param center.lv.dependent Logical. Whether to mean center dependent latent variables after generation.
#' @param center.lv.prod Logical. Whether to mean center interaction/product terms between latent variables.
#' @param center.indicators Logical. Whether to mean center indicator variables.
#' @param target.var Named list of target variances for dependent variables. Overrides R2 calculations if specified.
#' @param R2 Named list of R-squared values for dependent variables. Used to calculate residual variance when target.var not specified.
#' @param rel Reliability coefficient for indicators. Used to calculate measurement error variance.
#' @param seed Value passed into set.seed.
#' @param add.eta Logical. Whether to include latent variables in final dataset. 
#' @return A data.frame containing the generated data
#' 
##==============================================================================

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
                         rel = 0.90,
                         seed = NULL,
                         add.eta = FALSE,
                         pop_var_nozeta = NULL,
                         eta_pop_vars = NULL,
                         compute_pop_vars = FALSE,
                         N.pop = 500000L) {
  
  if(!is.null(seed)) set.seed(seed)
  
  # required population variances
  if(!compute_pop_vars && (is.null(pop_var_nozeta) || is.null(eta_pop_vars))) {
    stop("pop_var_nozeta and eta_pop_vars are required. Either provide them or set compute_pop_vars=TRUE.")
  }
  
  ##============================================================================
  ## MODEL INFORMATION
  ##============================================================================
  
  # Lavaan's parsed model information
  fit <- lavaan::sam(model)
  pt <- lavaan::parTable(fit)
  
  # Structural model information
  model_info <- list(structural = list(
    dependent = character(), 
    exogenous = character(),
    equations = list(),
    coefficients = list(),
    interactions = list(),
    generation_order = NULL
  ))
  
  # Process structural model
  model_info$structural$dependent <- unique(pt$lhs[pt$op == "~"])
  
  # Process equations, coefficients, and interactions 
  for(dv in model_info$structural$dependent) {
    eq_rows <- pt$op == "~" & pt$lhs == dv
    model_info$structural$equations[[dv]] <- pt$rhs[eq_rows]
    model_info$structural$coefficients[[dv]] <- pt$est[eq_rows]
    
    # Interactions 
    dv_interactions <- intersect(
      pt$rhs[eq_rows],
      c(fit@pta$vnames$lv.interaction[[1]], fit@pta$vnames$ov.interaction[[1]])
    )
    if (length(dv_interactions) > 0) {
      model_info$structural$interactions[[dv]] <- dv_interactions
    }
  }
  
  # All base variables
  base_vars <- c(
    fit@pta$vnames$lv.regular[[1]],
    setdiff(fit@pta$vnames$ov[[1]], fit@pta$vnames$ov.ind[[1]])
  )
  
  # "B" matrix
  B <- matrix(NA, nrow = length(base_vars), ncol = length(base_vars))
  rownames(B) <- colnames(B) <- base_vars
  
  # Process each regression equation to get dependencies
  for(dv in model_info$structural$dependent) {
    rhs_terms <- pt$rhs[pt$op == "~" & pt$lhs == dv]
    
    for(i in rhs_terms) {
      components <- unlist(strsplit(rhs_terms, ":"))
      for(comp in components) {
        if(comp %in% base_vars) {
          B[dv, comp] <- 1
        }
      }
    }
  }
  
  # Ancestors 
  ancestors <- lavaan:::lav_utils_get_ancestors(B)
  names(ancestors) <- rownames(B)
  ancestor_lengths <- sapply(ancestors, length)
  # Exogenous variables
  model_info$structural$exogenous <- names(ancestor_lengths)[ancestor_lengths == 0]
  # Generation order
  model_info$structural$generation_order <- names(ancestor_lengths)[order(ancestor_lengths)]
  
  ##============================================================================
  ## PART JUST FOR COMPUTING POPULATION VARIANCES (IF REQUESTED OR NOT PROVIDED)
  ##============================================================================
  
  if(compute_pop_vars || is.null(pop_var_nozeta) || is.null(eta_pop_vars)) {
    
    pop_var_nozeta <- list()
    eta_pop_vars <- list()
    
    # matrix for all variables
    all_vars_pop <- c(fit@pta$vnames$lv.regular[[1]], 
                      fit@pta$vnames$lv.interaction[[1]],
                      setdiff(fit@pta$vnames$ov[[1]], 
                              fit@pta$vnames$ov.ind[[1]]))
    
    Values_pop <- matrix(NA, nrow = N.pop, ncol = length(all_vars_pop))
    colnames(Values_pop) <- all_vars_pop
    
    # exo.vcov for population generation
    exo_vars_pop <- model_info$structural$exogenous
    psi_matrix_pop <- lavaan::lavInspect(fit, "est")$psi
    exo.vcov_pop <- psi_matrix_pop[exo_vars_pop, exo_vars_pop, drop = FALSE]
    
    for(i in 1:nrow(exo.vcov_pop)) { 
      if(is.na(exo.vcov_pop[i,i])) {
        exo.vcov_pop[i,i] <- 1
      }
    }
    
    # generate exogenous variables for population
    generate_exo_pop <- function() {
      if(distr.exo == "unif") {
        sd_target <- sqrt(diag(exo.vcov_pop))
        EXO <- as.matrix(faux::rmulti(
          n = N.pop,
          dist = setNames(rep("unif", length(exo_vars_pop)), exo_vars_pop),
          params = lapply(sd_target, function(sd) c(min = -sd * sqrt(3), max = sd * sqrt(3))),
          r = cov2cor(exo.vcov_pop),
          empirical = FALSE
        ))
      } else {
        EXO <- covsim::rIG(N.pop, sigma = exo.vcov_pop, skewness = skewness, 
                           excesskurtosis = excesskurtosis)[[1]]
      }
      EXO
    }
    
    EXO_pop <- generate_exo_pop()
    if(any(apply(EXO_pop, 2, stats::var) > 2 * diag(exo.vcov_pop))) {
      EXO_pop <- generate_exo_pop()
    }
    
    colnames(EXO_pop) <- exo_vars_pop
    
    # centering to population data
    eta_cols_pop <- grep("eta", exo_vars_pop, value = TRUE)
    if(center.exogenous.latent && length(eta_cols_pop) > 0) {
      EXO_pop[, eta_cols_pop] <- t(t(EXO_pop[, eta_cols_pop]) - colMeans(EXO_pop[, eta_cols_pop]))
    } else if(!is.null(exo.mean) && length(eta_cols_pop) > 0) {
      EXO_pop[, eta_cols_pop] <- t(t(EXO_pop[, eta_cols_pop]) + exo.mean[1:length(eta_cols_pop)])
    }
    
    manifest_cols_pop <- setdiff(exo_vars_pop, eta_cols_pop)
    if(center.exogenous.manifest && length(manifest_cols_pop) > 0) {
      EXO_pop[, manifest_cols_pop] <- t(t(EXO_pop[, manifest_cols_pop]) - colMeans(EXO_pop[, manifest_cols_pop]))
    } else if(!is.null(exo.mean) && length(exo.mean) > length(eta_cols_pop) && length(manifest_cols_pop) > 0) {
      EXO_pop[, manifest_cols_pop] <- t(t(EXO_pop[, manifest_cols_pop]) + 
                                          exo.mean[(length(eta_cols_pop) + 1):length(exo.mean)])
    }
    
    Values_pop[, exo_vars_pop] <- EXO_pop
    
    # generate dependent variables to compute their population variances
    intercepts_pop <- lavaan::lavInspect(fit, "est")$alpha
    for(var in model_info$structural$generation_order) {
      if(var %in% model_info$structural$dependent) {
        terms <- model_info$structural$equations[[var]]
        
        # handle interactions
        if(!is.null(model_info$structural$interactions[[var]])) {
          for(inter in model_info$structural$interactions[[var]]) {
            if(all(is.na(Values_pop[, inter]))) {
              components <- unlist(strsplit(inter, ":"))
              Values_pop[, inter] <- Values_pop[, components[1]] * Values_pop[, components[2]]
              
              if(center.lv.prod) {
                Values_pop[, inter] <- Values_pop[, inter] - mean(Values_pop[, inter])
              }
            }
          }
        }
        
        # calculate deterministic part
        equation_coefs <- model_info$structural$coefficients[[var]]
        deterministic_part <- Values_pop[, terms, drop = FALSE] %*% equation_coefs
        
        if (!is.null(intercepts_pop) && !is.na(intercepts_pop[var,1])) {
          deterministic_part <- intercepts_pop[var,1] + deterministic_part 
        }
        
        # store variance of deterministic part
        pop_var_nozeta[[var]] <- stats::var(deterministic_part) * (N.pop-1)/N.pop
        
        # generate full variable with residual for eta population variances
        if(!is.null(target.var[[var]])) {
          if(pop_var_nozeta[[var]] < target.var[[var]]) {
            current_target_var <- target.var[[var]] - pop_var_nozeta[[var]]
          } else {
            if(is.null(R2[[var]])) {
              stop(sprintf("No R2 value provided for %s. R2 must be specified.", var))
            }
            current_target_var <- pop_var_nozeta[[var]] * ((1-R2[[var]])/R2[[var]])
          }
        } else {
          if(pop_var_nozeta[[var]] < 0.0001) {
            current_target_var <- 1
          } else {
            if(is.null(R2[[var]])) {
              stop(sprintf("No R2 value provided for %s. R2 must be specified.", var))
            }
            current_target_var <- pop_var_nozeta[[var]] * ((1-R2[[var]])/R2[[var]])
          }
        }
        
        zeta_pop <- switch(distr.zeta,
                           "normal" = rnorm(N.pop, 0, sqrt(current_target_var)),
                           "exp.rate1" = rexp(N.pop, rate = 1/sqrt(current_target_var)) - sqrt(current_target_var),
                           stop(paste("wrong option for distr.zeta:", distr.zeta)))
        
        #zeta_pop <- zeta_pop - mean(zeta_pop)
        #zeta_pop <- (zeta_pop - mean(zeta_pop)) * sqrt(c(current_target_var) / 
                                                         #c(stats::var(zeta_pop - mean(zeta_pop)) * (N.pop-1)/N.pop))
        
        Values_pop[, var] <- deterministic_part + zeta_pop
        
        if(center.lv.dependent) {
          Values_pop[, var] <- t(t(Values_pop[, var]) - mean(Values_pop[, var]))
        }
      }
    }
    
    # population variances for all latent variables
    for(eta in fit@pta$vnames$lv.regular[[1]]) {
      if(eta %in% exo_vars_pop) {
        # For exogenous: use theoretical value from model
        eta_idx <- which(exo_vars_pop == eta)
        eta_pop_vars[[eta]] <- exo.vcov_pop[eta_idx, eta_idx]
      } else {
        # For dependent: compute empirically
        eta_pop_vars[[eta]] <- stats::var(Values_pop[, eta]) * (N.pop-1)/N.pop
      }
    }
  }
  
  # population variances for exogenous variables
  exo_vars <- model_info$structural$exogenous
  
  ##============================================================================
  ## DATA GENERATION (GIVEN THAT WE "KNOW") THE TRUE VARIANCES FROM ABOVE
  ##============================================================================
  
  ##============================================================================
  ## STRUCTURAL PART
  ##============================================================================
  
  # Matrix for all variables
  all_vars <- c(fit@pta$vnames$lv.regular[[1]],
                fit@pta$vnames$lv.interaction[[1]],
                setdiff(fit@pta$vnames$ov[[1]], 
                        fit@pta$vnames$ov.ind[[1]]))
  
  Values <- matrix(NA, nrow = N, ncol = length(all_vars))
  colnames(Values) <- all_vars
  
  # GENERATE ALL EXOGENOUS VARIABLES AT ONCE WITH rIG
  psi_matrix <- lavaan::lavInspect(fit, "est")$psi
  exo.vcov <- psi_matrix[exo_vars, exo_vars, drop = FALSE]
  
  for(i in 1:nrow(exo.vcov)) { 
    if(is.na(exo.vcov[i,i])) {
      exo.vcov[i,i] <- 1
    }
  }
  
  # function to generate exogenous variables
  generate_exo <- function() {
    if(distr.exo == "unif") {
      sd_target <- sqrt(diag(exo.vcov))
      EXO <- as.matrix(faux::rmulti(
        n = N,
        dist = setNames(rep("unif", length(exo_vars)), exo_vars),
        params = lapply(sd_target, function(sd) c(min = -sd * sqrt(3), max = sd * sqrt(3))),
        r = cov2cor(exo.vcov),
        empirical = FALSE # important to be false
      ))
    } else {
      EXO <- covsim::rIG(N, sigma = exo.vcov, skewness = skewness, 
                         excesskurtosis = excesskurtosis)[[1]]
    }
    EXO
  }
  
  # generate initial EXO
  EXO <- generate_exo()
  
  # check variances: if they are too large (>2times the intended value)
  if(any(apply(EXO, 2, stats::var) > 2 * diag(exo.vcov))) {
    var_ratios <- apply(EXO, 2, stats::var) / diag(exo.vcov)
    problem_vars <- exo_vars[var_ratios > 2]
    EXO <- generate_exo()
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
  
  # store zetas for potential covariance checking
  stored_zetas <- list()
  
  # GENERATE DEPENDENT VARIABLES ACCORDING TO THE DEPENDENCIES 
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
      
      # check coefficients
      if(!is.null(model_info$structural$coefficients[[var]])) { 
        equation_coefs <- model_info$structural$coefficients[[var]]
        if(length(equation_coefs) != length(terms)) {
          stop(sprintf("Error for %s: Number of coefficients (%d) does not match number of terms (%d)", 
                       var, length(equation_coefs), length(terms)))
        }
      } else {
        stop(sprintf("No coefficients provided for %s. Coefficients must be specified.", var))
      }
      
      # check if all terms exist in the matrix
      missing_terms <- terms[!terms %in% colnames(Values)]
      if(length(missing_terms) > 0) {
        stop(sprintf("Terms not found in data for variable '%s': %s", 
                     var, paste(missing_terms, collapse = ", ")))
      }
      
      # matrix algebra for the deterministic part
      deterministic_part <- Values[, terms, drop = FALSE] %*% equation_coefs
      
      if (!is.null(intercepts) && !is.na(intercepts[var,1])) {
        deterministic_part <- intercepts[var,1] + deterministic_part
      }
      
      # Use population variance (!!!)
      var.nozeta <- pop_var_nozeta[[var]]
      
      # calculate target variance
      if(!is.null(target.var) && !is.null(target.var[[var]])) {
        if(var.nozeta < target.var[[var]]) {
          current_target_var <- target.var[[var]] - var.nozeta
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
        if(var.nozeta < 0.0001) {
          current_target_var <- 1
        } else {
          if(is.null(R2) || is.null(R2[[var]])) {
            stop(sprintf("No R2 value provided for %s. R2 must be specified.", var))
          }
          current_target_var <- var.nozeta * ((1-R2[[var]])/R2[[var]])
        }
      }
      
      # generate zeta 
      zeta <- switch(distr.zeta,
                     "normal" = rnorm(N, 0, sqrt(current_target_var)),
                     "exp.rate1" = rexp(N, rate = 1/sqrt(current_target_var)) - sqrt(current_target_var),
                     stop(paste("wrong option for distr.zeta:", distr.zeta)))
      
      # not standardization anymore
      #zeta <- zeta - mean(zeta)
      #zeta <- (zeta - mean(zeta)) * sqrt(c(current_target_var) / c(stats::var(zeta - mean(zeta)) * (N-1)/N))
      
      # store zeta for checking covariances
      stored_zetas[[var]] <- zeta
      
      # entries for the variables with the residual
      Values[, var] <- deterministic_part + zeta
      
      if(center.lv.dependent) {
        Values[, var] <- t(t(Values[, var]) - mean(Values[, var]))
      }
    }
  }
  
  ##============================================================================
  ## MEASUREMENT PART
  ##============================================================================
  
  lambda <- lavaan::lavInspect(fit, "est")$lambda
  LAMBDA <- lambda[startsWith(rownames(lambda), "x"), fit@pta$vnames$lv.regular[[1]], 
                   drop = FALSE]
  
  # Use population variances
  eta_vars <- sapply(fit@pta$vnames$lv.regular[[1]], 
                     function(eta) eta_pop_vars[[eta]])
  
  target_var_indicators <- eta_vars * (1/rel - 1)
  
  # generate errors in theta 
  n_indicators <- nrow(LAMBDA)
  indicator_groups <- lapply(seq_along(fit@pta$vnames$lv.regular[[1]]),
                             function(i) which(LAMBDA[, i] == 1))
  
  THETA <- matrix(NA, nrow = N, ncol = n_indicators)
  colnames(THETA) <- rownames(LAMBDA)
  
  # THETA matrix by the indicators
  for(i in seq_along(indicator_groups)) {
    inds <- indicator_groups[[i]]
    t_var_ind <- target_var_indicators[i]
    
    if(distr.epsilon == "normal") {
      THETA[, inds] <- matrix(rnorm(N * length(inds), 0, sqrt(t_var_ind)), N, length(inds))
    } else {
      THETA[, inds] <- matrix(rexp(N * length(inds), rate = 1/sqrt(t_var_ind)) - sqrt(t_var_ind), N, length(inds))
    }
  }
  
  # to get the actual values for indicators: Y = eta %*% t(LAMBDA) + THETA
  Y <- as.matrix(Values[, fit@pta$vnames$lv.regular[[1]]]) %*% t(LAMBDA) + THETA 
  colnames(Y) <- rownames(LAMBDA)
  
  # eenter indicators 
  if(center.indicators) {
    Y <- t(t(Y) - colMeans(Y))
  }
  
  # keep etas or not in the final dataset
  if(add.eta) {
    Results <- cbind(Values, Y)
  } else {
    # get manifest variables from the model excluding indicators
    manifest_vars <- setdiff(fit@pta$vnames$ov[[1]], fit@pta$vnames$ov.ind[[1]])
    # remove interaction terms and variables latent variables
    manifest_vars <- manifest_vars[!grepl("eta", manifest_vars)]
    Values_indicators <- Values[, manifest_vars, drop = FALSE]
    Results <- cbind(Values_indicators, Y)
  }
  
  result <- as.data.frame(Results)
  
  # when population variances were computed, attach to result
  if(compute_pop_vars) {
    attr(result, "pop_var_nozeta") <- pop_var_nozeta
    attr(result, "eta_pop_vars") <- eta_pop_vars
    attr(result, "exo.vcov") <- if(exists("exo.vcov_pop")) exo.vcov_pop else exo.vcov
    attr(result, "model_info") <- model_info
  }
  
  result
}