############################ 1. General Information ############################
# GenerateData function - VITA only version
# Supports: normal, non-normal (via distribution types), and uniform

#' @param model A lavaan-based model syntax string specifying the structural equation model.
#' @param N Integer. Sample size for the generated dataset. Default is 1000L.
#' @param distr.exo Character. Distribution for exogenous variables. 
#'   Options: "normal" (default), "nonnormal", "uniform"
#' @param nonnormal.shape Numeric vector. Shape parameters for non-normal distributions.
#'   For gamma distribution: higher shape = less skewed. Default is c(1, 1, ...).
#'   Length must match number of exogenous variables.
#' @param nonnormal.rate Numeric vector. Rate parameters for non-normal distributions.
#'   For gamma distribution: controls scale. Default is c(1, 1, ...).
#'   Length must match number of exogenous variables.
#' @param exo.mean Numeric vector. Mean values for exogenous variables. Default is NULL.
#' @param distr.zeta Character. Distribution for structural residuals. Options: "normal" (default), "exp.rate1"
#' @param distr.epsilon Character. Distribution for measurement errors. Options: "normal" (default), "exp.rate1"
#' @param center.exogenous.latent Logical. Whether to center exogenous latent variables. Default is FALSE.
#' @param center.exogenous.manifest Logical. Whether to center exogenous manifest variables. Default is FALSE.
#' @param center.lv.dependent Logical. Whether to center dependent latent variables. Default is FALSE.
#' @param center.lv.prod Logical. Whether to center latent variable products/interactions. Default is FALSE.
#' @param center.indicators Logical. Whether to center indicator variables. Default is FALSE.
#' @param seed Integer. Random seed for reproducibility. Default is NULL.
#' @param add.eta Logical. Whether to include latent variables (eta) in the output dataset. Default is FALSE.
#' @param return.info Logical. Whether to return additional model information as attributes. Default is TRUE.

# Packages needed for this script:
# library(lavaan); library(covsim); library(rvinecopulib)

############################### 2. Function ####################################

GenerateData <- function(model, 
                         N = 1000L,
                         distr.exo = "normal",
                         nonnormal.shape = NULL,
                         nonnormal.rate = NULL,
                         exo.mean = NULL,
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
  
  # Check for required packages
  if(!requireNamespace("covsim", quietly = TRUE)) {
    stop("Package 'covsim' is required. Please install it.")
  }
  if(!requireNamespace("rvinecopulib", quietly = TRUE)) {
    stop("Package 'rvinecopulib' is required. Please install it.")
  }
  
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
  n_exo <- length(exo_vars)
  
  psi_matrix <- lavaan::lavInspect(fit, "est")$psi
  exo.vcov <- psi_matrix[exo_vars, exo_vars, drop = FALSE]
  
  for(i in 1:nrow(exo.vcov)) { 
    if(is.na(exo.vcov[i,i])) {
      exo.vcov[i,i] <- 1
    }
  }
  
  # helper to calculate variance from marginal
  calc_marginal_var <- function(margin) {
    distr <- margin$distr
    
    if(distr == "norm") {
      return(margin$sd^2)
    } else if(distr == "gamma") {
      return(margin$shape / margin$rate^2)
    } else if(distr == "unif") {
      return((margin$max - margin$min)^2 / 12)
    } else {
      stop(paste("Variance formula not implemented for distribution:", distr))
    }
  }
  
  # generate exogenous variables using VITA
  generate_exo <- function() {
    
    if(distr.exo == "normal") {
      # normal distributions with variance matching model
      margins_list <- lapply(1:n_exo, function(i) {
        list(
          distr = "norm",
          mean = 0,  # adjusted later if exo.mean provided
          sd = sqrt(exo.vcov[i, i])
        )
      })
      
    } else if(distr.exo == "nonnormal") {
      # gamma distributions (skewed, heavy-tailed)
      # default: shape=1, rate=1 gives skewness=2, excess kurtosis=6
      
      if(is.null(nonnormal.shape)) {
        nonnormal.shape <- rep(1, n_exo)
      }
      if(is.null(nonnormal.rate)) {
        nonnormal.rate <- rep(1, n_exo)
      }
      
      if(length(nonnormal.shape) != n_exo) {
        stop(paste("nonnormal.shape must have length", n_exo))
      }
      if(length(nonnormal.rate) != n_exo) {
        stop(paste("nonnormal.rate must have length", n_exo))
      }
      
      margins_list <- lapply(1:n_exo, function(i) {
        list(
          distr = "gamma",
          shape = nonnormal.shape[i],
          rate = nonnormal.rate[i]
        )
      })
      
    } else if(distr.exo == "uniform") {
      # uniform distributions with variance matching model
      # uniform(a, b) variance = (b-a)^2/12
      # for variance v: b - a = sqrt(12*v)
      # using symmetric bounds: [-sqrt(3*v), sqrt(3*v)]
      
      margins_list <- lapply(1:n_exo, function(i) {
        var_target <- exo.vcov[i, i]
        half_range <- sqrt(3 * var_target)
        list(
          distr = "unif",
          min = -half_range,
          max = half_range
        )
      })
      
    } else {
      stop(paste("Unknown distr.exo option:", distr.exo))
    }
    
    # calculate variances from marginals
    marginal_vars <- sapply(margins_list, calc_marginal_var)
    
    # get correlation structure from model covariance
    D_diag <- diag(exo.vcov)
    cor_matrix <- exo.vcov / sqrt(outer(D_diag, D_diag))
    
    # build covariance matrix with marginal variances
    D <- diag(sqrt(marginal_vars))
    adjusted_cov <- D %*% cor_matrix %*% D
    
    # create VITA distribution
    vitadist <- covsim::vita(margins_list, adjusted_cov, verbose = FALSE, Nmax = 10^6)
    
    # generate data
    EXO <- rvinecopulib::rvine(n = N, vine = vitadist)
    
    # SHIFTING HERE FOR RIGHT-SKEWED
    # population-level shift for nonnormal distributions (like CopSEM method)
    if(distr.exo == "nonnormal") {
      for(i in 1:ncol(EXO)) {
        natural_mean <- nonnormal.shape[i] / nonnormal.rate[i]
        EXO[, i] <- EXO[, i] - natural_mean
      }
    }
    # normal and uniform already have mean = 0, no shift needed
  
    EXO
  }
  
  # generate exogenous variables
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
      var_det <- var(as.vector(deterministic_part)) * (N-1)/N
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
      total_var <- var(Values[, var]) * (N-1)/N
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
      eta_var <- var(Values[, eta]) * (N-1)/N
      
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
