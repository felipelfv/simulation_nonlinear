############################ 1. General Information ############################
# Debug version of GenerateData function - VITA only version
# This version includes extensive debugging output to trace the data generation process
# Supports: normal, non-normal (via gamma distribution), and uniform

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
#' @param verbose Logical. Whether to print detailed debugging information. Default is FALSE.

# Packages needed for this script:
# library(lavaan); library(covsim); library(rvinecopulib)

############################### 2. Debug Function ##############################

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
                         return.info = TRUE,
                         verbose = FALSE) {
  
  if(!is.null(seed)) set.seed(seed)
  
  # debug function for tracking issues
  print_debug <- function(...) {
    if(verbose) {
      cat(..., "\n")
    }
  }
  
  print_debug("\n", paste(rep("=", 80), collapse = ""))
  print_debug("=== Starting VITA-Based Data Generation Process ===")
  print_debug(paste(rep("=", 80), collapse = ""))
  print_debug(sprintf("Sample size: %d", N))
  print_debug(sprintf("Distribution for exogenous: %s", distr.exo))
  print_debug(sprintf("Distribution for zeta: %s", distr.zeta))
  print_debug(sprintf("Distribution for epsilon: %s", distr.epsilon))
  print_debug(sprintf("Seed: %s", ifelse(is.null(seed), "NULL", as.character(seed))))
  
  # MODEL INFORMATION
  
  print_debug("\n", paste(rep("=", 80), collapse = ""))
  print_debug("=== Step 1: Parsing Model Information ===")
  print_debug(paste(rep("=", 80), collapse = ""))
  
  # lavaan's parsed model information
  fit <- lavaan::sam(model)
  pt <- lavaan::parTable(fit)
  
  print_debug("Model fitted successfully with sam()")
  print_debug(sprintf("Number of parameters in parTable: %d", nrow(pt)))
  
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
  print_debug("\nDependent variables identified:")
  print_debug(paste("  ", model_info$structural$dependent, collapse = "\n"))
  
  # process equations, coefficients, and interactions 
  for(dv in model_info$structural$dependent) {
    print_debug(sprintf("\nProcessing equations for: %s", dv))
    eq_rows <- pt$op == "~" & pt$lhs == dv
    model_info$structural$equations[[dv]] <- pt$rhs[eq_rows]
    model_info$structural$coefficients[[dv]] <- pt$est[eq_rows]
    
    print_debug("  Predictors:")
    for(i in seq_along(model_info$structural$equations[[dv]])) {
      print_debug(sprintf("    %s: coefficient = %.4f", 
                          model_info$structural$equations[[dv]][i], 
                          model_info$structural$coefficients[[dv]][i]))
    }
    
    # interactions 
    dv_interactions <- intersect(
      pt$rhs[eq_rows],
      c(fit@pta$vnames$lv.interaction[[1]], fit@pta$vnames$ov.interaction[[1]])
    )
    if (length(dv_interactions) > 0) {
      model_info$structural$interactions[[dv]] <- dv_interactions
      print_debug("  Interactions:")
      print_debug(paste("    ", dv_interactions, collapse = "\n"))
    }
  }
  
  # all base variables
  base_vars <- c(
    fit@pta$vnames$lv.regular[[1]],
    setdiff(fit@pta$vnames$ov[[1]], fit@pta$vnames$ov.ind[[1]])
  )
  
  print_debug("\nBase variables identified:")
  print_debug(paste("  ", base_vars, collapse = "\n"))
  
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
  
  print_debug("\nDependency matrix B created")
  if(verbose) {
    print_debug("B matrix (1 = dependency, NA = no dependency):")
    print(B)
  }
  
  # ancestors 
  ancestors <- lavaan:::lav_graph_get_ancestors(B)
  names(ancestors) <- rownames(B)
  ancestor_lengths <- sapply(ancestors, length)
  
  print_debug("\nAncestor analysis:")
  for(var in names(ancestors)) {
    if(length(ancestors[[var]]) > 0) {
      print_debug(sprintf("  %s has %d ancestor(s): %s", 
                          var, length(ancestors[[var]]), 
                          paste(ancestors[[var]], collapse = ", ")))
    } else {
      print_debug(sprintf("  %s has 0 ancestors (exogenous)", var))
    }
  }
  
  # exogenous variables
  model_info$structural$exogenous <- names(ancestor_lengths)[ancestor_lengths == 0]
  # generation order
  model_info$structural$generation_order <- names(ancestor_lengths)[order(ancestor_lengths)]
  
  print_debug("\nExogenous variables identified:")
  print_debug(paste("  ", model_info$structural$exogenous, collapse = "\n"))
  
  print_debug("\nGeneration order determined:")
  for(i in seq_along(model_info$structural$generation_order)) {
    var <- model_info$structural$generation_order[i]
    print_debug(sprintf("  %d. %s (%d ancestor%s)", 
                        i, var, ancestor_lengths[var],
                        ifelse(ancestor_lengths[var] == 1, "", "s")))
  }
  
  # EXTRACT ALL RESIDUAL VARIANCES FROM MODEL
  
  print_debug("\n", paste(rep("=", 80), collapse = ""))
  print_debug("=== Step 2: Extracting Fixed Residual Variances ===")
  print_debug(paste(rep("=", 80), collapse = ""))
  
  # residual variances specified in the model
  residual_rows <- pt[pt$op == "~~" & pt$lhs == pt$rhs, ]
  all_residual_vars <- setNames(residual_rows$est, residual_rows$lhs)
  
  print_debug(sprintf("Total residual variances found: %d", length(all_residual_vars)))
  
  # distinguish structural and indicator residual variances
  indicator_vars <- fit@pta$vnames$ov.ind[[1]]
  structural_residual_vars <- all_residual_vars[!names(all_residual_vars) %in% indicator_vars]
  indicator_residual_vars <- all_residual_vars[names(all_residual_vars) %in% indicator_vars]
  
  print_debug("\nStructural residual variances:")
  if(length(structural_residual_vars) > 0) {
    for(var in names(structural_residual_vars)) {
      print_debug(sprintf("  %s: %.4f", var, structural_residual_vars[var]))
    }
  } else {
    print_debug("  None specified")
  }
  
  print_debug("\nIndicator residual variances:")
  if(length(indicator_residual_vars) > 0) {
    for(var in names(indicator_residual_vars)) {
      print_debug(sprintf("  %s: %.4f", var, indicator_residual_vars[var]))
    }
  } else {
    print_debug("  None specified")
  }
  
  # DATA GENERATION

  print_debug("\n", paste(rep("=", 80), collapse = ""))
  print_debug("=== Step 3: Initializing Data Structures ===")
  print_debug(paste(rep("=", 80), collapse = ""))
  
  # matrix for all variables
  all_vars <- c(fit@pta$vnames$lv.regular[[1]],
                fit@pta$vnames$lv.interaction[[1]],
                setdiff(fit@pta$vnames$ov[[1]], 
                        fit@pta$vnames$ov.ind[[1]]))
  
  Values <- matrix(NA, nrow = N, ncol = length(all_vars))
  colnames(Values) <- all_vars
  
  print_debug(sprintf("Values matrix initialized: %d x %d", nrow(Values), ncol(Values)))
  print_debug(sprintf("Variables: %s", paste(all_vars, collapse = ", ")))
  
  # GENERATE EXOGENOUS VARIABLES USING VITA
 
  print_debug("\n", paste(rep("=", 80), collapse = ""))
  print_debug("=== Step 4: Generating Exogenous Variables (VITA) ===")
  print_debug(paste(rep("=", 80), collapse = ""))
  
  exo_vars <- model_info$structural$exogenous
  n_exo <- length(exo_vars)
  
  print_debug(sprintf("Number of exogenous variables: %d", n_exo))
  print_debug(sprintf("Exogenous variables: %s", paste(exo_vars, collapse = ", ")))
  
  # get covariance matrix from model
  psi_matrix <- lavaan::lavInspect(fit, "est")$psi
  exo.vcov <- psi_matrix[exo_vars, exo_vars, drop = FALSE]
  
  print_debug("\nModel-specified covariance matrix for exogenous variables:")
  if(verbose) {
    print(round(exo.vcov, 4))
  }
  
  # fill in missing diagonal elements
  for(i in 1:nrow(exo.vcov)) { 
    if(is.na(exo.vcov[i,i])) {
      exo.vcov[i,i] <- 1
      print_debug(sprintf("  WARNING: Missing variance for %s - set to 1", exo_vars[i]))
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
    
    print_debug(sprintf("\nGenerating exogenous variables with distribution: %s", distr.exo))
    
    if(distr.exo == "normal") {
      # normal distributions with variance matching model
      print_debug("  Setting up normal marginals...")
      margins_list <- lapply(1:n_exo, function(i) {
        margin <- list(
          distr = "norm",
          mean = 0,  # adjusted later if exo.mean provided
          sd = sqrt(exo.vcov[i, i])
        )
        print_debug(sprintf("    %s: mean=%.4f, sd=%.4f, var=%.4f", 
                            exo_vars[i], margin$mean, margin$sd, exo.vcov[i, i]))
        margin
      })
      
    } else if(distr.exo == "nonnormal") {
      # gamma distributions (skewed, heavy-tailed)
      print_debug("  Setting up gamma (non-normal) marginals...")
      
      if(is.null(nonnormal.shape)) {
        nonnormal.shape <- rep(1, n_exo)
        print_debug("    Using default shape parameters: all set to 1")
      }
      if(is.null(nonnormal.rate)) {
        nonnormal.rate <- rep(1, n_exo)
        print_debug("    Using default rate parameters: all set to 1")
      }
      
      if(length(nonnormal.shape) != n_exo) {
        stop(paste("nonnormal.shape must have length", n_exo))
      }
      if(length(nonnormal.rate) != n_exo) {
        stop(paste("nonnormal.rate must have length", n_exo))
      }
      
      print_debug("  Gamma distribution parameters:")
      margins_list <- lapply(1:n_exo, function(i) {
        margin <- list(
          distr = "gamma",
          shape = nonnormal.shape[i],
          rate = nonnormal.rate[i]
        )
        
        # theoretical moments
        gamma_mean <- margin$shape / margin$rate
        gamma_var <- margin$shape / margin$rate^2
        gamma_skew <- 2 / sqrt(margin$shape)
        gamma_kurt <- 6 / margin$shape
        
        print_debug(sprintf("    %s: shape=%.2f, rate=%.2f", exo_vars[i], margin$shape, margin$rate))
        print_debug(sprintf("      Natural mean=%.4f, var=%.4f", gamma_mean, gamma_var))
        print_debug(sprintf("      Skewness=%.4f, Excess kurtosis=%.4f", gamma_skew, gamma_kurt))
        
        margin
      })
      
    } else if(distr.exo == "uniform") {
      # uniform distributions with variance matching model
      print_debug("  Setting up uniform marginals...")
      print_debug("  Using symmetric bounds: [-sqrt(3*v), sqrt(3*v)] for variance v")
      
      margins_list <- lapply(1:n_exo, function(i) {
        var_target <- exo.vcov[i, i]
        half_range <- sqrt(3 * var_target)
        margin <- list(
          distr = "unif",
          min = -half_range,
          max = half_range
        )
        
        unif_var <- (margin$max - margin$min)^2 / 12
        
        print_debug(sprintf("    %s: min=%.4f, max=%.4f, range=%.4f", 
                            exo_vars[i], margin$min, margin$max, margin$max - margin$min))
        print_debug(sprintf("      Target var=%.4f, Actual var=%.4f", var_target, unif_var))
        
        margin
      })
      
    } else {
      stop(paste("Unknown distr.exo option:", distr.exo))
    }
    
    # calculate variances from marginals
    marginal_vars <- sapply(margins_list, calc_marginal_var)
    
    print_debug("\nMarginal variances from distribution specifications:")
    for(i in 1:n_exo) {
      print_debug(sprintf("  %s: %.4f", exo_vars[i], marginal_vars[i]))
    }
    
    # get correlation structure from model covariance
    D_diag <- diag(exo.vcov)
    cor_matrix <- exo.vcov / sqrt(outer(D_diag, D_diag))
    
    print_debug("\nCorrelation structure from model:")
    if(verbose) {
      print(round(cor_matrix, 4))
    }
    
    # build covariance matrix with marginal variances
    D <- diag(sqrt(marginal_vars))
    adjusted_cov <- D %*% cor_matrix %*% D
    
    print_debug("\nAdjusted covariance matrix for VITA:")
    if(verbose) {
      print(round(adjusted_cov, 4))
    }
    
    # create VITA distribution
    print_debug("\nCreating VITA distribution...")
    vitadist <- covsim::vita(margins_list, adjusted_cov, verbose = FALSE, Nmax = 10^6)
    print_debug("  VITA distribution created successfully")
    
    # generate data
    print_debug(sprintf("  Generating %d observations...", N))
    EXO <- rvinecopulib::rvine(n = N, vine = vitadist)
    print_debug("  Data generation complete")
    
    # check generated data before shifting
    print_debug("\nGenerated data (before shifting):")
    for(i in 1:ncol(EXO)) {
      print_debug(sprintf("  Variable %d: mean=%.4f, var=%.4f, min=%.4f, max=%.4f", 
                          i, mean(EXO[,i]), var(EXO[,i]), min(EXO[,i]), max(EXO[,i])))
    }
    
    # SHIFTING FOR NON-NORMAL DISTRIBUTIONS
    if(distr.exo == "nonnormal") {
      print_debug("\nShifting non-normal distributions to center at 0...")
      for(i in 1:ncol(EXO)) {
        natural_mean <- nonnormal.shape[i] / nonnormal.rate[i]
        mean_before <- mean(EXO[, i])
        EXO[, i] <- EXO[, i] - natural_mean
        mean_after <- mean(EXO[, i])
        
        print_debug(sprintf("  %s: natural mean=%.4f, shifted %.4f -> %.6f", 
                            exo_vars[i], natural_mean, mean_before, mean_after))
      }
    } else {
      print_debug("\nNo shifting needed (normal/uniform already centered at 0)")
    }
    
    # check generated data after shifting
    print_debug("\nGenerated data (after shifting):")
    for(i in 1:ncol(EXO)) {
      gen_skew <- moments::skewness(EXO[,i])
      gen_kurt <- moments::kurtosis(EXO[,i]) - 3  # excess kurtosis
      
      print_debug(sprintf("  %s: mean=%.6f, var=%.4f, sd=%.4f", 
                          exo_vars[i], mean(EXO[,i]), var(EXO[,i]), sd(EXO[,i])))
      print_debug(sprintf("         skewness=%.4f, excess kurtosis=%.4f", 
                          gen_skew, gen_kurt))
    }
    
    EXO
  }
  
  # generate exogenous variables
  EXO <- generate_exo()
  
  # ensure it's a matrix
  if(!is.matrix(EXO)) {
    EXO <- as.matrix(EXO)
  }
  
  colnames(EXO) <- exo_vars
  
  print_debug("\nExogenous variable matrix created:")
  print_debug(sprintf("  Dimensions: %d x %d", nrow(EXO), ncol(EXO)))
  print_debug(sprintf("  Variables: %s", paste(colnames(EXO), collapse = ", ")))
  
  # CENTER EXOGENOUS VARIABLES (IF REQUESTED)

  print_debug("\n", paste(rep("=", 80), collapse = ""))
  print_debug("=== Step 5: Centering/Adjusting Exogenous Variables ===")
  print_debug(paste(rep("=", 80), collapse = ""))
  
  # center exogenous variables
  eta_cols <- grep("eta", exo_vars, value = TRUE)
  manifest_cols <- setdiff(exo_vars, eta_cols)
  
  print_debug(sprintf("Latent exogenous variables (eta): %s", 
                      ifelse(length(eta_cols) > 0, paste(eta_cols, collapse = ", "), "none")))
  print_debug(sprintf("Manifest exogenous variables: %s", 
                      ifelse(length(manifest_cols) > 0, paste(manifest_cols, collapse = ", "), "none")))
  
  if(center.exogenous.latent && length(eta_cols) > 0) {
    print_debug("\nCentering latent exogenous variables:")
    for(col in eta_cols) {
      mean_before <- mean(EXO[, col])
      first_val_before <- EXO[1, col]
      
      EXO[, col] <- EXO[, col] - mean_before
      
      mean_after <- mean(EXO[, col])
      first_val_after <- EXO[1, col]
      
      print_debug(sprintf("  %s:", col))
      print_debug(sprintf("    Mean before: %.6f", mean_before))
      print_debug(sprintf("    Mean after: %.10f", mean_after))
      print_debug(sprintf("    First value: %.6f -> %.6f (diff: %.6f)", 
                          first_val_before, first_val_after, first_val_after - first_val_before))
    }
  } else if(!is.null(exo.mean) && length(eta_cols) > 0) {
    print_debug("\nAdjusting latent exogenous means:")
    for(i in seq_along(eta_cols)) {
      col <- eta_cols[i]
      mean_before <- mean(EXO[, col])
      first_val_before <- EXO[1, col]
      
      EXO[, col] <- EXO[, col] + exo.mean[i]
      
      mean_after <- mean(EXO[, col])
      first_val_after <- EXO[1, col]
      
      print_debug(sprintf("  %s:", col))
      print_debug(sprintf("    Added mean: %.6f", exo.mean[i]))
      print_debug(sprintf("    Mean before: %.6f", mean_before))
      print_debug(sprintf("    Mean after: %.6f", mean_after))
      print_debug(sprintf("    First value: %.6f -> %.6f", first_val_before, first_val_after))
    }
  } else {
    print_debug("\nNo centering/adjustment for latent exogenous variables")
  }
  
  if(center.exogenous.manifest && length(manifest_cols) > 0) {
    print_debug("\nCentering manifest exogenous variables:")
    for(col in manifest_cols) {
      mean_before <- mean(EXO[, col])
      first_val_before <- EXO[1, col]
      
      EXO[, col] <- EXO[, col] - mean_before
      
      mean_after <- mean(EXO[, col])
      first_val_after <- EXO[1, col]
      
      print_debug(sprintf("  %s:", col))
      print_debug(sprintf("    Mean before: %.6f", mean_before))
      print_debug(sprintf("    Mean after: %.10f", mean_after))
      print_debug(sprintf("    First value: %.6f -> %.6f", first_val_before, first_val_after))
    }
  } else if(!is.null(exo.mean) && length(manifest_cols) > 0) {
    print_debug("\nAdjusting manifest exogenous means:")
    offset <- length(eta_cols)
    for(i in seq_along(manifest_cols)) {
      col <- manifest_cols[i]
      mean_before <- mean(EXO[, col])
      first_val_before <- EXO[1, col]
      
      EXO[, col] <- EXO[, col] + exo.mean[offset + i]
      
      mean_after <- mean(EXO[, col])
      first_val_after <- EXO[1, col]
      
      print_debug(sprintf("  %s:", col))
      print_debug(sprintf("    Added mean: %.6f", exo.mean[offset + i]))
      print_debug(sprintf("    Mean before: %.6f", mean_before))
      print_debug(sprintf("    Mean after: %.6f", mean_after))
      print_debug(sprintf("    First value: %.6f -> %.6f", first_val_before, first_val_after))
    }
  } else {
    print_debug("\nNo centering/adjustment for manifest exogenous variables")
  }
  
  # store in Values matrix
  Values[, exo_vars] <- EXO
  
  print_debug("\nFinal exogenous variable statistics:")
  for(var in exo_vars) {
    print_debug(sprintf("  %s: mean=%.6f, var=%.4f, sd=%.4f", 
                        var, mean(Values[, var]), var(Values[, var]), sd(Values[, var])))
  }
  
  # GENERATE DEPENDENT VARIABLES WITH FIXED RESIDUAL VARIANCES
  
  print_debug("\n", paste(rep("=", 80), collapse = ""))
  print_debug("=== Step 6: Generating Dependent Variables ===")
  print_debug(paste(rep("=", 80), collapse = ""))
  
  observed_R2 <- list()
  deterministic_vars <- list()
  stored_zetas <- list()  # checking zeta independence
  
  intercepts <- lavaan::lavInspect(fit, "est")$alpha
  
  if(!is.null(intercepts)) {
    print_debug("Intercepts from model:")
    for(var in rownames(intercepts)) {
      if(!is.na(intercepts[var, 1])) {
        print_debug(sprintf("  %s: %.4f", var, intercepts[var, 1]))
      }
    }
  } else {
    print_debug("No intercepts specified in model")
  }
  
  for(var in model_info$structural$generation_order) { 
    if(var %in% model_info$structural$dependent) {
      
      print_debug(sprintf("\n--- Generating: %s ---", var))
      
      terms <- model_info$structural$equations[[var]]
      
      # handle interactions first
      if(!is.null(model_info$structural$interactions[[var]])) {
        print_debug("  Creating interaction terms:")
        for(inter in model_info$structural$interactions[[var]]) {
          if(all(is.na(Values[, inter]))) {
            components <- unlist(strsplit(inter, ":"))
            
            print_debug(sprintf("\n    --- Creating %s = %s * %s ---", inter, components[1], components[2]))
            
            # components exist and have values
            for(comp in components) {
              if(all(is.na(Values[, comp]))) {
                stop(sprintf("ERROR: Component %s has no values for interaction %s", comp, inter))
              }
              print_debug(sprintf("      Component %s:", comp))
              print_debug(sprintf("        Mean = %.6f", mean(Values[, comp])))
              print_debug(sprintf("        Var  = %.6f", var(Values[, comp])))
              print_debug(sprintf("        First 3 values: %.6f, %.6f, %.6f", 
                                  Values[1, comp], Values[2, comp], Values[3, comp]))
            }
            
            # the actual multiplication for first 3 observations
            print_debug(sprintf("\n      Multiplication verification (first 3 rows):"))
            for(i in 1:3) {
              val1 <- Values[i, components[1]]
              val2 <- Values[i, components[2]]
              product <- val1 * val2
              print_debug(sprintf("        Row %d: %.6f * %.6f = %.6f", i, val1, val2, product))
            }
            
            Values[, inter] <- Values[, components[1]] * Values[, components[2]]
            
            print_debug(sprintf("\n      Interaction %s created:", inter))
            print_debug(sprintf("        Mean = %.6f", mean(Values[, inter])))
            print_debug(sprintf("        Var  = %.6f", var(Values[, inter])))
            print_debug(sprintf("        First 3 values: %.6f, %.6f, %.6f", 
                                Values[1, inter], Values[2, inter], Values[3, inter]))
            
            if(center.lv.prod) {
              mean_before <- mean(Values[, inter])
              Values[, inter] <- Values[, inter] - mean_before
              mean_after <- mean(Values[, inter])
              print_debug(sprintf("\n      Centered interaction term:"))
              print_debug(sprintf("        Mean before centering: %.6f", mean_before))
              print_debug(sprintf("        Mean after centering:  %.10f", mean_after))
              print_debug(sprintf("        First 3 values after:  %.6f, %.6f, %.6f", 
                                  Values[1, inter], Values[2, inter], Values[3, inter]))
            }
          }
        }
      }
      
      # coefficients
      equation_coefs <- model_info$structural$coefficients[[var]]
      
      print_debug("\n  === Structural Equation ===")
      print_debug(sprintf("  %s = ", var))
      for(i in seq_along(terms)) {
        print_debug(sprintf("    + %.6f * %s", equation_coefs[i], terms[i]))
      }
      
      # all predictors exist?
      missing_terms <- terms[!terms %in% colnames(Values)]
      if(length(missing_terms) > 0) {
        stop(sprintf("ERROR: Terms not found for %s: %s", 
                     var, paste(missing_terms, collapse = ", ")))
      }
      
      # predictor values for first 3 observations
      print_debug("\n  Predictor values (first 3 rows):")
      for(term in terms) {
        print_debug(sprintf("    %s: %.6f, %.6f, %.6f", 
                            term, Values[1, term], Values[2, term], Values[3, term]))
      }
      
      # deterministic part
      deterministic_part <- Values[, terms, drop = FALSE] %*% equation_coefs
      
      print_debug("\n  === Deterministic Part Calculation (BEFORE Intercept) ===")
      
      # manual calculation verification for first 3 rows
      print_debug("  Manual calculation verification (first 3 rows):")
      for(i in 1:3) {
        calc_parts <- sprintf("%.6f*%.6f", equation_coefs, Values[i, terms])
        calc_string <- paste(calc_parts, collapse = " + ")
        calc_result <- sum(equation_coefs * Values[i, terms])
        print_debug(sprintf("    Row %d: %s = %.6f", i, calc_string, calc_result))
        print_debug(sprintf("           Matrix result: %.6f (diff: %.10f)", 
                            deterministic_part[i], deterministic_part[i] - calc_result))
      }
      
      print_debug("\n  Statistics (before intercept):")
      print_debug(sprintf("    Mean:     %.6f", mean(deterministic_part)))
      print_debug(sprintf("    Variance: %.6f", var(deterministic_part)))
      print_debug(sprintf("    SD:       %.6f", sd(deterministic_part)))
      print_debug(sprintf("    Min:      %.6f", min(deterministic_part)))
      print_debug(sprintf("    Max:      %.6f", max(deterministic_part)))
      
      # add intercept if present
      if (!is.null(intercepts) && !is.na(intercepts[var,1])) {
        intercept_val <- intercepts[var,1]
        print_debug(sprintf("\n  === Adding Intercept: %.6f ===", intercept_val))
        
        print_debug("  First 3 values before intercept:")
        for(i in 1:3) {
          print_debug(sprintf("    Row %d: %.6f", i, deterministic_part[i]))
        }
        
        deterministic_part <- intercept_val + deterministic_part
        
        print_debug("\n  First 3 values AFTER adding intercept:")
        for(i in 1:3) {
          print_debug(sprintf("    Row %d: %.6f + %.6f = %.6f", 
                              i, intercept_val, deterministic_part[i] - intercept_val, deterministic_part[i]))
        }
        
        print_debug("\n  Statistics (after intercept):")
        print_debug(sprintf("    Mean:     %.6f (should be ≈ %.6f)", 
                            mean(deterministic_part), mean(deterministic_part - intercept_val) + intercept_val))
        print_debug(sprintf("    Variance: %.6f (should be unchanged)", var(deterministic_part)))
        print_debug(sprintf("    SD:       %.6f", sd(deterministic_part)))
        print_debug(sprintf("    Min:      %.6f", min(deterministic_part)))
        print_debug(sprintf("    Max:      %.6f", max(deterministic_part)))
      } else {
        print_debug("\n  No intercept to add")
      }
      
      # variance of deterministic part
      var_det <- var(as.vector(deterministic_part)) * (N-1)/N
      deterministic_vars[[var]] <- var_det
      
      print_debug(sprintf("\n  Deterministic variance (population): %.6f", var_det))
      
      # fixed residual variance from model
      resid_var <- structural_residual_vars[var]
      if(is.na(resid_var)) {
        resid_var <- 1
        warning(paste("No residual variance found for", var, "- using 1"))
        print_debug(sprintf("  WARNING: No residual variance found for %s - using 1", var))
      } else {
        print_debug(sprintf("  Fixed residual variance from model: %.6f", resid_var))
      }
      
      # expected r-squared with fixed variance
      expected_r2 <- var_det / (var_det + resid_var)
      print_debug(sprintf("  Expected R² = %.6f / (%.6f + %.6f) = %.6f", 
                          var_det, var_det, resid_var, expected_r2))
      
      # generate residual with fixed variance
      print_debug(sprintf("\n  === Generating %s Residual (Zeta) ===", distr.zeta))
      print_debug(sprintf("  Target variance: %.6f", resid_var))
      print_debug(sprintf("  Target SD:       %.6f", sqrt(resid_var)))
      
      zeta <- switch(distr.zeta,
                     "normal" = rnorm(N, 0, sqrt(resid_var)),
                     "exp.rate1" = rexp(N, rate = 1/sqrt(resid_var)) - sqrt(resid_var),
                     stop(paste("wrong option for distr.zeta:", distr.zeta)))
      
      # zeta properties using moments package
      if(requireNamespace("moments", quietly = TRUE)) {
        zeta_skew <- moments::skewness(zeta)
        zeta_kurt <- moments::kurtosis(zeta) - 3  # excess kurtosis
      } else {
        zeta_skew <- NA
        zeta_kurt <- NA
      }
      
      print_debug("  Residual (zeta) diagnostics:")
      print_debug(sprintf("    Mean:     %.6f (target: 0.000000)", mean(zeta)))
      print_debug(sprintf("    Variance: %.6f (target: %.6f, diff: %.6f)", 
                          var(zeta), resid_var, var(zeta) - resid_var))
      print_debug(sprintf("    SD:       %.6f (target: %.6f)", sd(zeta), sqrt(resid_var)))
      if(!is.na(zeta_skew)) {
        print_debug(sprintf("    Skewness:       %.6f", zeta_skew))
        print_debug(sprintf("    Excess kurtosis: %.6f", zeta_kurt))
      }
      print_debug(sprintf("    First 3 values: %.6f, %.6f, %.6f", zeta[1], zeta[2], zeta[3]))
      
      # store for independence check
      stored_zetas[[var]] <- zeta
      
      # combine deterministic and stochastic parts
      print_debug("\n  === Combining Deterministic and Stochastic Parts ===")
      print_debug(sprintf("  Formula: %s = deterministic_part + zeta", var))
      
      print_debug("\n  Verification (first 5 rows):")
      print_debug("  Row | Deterministic |    Zeta    |   Combined  | Check")
      print_debug("  ----|---------------|------------|-------------|-------")
      for(i in 1:5) {
        det_val <- deterministic_part[i]
        zeta_val <- zeta[i]
        combined <- det_val + zeta_val
        print_debug(sprintf("  %3d | %13.6f | %10.6f | %11.6f | %.10f", 
                            i, det_val, zeta_val, combined, combined - (det_val + zeta_val)))
      }
      
      Values[, var] <- deterministic_part + zeta
      
      print_debug("\n  Combined variable statistics (before any centering):")
      print_debug(sprintf("    Mean:     %.6f", mean(Values[, var])))
      print_debug(sprintf("    Variance: %.6f", var(Values[, var])))
      print_debug(sprintf("    SD:       %.6f", sd(Values[, var])))
      
      # verify the combination
      manual_check <- deterministic_part + zeta
      max_diff <- max(abs(Values[, var] - manual_check))
      print_debug(sprintf("    Max difference from manual calculation: %.15f", max_diff))
      if(max_diff > 1e-10) {
        print_debug("    WARNING: Difference exceeds numerical precision threshold!")
      }
      
      if(center.lv.dependent) {
        mean_before <- mean(Values[, var])
        first_val_before <- Values[1, var]
        
        Values[, var] <- Values[, var] - mean_before
        
        mean_after <- mean(Values[, var])
        first_val_after <- Values[1, var]
        
        print_debug(sprintf("\n  Centered dependent variable %s:", var))
        print_debug(sprintf("    Mean before: %.6f", mean_before))
        print_debug(sprintf("    Mean after: %.10f", mean_after))
        print_debug(sprintf("    First value: %.6f -> %.6f", first_val_before, first_val_after))
      }
      
      # observed r-squared
      total_var <- var(Values[, var]) * (N-1)/N
      observed_R2[[var]] <- var_det / total_var
      
      print_debug(sprintf("\n  Final diagnostics for %s:", var))
      print_debug(sprintf("    Mean: %.6f", mean(Values[, var])))
      print_debug(sprintf("    Total variance: %.6f", total_var))
      print_debug(sprintf("    Deterministic variance: %.6f", var_det))
      print_debug(sprintf("    Residual variance: %.6f", total_var - var_det))
      print_debug(sprintf("    Observed R²: %.6f (expected: %.6f, diff: %.6f)", 
                          observed_R2[[var]], expected_r2, observed_R2[[var]] - expected_r2))
    }
  }
  
  # check zeta independence
  if(verbose && length(stored_zetas) > 1) {
    print_debug("\n", paste(rep("=", 80), collapse = ""))
    print_debug("=== Checking Residual (Zeta) Independence ===")
    print_debug(paste(rep("=", 80), collapse = ""))
    
    zeta_matrix <- do.call(cbind, stored_zetas)
    colnames(zeta_matrix) <- names(stored_zetas)
    
    zeta_cor <- cor(zeta_matrix)
    print_debug("Zeta correlations (should be near 0):")
    print(round(zeta_cor, 4))
    
    # flag large correlations
    for(i in 1:(ncol(zeta_cor)-1)) {
      for(j in (i+1):ncol(zeta_cor)) {
        if(abs(zeta_cor[i,j]) > 0.1) {
          print_debug(sprintf("  WARNING: Large correlation between residuals of %s and %s: %.4f",
                              colnames(zeta_cor)[i], colnames(zeta_cor)[j], zeta_cor[i,j]))
        }
      }
    }
    
    if(all(abs(zeta_cor[lower.tri(zeta_cor)]) < 0.1)) {
      print_debug("  ✓ All residual correlations are acceptably small")
    }
  }
  
  # MEASUREMENT PART - WITH FIXED VARIANCES FROM MODEL

  print_debug("\n", paste(rep("=", 80), collapse = ""))
  print_debug("=== Step 7: Generating Indicators (Measurement Model) ===")
  print_debug(paste(rep("=", 80), collapse = ""))
  
  lambda <- lavaan::lavInspect(fit, "est")$lambda
  LAMBDA <- lambda[startsWith(rownames(lambda), "x"), fit@pta$vnames$lv.regular[[1]], 
                   drop = FALSE]
  
  print_debug("Lambda matrix (factor loadings):")
  if(verbose) {
    print(round(LAMBDA, 4))
  }
  
  # analyze loading structure
  print_debug("\nLoading structure:")
  for(i in seq_along(fit@pta$vnames$lv.regular[[1]])) {
    eta <- fit@pta$vnames$lv.regular[[1]][i]
    indicators <- rownames(LAMBDA)[LAMBDA[, i] != 0]
    loadings <- LAMBDA[LAMBDA[, i] != 0, i]
    
    print_debug(sprintf("  %s has %d indicator(s):", eta, length(indicators)))
    for(j in seq_along(indicators)) {
      print_debug(sprintf("    %s: loading = %.4f", indicators[j], loadings[j]))
    }
  }
  
  # measurement errors with fixed variances from model
  n_indicators <- nrow(LAMBDA)
  THETA <- matrix(NA, nrow = N, ncol = n_indicators)
  colnames(THETA) <- rownames(LAMBDA)
  
  print_debug("\nGenerating measurement errors:")
  print_debug(sprintf("Distribution: %s", distr.epsilon))
  
  for(i in seq_len(n_indicators)) {
    ind_name <- rownames(LAMBDA)[i]
    # variance specified in the model
    ind_var <- indicator_residual_vars[ind_name]
    if(is.na(ind_var)) {
      ind_var <- 1
      warning(paste("No residual variance found for", ind_name, "- using 1"))
      print_debug(sprintf("  WARNING: No residual variance for %s - using 1", ind_name))
    } else {
      print_debug(sprintf("  %s: model-specified error variance = %.4f", ind_name, ind_var))
    }
    
    if(distr.epsilon == "normal") {
      THETA[, i] <- rnorm(N, 0, sqrt(ind_var))
    } else {
      THETA[, i] <- rexp(N, rate = 1/sqrt(ind_var)) - sqrt(ind_var)
    }
    
    # error properties
    print_debug(sprintf("    Generated error: mean=%.6f, var=%.6f", 
                        mean(THETA[, i]), var(THETA[, i])))
  }
  
  # indicators: Y = Lambda * Eta + Theta
  print_debug("\nCalculating indicator values: Y = Lambda * Eta + Theta")
  
  Y <- as.matrix(Values[, fit@pta$vnames$lv.regular[[1]]]) %*% t(LAMBDA) + THETA 
  colnames(Y) <- rownames(LAMBDA)
  
  # diagnostics before centering
  print_debug("\nIndicator diagnostics (before centering):")
  for(ind in colnames(Y)) {
    print_debug(sprintf("  %s: mean=%.6f, var=%.4f, sd=%.4f", 
                        ind, mean(Y[, ind]), var(Y[, ind]), sd(Y[, ind])))
  }
  
  # center indicators 
  if(center.indicators) {
    print_debug("\nCentering indicators:")
    for(ind in colnames(Y)) {
      mean_before <- mean(Y[, ind])
      Y[, ind] <- Y[, ind] - mean_before
      mean_after <- mean(Y[, ind])
      
      print_debug(sprintf("  %s: %.6f -> %.10f", ind, mean_before, mean_after))
    }
  } else {
    print_debug("\nNo centering applied to indicators")
  }
  
  # VERIFY RELIABILITIES

  print_debug("\n", paste(rep("=", 80), collapse = ""))
  print_debug("=== Step 8: Verifying Indicator Reliabilities ===")
  print_debug(paste(rep("=", 80), collapse = ""))
  
  # calculate observed reliabilities
  observed_reliabilities <- list()
  
  for(i in seq_along(fit@pta$vnames$lv.regular[[1]])) {
    eta <- fit@pta$vnames$lv.regular[[1]][i]
    eta_var <- var(Values[, eta]) * (N-1)/N
    
    print_debug(sprintf("\nReliabilities for indicators of %s:", eta))
    print_debug(sprintf("  Latent variable variance: %.4f", eta_var))
    
    ind_idx <- which(LAMBDA[, i] != 0)
    if(length(ind_idx) > 0) {
      reliabilities <- numeric(length(ind_idx))
      for(j in seq_along(ind_idx)) {
        ind_name <- rownames(LAMBDA)[ind_idx[j]]
        loading <- LAMBDA[ind_idx[j], i]
        error_var <- indicator_residual_vars[ind_name]
        
        # theoretical reliability from model specification
        theoretical_rel <- (loading^2 * eta_var) / (loading^2 * eta_var + error_var)
        
        # observed reliability
        cor_with_lv <- cor(Y[, ind_name], Values[, eta])
        observed_rel <- cor_with_lv^2
        
        reliabilities[j] <- theoretical_rel
        
        print_debug(sprintf("  %s:", ind_name))
        print_debug(sprintf("    Loading: %.4f", loading))
        print_debug(sprintf("    Error variance: %.4f", error_var))
        print_debug(sprintf("    True score variance: %.4f", loading^2 * eta_var))
        print_debug(sprintf("    Observed variance: %.4f", var(Y[, ind_name])))
        print_debug(sprintf("    Theoretical reliability: %.4f", theoretical_rel))
        print_debug(sprintf("    Observed reliability (r²): %.4f", observed_rel))
        print_debug(sprintf("    Difference: %.4f", observed_rel - theoretical_rel))
        
        if(abs(observed_rel - theoretical_rel) > 0.05) {
          print_debug("    WARNING: Large discrepancy between theoretical and observed reliability!")
        }
      }
      names(reliabilities) <- rownames(LAMBDA)[ind_idx]
      observed_reliabilities[[eta]] <- reliabilities
    }
  }
  
  # FINALIZE DATASET
 
  print_debug("\n", paste(rep("=", 80), collapse = ""))
  print_debug("=== Step 9: Finalizing Dataset ===")
  print_debug(paste(rep("=", 80), collapse = ""))
  
  if(add.eta) {
    print_debug("Including latent variables in final dataset")
    Results <- cbind(Values, Y)
  } else {
    print_debug("Excluding latent variables from final dataset")
    manifest_vars <- setdiff(fit@pta$vnames$ov[[1]], fit@pta$vnames$ov.ind[[1]])
    manifest_vars <- manifest_vars[!grepl("eta", manifest_vars)]
    
    print_debug(sprintf("Manifest variables to include: %s", 
                        ifelse(length(manifest_vars) > 0, 
                               paste(manifest_vars, collapse = ", "), 
                               "none")))
    
    if(length(manifest_vars) > 0) {
      Values_indicators <- Values[, manifest_vars, drop = FALSE]
      Results <- cbind(Values_indicators, Y)
    } else {
      Results <- Y
    }
  }
  
  # final validation
  print_debug("\nFinal dataset validation:")
  print_debug(sprintf("  Dimensions: %d x %d", nrow(Results), ncol(Results)))
  print_debug(sprintf("  Variables: %s", paste(colnames(Results), collapse = ", ")))
  
  # check for NAs
  na_count <- colSums(is.na(Results))
  if(any(na_count > 0)) {
    print_debug("  WARNING: NA values detected in final dataset:")
    for(var in names(na_count)[na_count > 0]) {
      print_debug(sprintf("    %s: %d NAs (%.2f%%)", var, na_count[var], 100*na_count[var]/nrow(Results)))
    }
  } else {
    print_debug("  ✓ No NA values detected")
  }
  
  # check for infinite values
  inf_count <- colSums(is.infinite(as.matrix(Results)))
  if(any(inf_count > 0)) {
    print_debug("  WARNING: Infinite values detected in final dataset:")
    for(var in names(inf_count)[inf_count > 0]) {
      print_debug(sprintf("    %s: %d Inf values", var, inf_count[var]))
    }
  } else {
    print_debug("  ✓ No infinite values detected")
  }
  
  # summary stats
  print_debug("\nFinal variable summary statistics:")
  for(var in colnames(Results)) {
    var_mean <- mean(Results[, var])
    var_sd <- sd(Results[, var])
    var_min <- min(Results[, var])
    var_max <- max(Results[, var])
    
    print_debug(sprintf("  %s: mean=%.6f, sd=%.4f, range=[%.4f, %.4f]", 
                        var, var_mean, var_sd, var_min, var_max))
  }
  
  result <- as.data.frame(Results)
  
  if(return.info) {
    print_debug("\nAdding attributes to result:")
    attr(result, "observed_R2") <- observed_R2
    attr(result, "observed_reliabilities") <- observed_reliabilities
    attr(result, "fixed_residual_variances") <- all_residual_vars
    attr(result, "model_info") <- model_info
    
    print_debug("  ✓ observed_R2")
    print_debug("  ✓ observed_reliabilities")
    print_debug("  ✓ fixed_residual_variances")
    print_debug("  ✓ model_info")
  }
  
  # SUMMARY
  
  print_debug("\n", paste(rep("=", 80), collapse = ""))
  print_debug("=== GENERATION SUMMARY ===")
  print_debug(paste(rep("=", 80), collapse = ""))
  print_debug(sprintf("Dataset dimensions: %d x %d", nrow(result), ncol(result)))
  print_debug("\nVariable types generated:")
  print_debug(sprintf("  - Exogenous: %d", length(model_info$structural$exogenous)))
  print_debug(sprintf("  - Dependent: %d", length(model_info$structural$dependent)))
  print_debug(sprintf("  - Indicators: %d", ncol(Y)))
  
  n_interactions <- sum(sapply(model_info$structural$interactions, length))
  if(n_interactions > 0) {
    print_debug(sprintf("  - Interactions: %d", n_interactions))
  }
  
  if(length(observed_R2) > 0) {
    print_debug("\nR² values achieved:")
    for(var in names(observed_R2)) {
      print_debug(sprintf("  %s: %.4f", var, observed_R2[[var]]))
    }
  }
  
  if(length(observed_reliabilities) > 0) {
    print_debug("\nAverage reliabilities by factor:")
    for(eta in names(observed_reliabilities)) {
      avg_rel <- mean(observed_reliabilities[[eta]])
      print_debug(sprintf("  %s: %.4f", eta, avg_rel))
    }
  }
  
  print_debug("\n", paste(rep("=", 80), collapse = ""))
  print_debug("=== Data Generation Complete ===")
  print_debug(paste(rep("=", 80), collapse = ""))
  
  result
}


# example:
# 
# 
#data <- GenerateData(
#   model = all_models[["normal_rel08"]],
#   N = 10000,
#   distr.exo = "nonnormal",
#   nonnormal.shape = c(1, 1),
#   nonnormal.rate = c(1, 1),
#   verbose = TRUE  
#)
# 
# attr(data, "observed_R2")
# attr(data, "observed_reliabilities")
