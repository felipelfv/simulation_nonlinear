############################ 1. General Information ############################

# See README file for more information concerning this file. 

# This file contains the code used to estimate the different approaches. We use
# modsem for all except the SAM approach. Thus, we have the DBLCENT, LMS, and 
# QML estimated through modsem and the SAM approach with lavaan. 
# All functions return the full parameter table.

# Importantly, the helper function is used to extract only the relevant 
# estimates. In our case, the structural model. 

############################### 2. Functions ###################################

#### 2.1. Product Indicator Approach with double mean centering (DBLCENT) ####
method_dblcent <- function(
    Data = NULL,
    model.fit = NULL,
    robust.se = FALSE,
    match = NULL,
    bounds = FALSE
) {
  # for matching, robust, and bounds
  args <- list(
    model.syntax = model.fit,
    data = Data,
    method = "dblcent",
    bounds = bounds
  )
  
  # robust SE if requested
  if (robust.se) {
    args$se <- "robust.huber.white"
  }
  
  if (!is.null(match)) {
    args$match <- match
    out <- do.call(modsem::modsem_pi, args)
  } else {
    out <- do.call(modsem::modsem, args)
  }
  
  # coefficient parameter table
  out$coefParTable
}

#### 2.2 Structural-After-Measurement (SAM) Approach ####
method_sam <- function(Data = NULL, estimator = "ML",
                       joint = FALSE, add.attr = FALSE, 
                       model.fit = NULL,
                       mm.list = NULL, 
                       robust.se = FALSE) {
  
  # Adjust estimator if robust SE requested
  if (robust.se) {
    estimator <- "MLR"  # or "MLM" depending on your needs
  }
  
  out <- lavaan::sam(model.fit, data = Data, se = "local",
                     mm.args = list(estimator = estimator),
                     mm.list = mm.list)
  
  # parameter estimates table
  parameterEstimates(out, remove.step1 = FALSE)
}

#### 2.3 Distribution Analytic Approaches (LMS and QML) ####
method_analytic <- function(Data = NULL, model.fit = NULL, 
                            standardized = FALSE, method = "lms", robust.se = FALSE) {
  
  # Build argument list
  args <- list(model.syntax = model.fit, data = Data, method = method)
  
  # Add robust SE if requested
  if (robust.se) {
    args$robust.se <- TRUE
  }
  
  out <- do.call(modsem::modsem, args)
  
  # entire parameter table
  out$parTable
}

##==============================================================================
## 2. Helper function to extract specific parameters from stored tables
##==============================================================================
extract_eta3_parameters <- function(table, method_type) {
  if(is.null(table)) return(NULL)
  
  if(method_type == "dblcent") {
    # for dblcent: use coefParTable structure
    rows <- table[table$lhs == "eta3" & table$op == "~", ]
    
    params <- if(nrow(rows) == 3) {
      c("eta1", "eta2", "eta1:eta2")
    } else {
      c("eta1", "eta2", "eta1:eta2", "eta1:eta1", "eta2:eta2")
    }
    
    RESULTS <- list(
      "Estimates" = setNames(rows$est, params),
      "Standard Errors" = setNames(rows$se, params),
      "P-values" = setNames(rows$pvalue, params),
      "CI_lower" = setNames(rows$ci.lower, params),
      "CI_upper" = setNames(rows$ci.upper, params)
    )
    
  } else if(method_type == "sam") {
    # for SAM: use parameterEstimates structure
    coefs <- table[table$lhs == "eta3" & table$op == "~", ]
    
    params <- if(nrow(coefs) == 3) {
      c("eta1", "eta2", "eta1:eta2")
    } else {
      c("eta1", "eta2", "eta1:eta2", "eta1:eta1", "eta2:eta2")
    }
    
    rows <- table$lhs == "eta3" & table$op == "~" & table$rhs %in% params
    
    ests <- setNames(table$est[rows], table$rhs[rows])
    ses  <- setNames(table$se[rows], table$rhs[rows])
    pval <- setNames(table$pvalue[rows], table$rhs[rows])
    ci_lower <- setNames(table$ci.lower[rows], table$rhs[rows])
    ci_upper <- setNames(table$ci.upper[rows], table$rhs[rows])
    
    # order according to `params`
    RESULTS <- list(
      "Estimates"        = ests[params],
      "Standard Errors"  = ses[params],
      "P-values"         = pval[params],
      "CI_lower"         = ci_lower[params],
      "CI_upper"         = ci_upper[params]
    )
    
  } else if(method_type %in% c("lms", "qml")) {
    # for LMS/QML: use parTable structure
    rows <- table[table$lhs == "eta3" & table$op == "~", ]
    
    # 5 parameters; we need to reorder because eta1:eta1 before eta1:eta2
    if(nrow(rows) == 5) {
      # indices of the parameters we want to swap
      eta1eta1_idx <- which(rows$rhs == "eta1:eta1")
      eta1eta2_idx <- which(rows$rhs == "eta1:eta2")
      
      # reordered index vector
      idx <- 1:nrow(rows)
      idx[eta1eta1_idx] <- eta1eta2_idx
      idx[eta1eta2_idx] <- eta1eta1_idx
      
      # swapped parameter order
      RESULTS <- list(
        "Estimates" = setNames(rows$est[idx], rows$rhs[idx]),
        "Standard Errors" = setNames(rows$std.error[idx], rows$rhs[idx]),
        "P-values" = setNames(rows$p.value[idx], rows$rhs[idx]),
        "CI_lower" = setNames(rows$ci.lower[idx], rows$rhs[idx]),
        "CI_upper" = setNames(rows$ci.upper[idx], rows$rhs[idx])
      )
    } else {
      # 3 parameters, keep as it was before
      RESULTS <- list(
        "Estimates" = setNames(rows$est, rows$rhs),
        "Standard Errors" = setNames(rows$std.error, rows$rhs),
        "P-values" = setNames(rows$p.value, rows$rhs),
        "CI_lower" = setNames(rows$ci.lower, rows$rhs),
        "CI_upper" = setNames(rows$ci.upper, rows$rhs)
      )
    }
  }
  
  RESULTS
}

