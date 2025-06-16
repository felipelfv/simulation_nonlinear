##==============================================================================
## General Information 
##
## Script name: Methods(1).R
##
## Purpose of script: functions used for the four different approaches. 
## 
#' @param Data Data frame.
#' @param model.fit Model used for the estimation.
#' 
#' 
#' 
##==============================================================================

##==============================================================================
## 1. Methods: UCA, SAM, LMS, and QML 
##==============================================================================

#### 1.1. Product Indicator Approaches (UCA) ####

# Note: On page 186, Brandt et al. (2014) claim that for the interaction model
# the term unconstrained approach (as in Marsh et al., 2004, 2006) would be 
# more adequate, but for simplicity they refer to the approach as ExUC.

method_uca <- function(Data = NULL, model.fit = NULL) {
  out <- modsem::modsem(model.syntax = model.fit, data = Data, method = "uca")
  rows <- out$coefParTable[out$coefParTable$lhs == "eta3" & out$coefParTable$op == "~", ]
  
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
  
  RESULTS
}

#### 1.2 Structural-After-Measurement (SAM) Approach ####

method_sam <- function(Data = NULL, estimator = "ML",
                       joint = FALSE, add.attr = FALSE, 
                       model.fit = NULL,
                       mm.list = NULL) {
  out <- lavaan::sam(model.fit, data = Data, se = "local",
                     mm.args = list(estimator = estimator),
                     mm.list = mm.list)
  
  # parameter estimates 
  est_table <- parameterEstimates(out)
  coefs <- est_table[est_table$lhs == "eta3" & est_table$op == "~", ]
  
  params <- if(nrow(coefs) == 3) {
    c("eta1", "eta2", "eta1:eta2")
  } else {
    c("eta1", "eta2", "eta1:eta2", "eta1:eta1", "eta2:eta2")
  }
  
  rows <- est_table$lhs == "eta3" & est_table$op == "~" & est_table$rhs %in% params
  
  ests <- setNames(est_table$est[rows], est_table$rhs[rows])
  ses  <- setNames(est_table$se[rows], est_table$rhs[rows])
  pval <- setNames(est_table$pvalue[rows], est_table$rhs[rows])
  ci_lower <- setNames(est_table$ci.lower[rows], est_table$rhs[rows])
  ci_upper <- setNames(est_table$ci.upper[rows], est_table$rhs[rows])
  
  # order according to `params`
  RESULTS <- list(
    "Estimates"        = ests[params],
    "Standard Errors"  = ses[params],
    "P-values"         = pval[params],
    "CI_lower"         = ci_lower[params],
    "CI_upper"         = ci_upper[params]
  )
  
  RESULTS
}


#### 1.3 Analytic Approaches (LMS and QML) ####

# Note: both approaches are estimated using the same function. 
# The argument "method" used in the Simulation script that calls each.

method_analytic <- function(Data = NULL, model.fit = NULL, 
                            standardized = FALSE, method = "lms") {
  out <- modsem::modsem(model.syntax = model.fit, data = Data, method = method)
  rows <- out$parTable[out$parTable$lhs == "eta3" & out$parTable$op == "~", ]
  
  # 5 parameters; we need to reorder because eta1:eta1 before eta1:eta2
  if(nrow(rows) == 5) {
    # indices of the parameters we want to swap
    eta1eta1_idx <- which(rows$rhs == "eta1:eta1")
    eta1eta2_idx <- which(rows$rhs == "eta1:eta2")
    
    # reordered index vector
    idx <- 1:nrow(rows)
    idx[eta1eta1_idx] <- eta1eta2_idx
    idx[eta1eta2_idx] <- eta1eta1_idx
    
    # Swapped parameter order
    # "eta1", "eta2", "eta1:eta2", "eta1:eta1", "eta2:eta2"
    RESULTS <- list(
      "Estimates" = setNames(rows$est[idx], rows$rhs[idx]),
      "Standard Errors" = setNames(rows$std.error[idx], rows$rhs[idx]),
      "P-values" = setNames(rows$p.value[idx], rows$rhs[idx]),
      "CI_lower" = setNames(rows$ci.lower[idx], rows$rhs[idx]),
      "CI_upper" = setNames(rows$ci.upper[idx], rows$rhs[idx])
    )
  } else {
    # 3 parameters, keep as it was before:
    # "eta1", "eta2", "eta1:eta2"
    RESULTS <- list(
      "Estimates" = setNames(rows$est, rows$rhs),
      "Standard Errors" = setNames(rows$std.error, rows$rhs),
      "P-values" = setNames(rows$p.value, rows$rhs),
      "CI_lower" = setNames(rows$ci.lower, rows$rhs),
      "CI_upper" = setNames(rows$ci.upper, rows$rhs)
    )
  }
  
  RESULTS
}

##### 2.3.1 NSEMM Approach #####

# Note: as of (26/02/25), not much advantage over LMS/QML given the current 
# models being estimated. 
