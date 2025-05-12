#### 1. General Session Info ####

library(modsem) # For QML, LMS, and UCA
library(lavaan) # For SAM

#### 2. Methods: UCA, SAM, LMS, and QML ####

#### 2.1. Product Indicator Approaches (UCA) ####

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
    "P-values" = setNames(rows$pvalue, params)
  )
  
  RESULTS
}

#### 2.2 Structural-After-Measurement (SAM) Approach ####

method_sam <- function(Data = NULL, estimator = "ML",
                       joint = FALSE, add.attr = FALSE, 
                       model.fit = NULL) {
  out <- lavaan::sam(model.fit, data = Data, se = "none",
                     mm.args = list(estimator = estimator))
  
  all_coefs <- names(coef(out))[grepl("eta3~", names(coef(out)))]
  all_coefs <- all_coefs[!grepl("eta3~~|eta3~1", all_coefs)] # just for now
  
  params <- if(length(all_coefs) == 3) {
    c("eta1", "eta2", "eta1:eta2")
  } else {
    c("eta1", "eta2", "eta1:eta2", "eta1:eta1", "eta2:eta2")
  }
  
  coefs <- paste0("eta3~", params)
  RESULTS <- list("Estimates" = setNames(coef(out)[coefs], params))
  RESULTS
}


#### 2.3 Analytic Approaches (LMS and QML) ####

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
    RESULTS <- list(
      "Estimates" = setNames(rows$est[idx], rows$rhs[idx]),
      "Standard Errors" = setNames(rows$std.error[idx], rows$rhs[idx]),
      "P-values" = setNames(rows$p.value[idx], rows$rhs[idx])
    )
  } else {
    # 3 parameters, keep as it was before:
    RESULTS <- list(
      "Estimates" = setNames(rows$est, rows$rhs),
      "Standard Errors" = setNames(rows$std.error, rows$rhs),
      "P-values" = setNames(rows$p.value, rows$rhs)
    )
  }
  
  RESULTS
}

##### 2.3.1 NSEMM Approach #####

# Note: as of (26/02/25), not much advantage over LMS/QML given the current 
# models being estimated. 
