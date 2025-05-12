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
  rows <- out$coefParTable[out$coefParTable$lhs == "eta4" & out$coefParTable$op == "~", ]
  
  params <- c("eta1", "eta2", "eta3", "eta2:eta3", "eta1:eta3", "eta6:eta2", "eta5:eta6", 
              "eta5:eta5", "eta1:eta2") # Works according to the model order
  
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
  
  all_coefs <- names(coef(out))[grepl("eta4~", names(coef(out)))]
  all_coefs <- all_coefs[!grepl("eta4~~|eta4~1", all_coefs)] # just for now
  
  params <- c("eta1", "eta2", "eta3", "eta2:eta3", "eta1:eta3", "eta6:eta2", "eta5:eta6", 
              "eta5:eta5", "eta1:eta2") # Works well with the model order
  
  coefs <- paste0("eta4~", params)
  RESULTS <- list("Estimates" = setNames(coef(out)[coefs], params))
  RESULTS
}


#### 2.3 Analytic Approaches (LMS and QML) ####

# Note: both approaches are estimated using the same function. 
# The argument "method" used in the Simulation script that calls each.

method_analytic <- function(Data = NULL, model.fit = NULL, 
                            standardized = FALSE, method = "lms") {
  out <- modsem::modsem(model.syntax = model.fit, data = Data, method = method)
  rows <- out$parTable[out$parTable$lhs == "eta4" & out$parTable$op == "~", ]
  
  # Need to change for obtaining the estimate values in a more accurate way. Right now it doesnt work well (!)
  # Order in QML: "eta2", "eta1", "eta3", "eta2:eta3", "eta2:eta6", "eta1:eta2", "eta1:eta3", "eta5:eta5", "eta5:eta6"
  # Order in LMS: (...)
  params <- c("eta1", "eta2", "eta3", "eta2:eta3", "eta1:eta3", "eta6:eta2", "eta5:eta6", 
              "eta5:eta5", "eta1:eta2") 
  
  RESULTS <- list(
    "Estimates" = setNames(rows$est, params),
    "Standard Errors" = setNames(rows$std.error, params),
    "P-values" = setNames(rows$p.value, params)
  )
  
  RESULTS
}

##### 2.3.1 NSEMM Approach #####

# Note: as of (26/02/25), not much advantage over LMS/QML given the current 
# models being estimated. 
