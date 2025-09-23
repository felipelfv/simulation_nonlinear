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
  # for auto.split.syntax, the default is therefore TRUE for the QML approach
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

