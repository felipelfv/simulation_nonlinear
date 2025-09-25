############################ 1. General Information ############################
# See README file for more information concerning this file. 
# This file contains the code used to estimate the different approaches. We use
# modsem for all except the SAM approach. Thus, we have the DBLCENT, LMS, and 
# QML estimated through modsem and the SAM approach with lavaan. 
# All functions return the full parameter table.
# Importantly, the helper function is used to extract only the relevant 
# estimates. In our case, the structural model. 

############################### 2. Function Documentation #######################
#' method_dblcent: Product Indicator Approach (UPI) with Double Mean Centering
#' 
#' @param Data          Data.frame. The dataset containing observed variables.
#' @param model.fit     Character. Model syntax in lavaan format specifying the SEM with interactions.
#' @param robust.se     Logical. Whether to use robust (Huber-White) standard errors (default = FALSE).
#' @param match         Character or NULL. Matching specification for product indicators (default = NULL).
#' @param bounds        Logical. Whether to apply bounds to parameter estimates (default = FALSE).

#' method_sam: Local Structural-After-Measurement (LSAM) Approach
#' 
#' @param Data          Data.frame. The dataset containing observed variables.
#' @param estimator     Character. Estimator to use (default = "ML"), switches to "MLR" if robust.se = TRUE.
#' @param joint         Logical. Whether to use joint estimation (default = FALSE).
#' @param add.attr      Logical. Whether to add additional attributes to output (default = FALSE).
#' @param model.fit     Character. Model syntax in lavaan format specifying the SEM with interactions.
#' @param mm.list       List or NULL. Additional arguments for the measurement model step (default = NULL).

#' method_analytic: Distribution Analytic Approaches (LMS or QML)
#' 
#' @param Data          Data.frame. The dataset containing observed variables.
#' @param model.fit     Character. Model syntax in lavaan format specifying the SEM with interactions.
#' @param standardized  Logical. Whether to return standardized estimates (default = FALSE).
#' @param method        Character. Distribution analytic method to use: "lms" (default) or "qml".
#' @param robust.se     Logical. Whether to compute robust standard errors (default = FALSE).

############################### 3. Functions ###################################

#### 3.1. Product Indicator Approach with double mean centering (DBLCENT) ####
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

#### 3.2 Structural-After-Measurement (SAM) Approach ####
method_sam <- function(Data = NULL, estimator = "ML",
                       joint = FALSE, add.attr = FALSE, 
                       model.fit = NULL,
                       mm.list = NULL) {
  
  out <- lavaan::sam(model.fit, data = Data, se = "local",
                     mm.args = list(estimator = estimator),
                     mm.list = mm.list)
  
  # parameter estimates table
  parameterEstimates(out, remove.step1 = FALSE)
}

#### 3.3 Distribution Analytic Approaches (LMS and QML) ####
method_analytic <- function(Data = NULL, model.fit = NULL, 
                            standardized = FALSE, method = "lms", robust.se = FALSE) {
  # for auto.split.syntax, the default is therefore TRUE for the QML approach
  # Build argument list
  # default is method lms, but change method to "qml" for qml
  args <- list(model.syntax = model.fit, data = Data, method = method)
  
  # robust SE if requested
  if (robust.se) {
    args$robust.se <- TRUE
  }
  
  out <- do.call(modsem::modsem, args)
  
  # entire parameter table
  out$parTable
}

