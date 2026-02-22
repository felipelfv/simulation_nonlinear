################################################################################
# File:         Methods.R
# Description:  Estimation functions for nonlinear SEM. Wraps modsem for UPI,
#               LMS, and QML, and lavaan for LSAM. All functions return the
#               full parameter table.
# Dependencies: lavaan, modsem
# Used by:      Simulation(1).R, Simulation(2).R

############################### Function Documentation #########################

#' method_upi: Product Indicator Approach (UPI) with Double Mean Centering
#' @param Data          Data.frame. The dataset containing observed variables.
#' @param model.fit     Character. Model syntax in lavaan format specifying the SEM with interactions.
#' @param robust.se     Logical. Whether to use robust (Huber-White) standard errors (default = FALSE).
#' @param match         Character or NULL. Matching specification for product indicators (default = NULL).
#' @param bounds        Logical. Whether to apply bounds to parameter estimates (default = FALSE).
#'
#' method_lsam: Local Structural-After-Measurement (LSAM) Approach
#' @param Data          Data.frame. The dataset containing observed variables.
#' @param estimator     Character. Estimator to use (default = "ML").
#' @param joint         Logical. Whether to use joint estimation (default = TRUE).
#' @param model.fit     Character. Model syntax in lavaan format specifying the SEM with interactions.
#' @param mm.list       List or NULL. Measurement model blocks. If NULL and joint = TRUE,
#'                      all latent variables are placed in one block (default = NULL).
#'
#' method_analytic: Distribution Analytic Approaches (LMS or QML)
#' @param Data          Data.frame. The dataset containing observed variables.
#' @param model.fit     Character. Model syntax in lavaan format specifying the SEM with interactions.
#' @param standardized  Logical. Whether to return standardized estimates (default = FALSE).
#' @param method        Character. Distribution analytic method to use: "lms" (default) or "qml".
#' @param robust.se     Logical. Whether to compute robust standard errors (default = FALSE).

################################ Functions #####################################

# 1.1. Extended Product Indicator Approach with double mean centering (UPI)
method_upi <- function(Data = NULL, model.fit = NULL,
                       robust.se = FALSE, match = NULL, bounds = FALSE) {
  
  args <- list(
    model.syntax = model.fit,
    data = Data,
    method = "dblcent",
    bounds = bounds,
    match = match
  )
  
  if (robust.se) {
    args$se <- "robust.huber.white"
  }
  
  out <- do.call(modsem::modsem_pi, args) 
  out$coefParTable
}

# 1.2 Local Structural-After-Measurement (LSAM) Approach 
method_lsam <- function(Data = NULL, estimator = "ML",
                        joint = TRUE,
                        model.fit = NULL,
                        mm.list = NULL) {
  
  # joint=TRUE and mm.list not provided, put all LVs in one block
  if(joint && is.null(mm.list)) {
    if(is.character(model.fit)) {
      pt <- lavaan::lavaanify(model.fit)
    } else {
      pt <- lavaan::parTable(model.fit)
    }
    
    lv_names <- unique(pt$lhs[pt$op == "=~"])
    mm.list <- list(lv_names)
  }
  
  out <- lavaan::sam(model.fit, data = Data, se = "local",
                     mm.args = list(estimator = estimator),
                     mm.list = mm.list) # note: the method for lsam will be local once interactions are detected
  
  parameterEstimates(out, remove.step1 = FALSE)
}

# 1.3 Distribution Analytic Approaches (LMS and QML) 
method_analytic <- function(Data = NULL, model.fit = NULL, 
                            standardized = FALSE, method = "lms", robust.se = FALSE) {
  # for auto.split.syntax, the default is therefore TRUE for the QML approach
  # build argument list
  # default is method lms, but change method to "qml" for qml
  args <- list(
    model.syntax = model.fit, 
    data = Data, 
    method = method,
    robust.se = robust.se 
  )
  
  out <- do.call(modsem::modsem, args)
  out$parTable
}
