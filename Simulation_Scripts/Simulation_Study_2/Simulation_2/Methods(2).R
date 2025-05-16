#### 2.1. Product Indicator Approaches (UCA) ####
method_uca <- function(Data = NULL, model.fit = NULL) {
  out <- modsem::modsem(model.syntax = model.fit, data = Data, method = "uca")
  rows <- out$coefParTable[out$coefParTable$lhs == "eta4" & out$coefParTable$op == "~", ]
  
  # parameter names based on how many parameters are in the model
  params <- switch(as.character(nrow(rows)),
                   "4" = c("eta1", "eta2", "eta3", "eta2:eta3"),
                   "5" = c("eta1", "eta2", "eta3", "eta2:eta3", "eta1:eta3"),
                   "6" = c("eta1", "eta2", "eta3", "eta2:eta3", "eta1:eta3", "eta6:eta2"),
                   "7" = c("eta1", "eta2", "eta3", "eta2:eta3", "eta1:eta3", "eta6:eta2", "eta5:eta6"),
                   "8" = c("eta1", "eta2", "eta3", "eta2:eta3", "eta1:eta3", "eta6:eta2", "eta5:eta6", "eta5:eta5"),
                   "9" = c("eta1", "eta2", "eta3", "eta2:eta3", "eta1:eta3", "eta6:eta2", "eta5:eta6", "eta5:eta5", "eta1:eta2"),
                   "10" = c("eta1", "eta2", "eta3", "eta2:eta3", "eta1:eta3", "eta6:eta2", "eta5:eta6", "eta5:eta5", "eta1:eta2", "eta6:eta6"),
                   rows$rhs) # use actual names if not in predefined set
  
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
  all_coefs <- all_coefs[!grepl("eta4~~|eta4~1", all_coefs)]
  
  # parameter names based on how many parameters are in the model
  params <- switch(as.character(length(all_coefs)),
                   "4" = c("eta1", "eta2", "eta3", "eta2:eta3"),
                   "5" = c("eta1", "eta2", "eta3", "eta2:eta3", "eta1:eta3"),
                   "6" = c("eta1", "eta2", "eta3", "eta2:eta3", "eta1:eta3", "eta6:eta2"),
                   "7" = c("eta1", "eta2", "eta3", "eta2:eta3", "eta1:eta3", "eta6:eta2", "eta5:eta6"),
                   "8" = c("eta1", "eta2", "eta3", "eta2:eta3", "eta1:eta3", "eta6:eta2", "eta5:eta6", "eta5:eta5"),
                   "9" = c("eta1", "eta2", "eta3", "eta2:eta3", "eta1:eta3", "eta6:eta2", "eta5:eta6", "eta5:eta5", "eta1:eta2"),
                   "10" = c("eta1", "eta2", "eta3", "eta2:eta3", "eta1:eta3", "eta6:eta2", "eta5:eta6", "eta5:eta5", "eta1:eta2", "eta6:eta6"),
                   sub("eta4~", "", all_coefs)) # fallback as above
  
  coefs <- paste0("eta4~", params)
  RESULTS <- list("Estimates" = setNames(coef(out)[coefs], params))
  RESULTS
}

#### 2.3 Analytic Approaches (LMS and QML) ####
method_analytic <- function(Data = NULL, model.fit = NULL, 
                            standardized = FALSE, method = "lms") {
  out <- modsem::modsem(model.syntax = model.fit, data = Data, method = method)
  rows <- out$parTable[out$parTable$lhs == "eta4" & out$parTable$op == "~", ]
  
  # parameter orders for analytical methods (same for QML and LMS)
  analytic_params <- switch(as.character(nrow(rows)),
                            "4" = c("eta2", "eta1", "eta3", "eta2:eta3"),
                            "5" = c("eta2", "eta1", "eta3", "eta2:eta3", "eta1:eta3"),
                            "6" = c("eta2", "eta1", "eta3", "eta2:eta3", "eta1:eta2", "eta1:eta3"),
                            "7" = c("eta2", "eta1", "eta3", "eta2:eta3", "eta1:eta2", "eta1:eta3", "eta6:eta2"),
                            "8" = c("eta2", "eta1", "eta3", "eta2:eta3", "eta1:eta2", "eta1:eta3", "eta6:eta2", "eta6:eta6"),
                            "9" = c("eta2", "eta1", "eta3", "eta2:eta3", "eta1:eta2", "eta1:eta3", "eta6:eta2", "eta6:eta6", "eta5:eta6"),
                            "10" = c("eta2", "eta1", "eta3", "eta2:eta3", "eta1:eta2", "eta1:eta3", "eta6:eta2", "eta6:eta6", "eta5:eta6", "eta5:eta5"))
  
  # desired output order (matching other methods)
  desired_output_order <- switch(as.character(nrow(rows)),
                                 "4" = c("eta1", "eta2", "eta3", "eta2:eta3"),
                                 "5" = c("eta1", "eta2", "eta3", "eta2:eta3", "eta1:eta3"),
                                 "6" = c("eta1", "eta2", "eta3", "eta2:eta3", "eta1:eta3", "eta6:eta2"),
                                 "7" = c("eta1", "eta2", "eta3", "eta2:eta3", "eta1:eta3", "eta6:eta2", "eta5:eta6"),
                                 "8" = c("eta1", "eta2", "eta3", "eta2:eta3", "eta1:eta3", "eta6:eta2", "eta5:eta6", "eta5:eta5"),
                                 "9" = c("eta1", "eta2", "eta3", "eta2:eta3", "eta1:eta3", "eta6:eta2", "eta5:eta6", "eta5:eta5", "eta1:eta2"),
                                 "10" = c("eta1", "eta2", "eta3", "eta2:eta3", "eta1:eta3", "eta6:eta2", "eta5:eta6", "eta5:eta5", "eta1:eta2", "eta6:eta6"))
  
  # if we have parameter orders defined for this model size
  if(is.null(analytic_params)) {
    # fallback: just use actual parameter names
    warning("No predefined parameter order for model with ", nrow(rows), " parameters.")
    RESULTS <- list(
      "Estimates" = setNames(rows$est, rows$rhs),
      "Standard Errors" = setNames(rows$std.error, rows$rhs),
      "P-values" = setNames(rows$p.value, rows$rhs)
    )
    RESULTS
  }
  
  # For debugging - check if the parameter sets match
  if(length(rows$rhs) != length(analytic_params) || !all(sort(rows$rhs) == sort(analytic_params))) {
    warning("Parameter names in model output don't match expected parameters for ", method, 
            ". Expected: ", paste(analytic_params, collapse=", "), 
            ". Got: ", paste(rows$rhs, collapse=", "))
    RESULTS <- list(
      "Estimates" = setNames(rows$est, rows$rhs),
      "Standard Errors" = setNames(rows$std.error, rows$rhs),
      "P-values" = setNames(rows$p.value, rows$rhs)
    )
    RESULTS
  }
  
  # mapping index to reorder from actual order to desired order
  reorder_idx <- match(desired_output_order, rows$rhs)
  
  # reordering
  RESULTS <- list(
    "Estimates" = setNames(rows$est[reorder_idx], desired_output_order),
    "Standard Errors" = setNames(rows$std.error[reorder_idx], desired_output_order),
    "P-values" = setNames(rows$p.value[reorder_idx], desired_output_order)
  )
  
  RESULTS
}