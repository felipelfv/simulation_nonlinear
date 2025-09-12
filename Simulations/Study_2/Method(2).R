method_sam <- function(Data = NULL, estimator = "ML",
                       joint = FALSE, add.attr = FALSE, 
                       model.fit = NULL) {
  out <- lavaan::sam(model.fit, data = Data, se = "none",
                     mm.args = list(estimator = estimator))
  
  # parameter estimates 
  est_table <- parameterEstimates(out)
  
  # Get all regression parameters (all paths with operator "~")
  reg_rows <- est_table$op == "~"
  
  # Extract values
  # Use path notation (lhs~rhs) for parameter names to distinguish between equations
  param_names <- paste0(est_table$lhs[reg_rows], "~", est_table$rhs[reg_rows])
  ests <- setNames(est_table$est[reg_rows], param_names)
  #ses  <- setNames(est_table$se[reg_rows], param_names)
  #pval <- setNames(est_table$pvalue[reg_rows], param_names)
  
  # Return results
  RESULTS <- list(
    "Estimates"        = ests
    #"Standard Errors"  = ses,
    #"P-values"         = pval
  )
  
  RESULTS
}