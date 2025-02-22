#### 1. General Session Info ####

#### 2. Product Indicator Approaches (UCA) ####

library(modsem)

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


