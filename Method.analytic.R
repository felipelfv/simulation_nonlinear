#### 1. General Session Info ####

#### 2. Analytic Approaches (LMS and QML) ####

library(modsem)

method_analytic <- function(Data = NULL, model.fit = NULL, 
                            standardized = FALSE, method = "lms") {
  out <- modsem::modsem(model.syntax = model.fit, data = Data, method = method)
  rows <- out$parTable[out$parTable$lhs == "eta3" & out$parTable$op == "~", ]
  
  params <- if(nrow(rows) == 3) {
    c("eta1", "eta2", "eta1:eta2")
  } else {
    c("eta1", "eta2", "eta1:eta2", "eta1:eta1", "eta2:eta2")
  }
  
  RESULTS <- list(
    "Estimates" = setNames(rows$est, params),
    "Standard Errors" = setNames(rows$std.error, params),
    "P-values" = setNames(rows$p.value, params)
  )
  
  RESULTS
}

#### 2.1 NSEMM Approach ####

