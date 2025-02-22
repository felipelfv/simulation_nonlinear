#### 1. General Session Info ####

#### 2 Structural-After-Measurement (SAM) Approach ####

library(lavaan)

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
