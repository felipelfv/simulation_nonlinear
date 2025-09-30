############################ 1. General Information ############################
# See README file for more information concerning this file. 

# This file contains the code necessary to calculate the performance results 
# reported in the manuscript. The final dataset contains all the necessary
# metrics to plot. These are also used in the tables reported. 

############################### 2. Documentation ################################

#' extract_eta3_parameters: Extract Structural Parameters from Method Output
#' 
#' @param table         Data.frame. Parameter table from estimation method output
#' @param method_type   Character. Type of method: "dblcent", "sam", "lms", or "qml"
#' 
#' @return List containing:
#'   - Estimates:       Named vector of parameter estimates
#'   - Standard Errors: Named vector of standard errors
#'   - P-values:        Named vector of p-values
#'   - CI_lower:        Named vector of confidence interval lower bounds
#'   - CI_upper:        Named vector of confidence interval upper bounds
#'
#' ExtractConvergenceOutliers: Process Convergence and Identify Outliers
#' 
#' @param all_results              List. Complete simulation results from all conditions
#' @param parameters_of_interest   Character vector. Parameters to extract (default = c("eta1","eta2","eta1:eta2","eta1:eta1","eta2:eta2"))
#' @param remove_outliers          Logical. Whether to identify and exclude outliers (default = TRUE)
#' @param outlier_threshold        Numeric. IQR multiplier for outlier detection (default = 3)
#' @param min_reps                 Integer. Minimum replications required per condition (default = 10)
#' @param exclude_warnings         Logical. Whether to exclude iterations with warnings (default = FALSE)
#' 
#' @return List containing:
#'   - convergence_outliers_summary: Tibble with convergence/outlier statistics per condition/method/parameter
#'   - convergence_outliers_details: List with detailed iteration-level information
#'   - filtered_data:               List with clean data ready for performance calculations
#' 
#' @details Processing steps:
#'   1. Check convergence for each method independently (non-missing, finite values)
#'   2. Track warnings in converged iterations
#'   3. Identify valid replications (converged across all methods)
#'   4. Optionally exclude iterations with warnings
#'   5. Detect outliers using IQR method (Q1 - threshold*IQR, Q3 + threshold*IQR)
#'   6. Apply global synchronization (exclude if outlier in any method)
#'   7. Extract filtered data for performance metrics
#' 
#' @note Convergence criteria:
#'   - All required parameters present
#'   - No NA, NaN, or Inf values in estimates, SEs, p-values, or CIs
#'   - Method-specific convergence tracked independently
#'   - Global convergence requires all methods converged
#' 
#' @note Summary statistics tracked:
#'   - N_Total:                     Total replications attempted
#'   - N_Converged_This_Method:     Converged for this specific method
#'   - N_Warnings_This_Method:      Converged iterations with warnings for this method
#'   - Convergence_Rate_This_Method: Percentage converged for this method
#'   - N_After_Global_Convergence:  Converged across all methods
#'   - Global_Convergence_Rate:     Percentage converged across all methods
#'   - N_Excluded_By_Warnings:      Number excluded due to warnings (if exclude_warnings = TRUE)
#'   - N_Outliers_This_Method:      Outliers detected in this method
#'   - Outlier_Rate_This_Method:    Percentage of outliers among converged
#'   - N_After_Global_Outlier:      Valid after global outlier removal
#'   - Global_Outlier_Rate:         Percentage valid after outlier removal
#'   - N_Final:                     Final sample size for analysis

#' CalculatePerformanceMetrics: Calculate All Performance Metrics
#' 
#' @param filtered_data                 List. Filtered data from ExtractConvergenceOutliers
#' @param convergence_outliers_summary  Tibble. Summary statistics from ExtractConvergenceOutliers
#' @param alpha                         Numeric. Significance level for tests (default = 0.05)
#' 
#' @return Tibble containing all performance metrics per condition/method/parameter:
#' 
#' @details Mean-based metrics:
#'   - MeanEstimate:                Average of parameter estimates
#'   - Bias_Mean:                   Mean estimate - true value
#'   - RelativeBias_Mean:            Bias / true value - 1
#'   - PercentRelativeBias_Mean:     Relative bias * 100
#'   - Variance:                     Variance of estimates
#'   - SD:                          Standard deviation of estimates
#'   - MSE_Mean:                    Mean squared error
#'   - RMSE_Mean:                   Root mean squared error
#'   
#' @details Median-based metrics:
#'   - MedianEstimate:              Median of parameter estimates
#'   - Bias_Median:                 Median estimate - true value  
#'   - RelativeBias_Median:         Median bias / true value
#'   - PercentRelativeBias_Median:  Relative median bias * 100
#'   - MAD:                         Median absolute deviation
#'   - RMSE_Median:                 sqrt(Bias_Median^2 + MAD^2)
#'   
#' @details Standard error metrics:
#'   - MeanSE:                      Average of standard errors
#'   - SE_SD_Ratio:                 Mean SE / SD of estimates
#'   
#' @details Coverage and confidence interval metrics:
#'   - CoverageRate:                Percentage of CIs containing true value
#'   - CI_Width:                    Average width of confidence intervals
#'   
#' @details Hypothesis testing metrics:
#'   - RejectionRate:               Percentage of p-values < alpha
#'   - TypeI_Error:                 Rejection rate when true value equal to 0
#'   - Power:                       Rejection rate when true value not equal to 0
#'   
#' @details Monte Carlo standard errors (MCSE):
#'   - Bias_Mean_MCSE, RelativeBias_MCSE, Variance_MCSE, SD_MCSE
#'   - MSE_Mean_MCSE, RMSE_Mean_MCSE, Relative_MSE_MCSE, Relative_RMSE_MCSE
#'   - CoverageRate_MCSE, CI_Width_MCSE, RejectionRate_MCSE
#'   - TypeI_Error_MCSE, Power_MCSE
#'   
#' @details Sample size tracking (from convergence_outliers_summary):
#'   - All N_* and convergence/outlier statistics
#'   - N_Excluded_* breakdown by reason
#'   - Percent_Total_Excluded

############################### 3. Functions ####################################

## helper function to extract specific parameters from stored tables
extract_eta3_parameters <- function(table, method_type) {
  if(is.null(table)) return(NULL)
  
  if(method_type == "dblcent") {
    # for dblcent: use coefParTable structure
    rows <- table[table$lhs == "eta3" & table$op == "~", ]
    
    params <- if(nrow(rows) == 3) {
      c("eta1", "eta2", "eta1:eta2")
    } else {
      c("eta1", "eta2", "eta1:eta2", "eta1:eta1", "eta2:eta2")
    }
    
    RESULTS <- list(
      "Estimates" = setNames(rows$est, params),
      "Standard Errors" = setNames(rows$se, params),
      "P-values" = setNames(rows$pvalue, params),
      "CI_lower" = setNames(rows$ci.lower, params),
      "CI_upper" = setNames(rows$ci.upper, params)
    )
    
  } else if(method_type == "sam") {
    # for SAM: use parameterEstimates structure
    coefs <- table[table$lhs == "eta3" & table$op == "~", ]
    
    params <- if(nrow(coefs) == 3) {
      c("eta1", "eta2", "eta1:eta2")
    } else {
      c("eta1", "eta2", "eta1:eta2", "eta1:eta1", "eta2:eta2")
    }
    
    rows <- table$lhs == "eta3" & table$op == "~" & table$rhs %in% params
    
    ests <- setNames(table$est[rows], table$rhs[rows])
    ses  <- setNames(table$se[rows], table$rhs[rows])
    pval <- setNames(table$pvalue[rows], table$rhs[rows])
    ci_lower <- setNames(table$ci.lower[rows], table$rhs[rows])
    ci_upper <- setNames(table$ci.upper[rows], table$rhs[rows])
    
    # order according to `params`
    RESULTS <- list(
      "Estimates"        = ests[params],
      "Standard Errors"  = ses[params],
      "P-values"         = pval[params],
      "CI_lower"         = ci_lower[params],
      "CI_upper"         = ci_upper[params]
    )
    
  } else if(method_type %in% c("lms", "qml")) {
    # for LMS/QML: use parTable structure
    rows <- table[table$lhs == "eta3" & table$op == "~", ]
    
    # 5 parameters; we need to reorder because eta1:eta1 before eta1:eta2
    if(nrow(rows) == 5) {
      # indices of the parameters we want to swap
      eta1eta1_idx <- which(rows$rhs == "eta1:eta1")
      eta1eta2_idx <- which(rows$rhs == "eta1:eta2")
      
      # reordered index vector
      idx <- 1:nrow(rows)
      idx[eta1eta1_idx] <- eta1eta2_idx
      idx[eta1eta2_idx] <- eta1eta1_idx
      
      # swapped parameter order
      RESULTS <- list(
        "Estimates" = setNames(rows$est[idx], rows$rhs[idx]),
        "Standard Errors" = setNames(rows$std.error[idx], rows$rhs[idx]),
        "P-values" = setNames(rows$p.value[idx], rows$rhs[idx]),
        "CI_lower" = setNames(rows$ci.lower[idx], rows$rhs[idx]),
        "CI_upper" = setNames(rows$ci.upper[idx], rows$rhs[idx])
      )
    } else {
      # 3 parameters, keep as it was before
      RESULTS <- list(
        "Estimates" = setNames(rows$est, rows$rhs),
        "Standard Errors" = setNames(rows$std.error, rows$rhs),
        "P-values" = setNames(rows$p.value, rows$rhs),
        "CI_lower" = setNames(rows$ci.lower, rows$rhs),
        "CI_upper" = setNames(rows$ci.upper, rows$rhs)
      )
    }
  }
  
  RESULTS
}

## Function 1: Convergence and Outliers
ExtractConvergenceOutliers <- function(all_results,
                                       parameters_of_interest = c("eta1","eta2","eta1:eta2","eta1:eta1","eta2:eta2"),
                                       remove_outliers = TRUE,
                                       outlier_threshold = 3,
                                       min_reps = 10,
                                       exclude_warnings = FALSE) {
  
  options(scipen = 999)
  
  methods <- c("sam", "lms", "qml", "dblcent")
  convergence_outliers_summary <- dplyr::tibble()
  convergence_outliers_details <- list()
  filtered_data <- list()
  param_index <- c("eta1"=1, "eta2"=2, "eta1:eta2"=3, "eta1:eta1"=4, "eta2:eta2"=5)
  
  for (i in seq_along(all_results)) {
    condition <- all_results[[i]]$condition
    res_list  <- all_results[[i]]$results
    
    # Get true values
    true_values <- all_results[[i]]$true_parameters
    if (is.null(true_values)) {
      if (condition$Model_Type == "linear") {
        true_values <- c(0.316, 0.316, 0, 0, 0)
      } else {  # Must be "full"
        true_values <- c(0.316, 0.316, 0.139, 0.101, 0.101)
      }
    }
    
    model_name <- ifelse(condition$Model_Type == "linear", "Linear", "Full")
    distribution <- as.character(condition$Distribution)
    sample_size <- condition$N
    reliability <- condition$Rel
    
    # TRACK CONVERGENCE FOR EACH METHOD SEPARATELY
    max_reps <- 0
    # number of replications
    # could be a hard coded 1000 but as for testing before we have this 
    for (method in methods) {
      tbls <- res_list[[paste0(method, "_tables")]]
      if (!is.null(tbls)) {
        max_reps <- max(max_reps, length(tbls))
      }
    }
    
    if (max_reps == 0) next
    
    # convergence for each method separately
    method_convergence <- list()
    for (method in methods) {
      method_convergence[[method]] <- rep(FALSE, max_reps)  
    }
    
    # each method's convergence independently
    for (method in methods) {
      tbls <- res_list[[paste0(method, "_tables")]]
      if (is.null(tbls)) {
        next
      }
      
      for (rep in seq_len(max_reps)) {
        if (rep > length(tbls) || is.null(tbls[[rep]])) {
          next
        }
        
        # parameters for this replication
        extr <- extract_eta3_parameters(tbls[[rep]], method)
        if (is.null(extr)) {
          next  
        }
        
        # check if all required parameters are present and complete
        converged <- TRUE
        for (p in parameters_of_interest) {
          if (p %in% names(extr$Estimates)) {
            if (is.na(extr$Estimates[p]) || is.nan(extr$Estimates[p]) || is.infinite(extr$Estimates[p]) ||
                is.na(extr$`Standard Errors`[p]) || is.nan(extr$`Standard Errors`[p]) || is.infinite(extr$`Standard Errors`[p]) ||
                is.na(extr$`P-values`[p]) || is.nan(extr$`P-values`[p]) || is.infinite(extr$`P-values`[p]) ||
                is.na(extr$CI_lower[p]) || is.nan(extr$CI_lower[p]) || is.infinite(extr$CI_lower[p]) ||
                is.na(extr$CI_upper[p]) || is.nan(extr$CI_upper[p]) || is.infinite(extr$CI_upper[p])) {
              converged <- FALSE
              break
            }
          } else {
            converged <- FALSE
            break
          }
        }
        method_convergence[[method]][rep] <- converged
      }
    }
    
    # TRACK WARNINGS FOR CONVERGED ITERATIONS
    method_warnings <- list()
    method_warning_counts <- list()
    
    for (method in methods) {
      method_warnings[[method]] <- list()
      method_warning_counts[[method]] <- 0
      
      # if warnings exist for this method
      method_warns <- res_list[["warnings"]][[method]]
      if (!is.null(method_warns)) {
        # through each converged replication
        for (rep in which(method_convergence[[method]])) {
          if (rep <= length(method_warns) && !is.null(method_warns[[rep]]) && length(method_warns[[rep]]) > 0) {
            method_warnings[[method]][[as.character(rep)]] <- method_warns[[rep]]
            method_warning_counts[[method]] <- method_warning_counts[[method]] + 1
          }
        }
      }
    }
    
    # store non-converged iterations and warnings
    non_converged_iterations <- list()
    for (method in methods) {
      non_converged_iterations[[method]] <- which(!method_convergence[[method]])
    }
    
    convergence_outliers_details[[paste0("Condition_", i)]] <- list(
      non_converged_iterations = non_converged_iterations,
      warning_iterations = method_warnings,  # which converged iterations have warnings
      warning_counts = method_warning_counts,  # count of converged iterations with warnings
      outlier_iterations = list()  # will be populated later
    )
    
    # individual method convergence statistics
    method_convergence_counts <- list()
    method_convergence_rates <- list()
    
    for (method in methods) {
      method_convergence_counts[[method]] <- sum(method_convergence[[method]])
      method_convergence_rates[[method]] <- sum(method_convergence[[method]]) / max_reps * 100
    }
    
    # IDENTIFY VALID REPLICATIONS ACROSS ALL METHODS (CONVERGENCE + OPTIONAL WARNINGS)
    valid_reps <- rep(TRUE, max_reps)
    for (method in methods) {
      valid_reps <- valid_reps & method_convergence[[method]]
    }
    
    # optionally exclude iterations with warnings
    n_excluded_by_warnings <- 0
    if (exclude_warnings) {
      for (method in methods) {
        for (rep in seq_len(max_reps)) {
          if (as.character(rep) %in% names(method_warnings[[method]])) {
            valid_reps[rep] <- FALSE
          }
        }
      }
      converged_indices <- which(Reduce("&", method_convergence))
      n_excluded_by_warnings <- length(converged_indices) - sum(valid_reps)
    }
    
    valid_rep_indices <- which(valid_reps)
    n_total <- max_reps
    n_after_global_convergence <- length(valid_rep_indices)
    global_convergence_rate <- n_after_global_convergence / n_total * 100
    
    if (n_after_global_convergence < min_reps) {
      message(paste("Condition", i, ": Only", n_after_global_convergence, 
                    "valid replications across all methods",
                    ifelse(exclude_warnings, "(after excluding warnings).", "."),
                    "Skipping."))
      next
    }
    
    # TRACK OUTLIERS FOR EACH METHOD INDEPENDENTLY (AFTER CONVERGENCE)
    method_outliers <- list()
    method_outlier_counts <- list()
    
    for (method in methods) {
      method_outliers[[method]] <- rep(FALSE, max_reps)
      method_outlier_counts[[method]] <- 0
    }
    
    if (remove_outliers) {
      for (method in methods) {
        tbls <- res_list[[paste0(method, "_tables")]]
        if (is.null(tbls)) next
        
        # collect all estimates from valid (converged) replications
        all_estimates <- list()
        for (p in parameters_of_interest) {
          all_estimates[[p]] <- numeric()
        }
        
        for (rep in valid_rep_indices) {
          extr <- extract_eta3_parameters(tbls[[rep]], method)
          if (!is.null(extr)) {
            for (p in names(extr$Estimates)) {
              if (p %in% parameters_of_interest) {
                all_estimates[[p]] <- c(all_estimates[[p]], extr$Estimates[p])
              }
            }
          }
        }
        
        # calculate outlier bounds per parameter
        outlier_bounds <- list()
        for (p in parameters_of_interest) {
          if (length(all_estimates[[p]]) > 0) {
            q1 <- stats::quantile(all_estimates[[p]], 0.25, na.rm = TRUE)
            q3 <- stats::quantile(all_estimates[[p]], 0.75, na.rm = TRUE)
            iqr <- q3 - q1
            outlier_bounds[[p]] <- list(
              lower = q1 - outlier_threshold * iqr,
              upper = q3 + outlier_threshold * iqr
            )
          }
        }
        
        # mark outliers for each replication
        for (rep in valid_rep_indices) {
          extr <- extract_eta3_parameters(tbls[[rep]], method)
          if (!is.null(extr)) {
            is_outlier <- FALSE
            for (p in parameters_of_interest) {
              if (p %in% names(extr$Estimates) && p %in% names(outlier_bounds)) {
                est_val <- extr$Estimates[p]
                if (!is.na(est_val) && !is.nan(est_val) && !is.infinite(est_val)) {
                  if (est_val < outlier_bounds[[p]]$lower || est_val > outlier_bounds[[p]]$upper) {
                    is_outlier <- TRUE
                    break
                  }
                }
              }
            }
            method_outliers[[method]][rep] <- is_outlier
          }
        }
        
        method_outlier_counts[[method]] <- sum(method_outliers[[method]][valid_rep_indices])
      }
    }
    
    # store outlier iterations
    outlier_iterations <- list()
    for (method in methods) {
      outlier_iterations[[method]] <- which(method_outliers[[method]])
    }
    convergence_outliers_details[[paste0("Condition_", i)]]$outlier_iterations <- outlier_iterations
    
    # GLOBAL OUTLIER SYNCHRONIZATION
    valid_after_outliers <- valid_reps
    if (remove_outliers) {
      for (method in methods) {
        valid_after_outliers <- valid_after_outliers & !method_outliers[[method]]
      }
    }
    
    valid_final_indices <- which(valid_after_outliers)
    n_after_global_outlier <- length(valid_final_indices)
    global_outlier_rate <- n_after_global_outlier / n_total * 100
    
    if (n_after_global_outlier < min_reps) {
      message(paste("Condition", i, ": Only", n_after_global_outlier, 
                    "valid replications after outlier removal. Skipping."))
      next
    }
    
    # EXTRACT FILTERED DATA FOR EACH METHOD
    for (method in methods) {
      tbls <- res_list[[paste0(method, "_tables")]]
      if (is.null(tbls)) next
      
      # store convergence and outlier stats
      n_converged_this_method <- method_convergence_counts[[method]]
      convergence_rate_this_method <- method_convergence_rates[[method]]
      n_outliers_this_method <- method_outlier_counts[[method]]
      outlier_rate_this_method <- round(n_outliers_this_method / n_after_global_convergence * 100, 2)
      
      # extract data for each parameter
      params <- parameters_of_interest
      
      for (p_name in params) {
        # collect data from globally valid replications
        ests <- numeric()
        ses <- numeric()
        pvs <- numeric()
        los <- numeric()
        his <- numeric()
        
        for (rep in valid_final_indices) {
          tab <- tbls[[rep]]
          if (is.null(tab)) next
          extr <- extract_eta3_parameters(tab, method)
          if (is.null(extr)) next
          
          if (p_name %in% names(extr$Estimates)) {
            ests <- c(ests, unname(extr$Estimates[p_name]))
            ses <- c(ses, unname(extr$`Standard Errors`[p_name]))
            pvs <- c(pvs, unname(extr$`P-values`[p_name]))
            los <- c(los, unname(extr$CI_lower[p_name]))
            his <- c(his, unname(extr$CI_upper[p_name]))
          }
        }
        
        if (length(ests) == 0) next
        
        tv <- true_values[param_index[p_name]]
        
        # store summary information
        convergence_outliers_summary <- dplyr::bind_rows(
          convergence_outliers_summary,
          tibble::tibble(
            Condition = i,
            Method = toupper(method),
            Model = model_name,
            Distribution = distribution,
            SampleSize = sample_size,
            Reliability = reliability,
            Parameter = p_name,
            TrueValue = tv,
            
            # sample size tracking
            N_Total = n_total,
            N_Converged_This_Method = n_converged_this_method,
            N_Warnings_This_Method = method_warning_counts[[method]],  # warnings in converged iterations
            Convergence_Rate_This_Method = convergence_rate_this_method,
            N_After_Global_Convergence = n_after_global_convergence,
            Global_Convergence_Rate = global_convergence_rate,
            N_Excluded_By_Warnings = n_excluded_by_warnings,
            N_Outliers_This_Method = n_outliers_this_method,
            Outlier_Rate_This_Method = outlier_rate_this_method,
            N_After_Global_Outlier = n_after_global_outlier,
            Global_Outlier_Rate = global_outlier_rate,
            N_Final = length(ests),
            
            # exclusion summary
            N_Excluded_Convergence = n_total - length(which(Reduce("&", method_convergence))),
            N_Excluded_Warnings = n_excluded_by_warnings,
            N_Excluded_Outliers = n_after_global_convergence - n_after_global_outlier,
            N_Total_Excluded = n_total - length(ests),
            Percent_Total_Excluded = round((n_total - length(ests)) / n_total * 100, 2)
          )
        )
        
        # store filtered data for performance calculations
        key <- paste(i, toupper(method), p_name, sep = "_")
        filtered_data[[key]] <- list(
          condition = i,
          method = toupper(method),
          model = model_name,
          distribution = distribution,
          sample_size = sample_size,
          reliability = reliability,
          parameter = p_name,
          true_value = tv,
          estimates = ests,
          standard_errors = ses,
          p_values = pvs,
          ci_lower = los,
          ci_upper = his
        )
      }
    }
  }
  
  list(
    convergence_outliers_summary = convergence_outliers_summary,
    convergence_outliers_details = convergence_outliers_details,
    filtered_data = filtered_data
  )
}

## Function 2: Performance Metrics
CalculatePerformanceMetrics <- function(filtered_data, 
                                        convergence_outliers_summary,
                                        alpha = 0.05) {
  
  library(simhelpers)
  results_summary <- tibble::tibble()
  
  for (key in names(filtered_data)) {
    data_item <- filtered_data[[key]]
    
    # extract the stored information
    condition_num <- data_item$condition
    method_name <- data_item$method
    model_name <- data_item$model
    distribution <- data_item$distribution
    sample_size <- data_item$sample_size
    reliability <- data_item$reliability
    p_name <- data_item$parameter
    tv <- data_item$true_value
    
    # create dataframe from filtered data
    df <- tibble::tibble(
      est = data_item$estimates,
      se = data_item$standard_errors,
      p_val = data_item$p_values,
      lower_bound = data_item$ci_lower,
      upper_bound = data_item$ci_upper,
      true_param = tv
    )
    
    if (nrow(df) == 0) next
    
    # CALCULATE MEDIAN-BASED METRICS
    median_est <- stats::median(df$est)
    bias_median <- median_est - tv
    mad_est <- stats::mad(df$est, constant = 1)
    rmse_median <- sqrt(bias_median^2 + mad_est^2)
    
    if (tv != 0) {
      rel_bias_median <- bias_median / tv
      rel_bias_median_pct <- rel_bias_median * 100
    } else {
      rel_bias_median <- NA_real_
      rel_bias_median_pct <- NA_real_
    }
    
    # CALCULATE ALL MEAN-BASED METRICS USING SIMHELPERS
    # ABSOLUTE METRICS
    abs_metrics <- simhelpers::calc_absolute(
      data = df,
      estimates = est, 
      true_param = true_param,
      criteria = c("bias", "variance", "stddev", "mse", "rmse"),
      winz = Inf  # no winsorization
    )
    
    # RELATIVE METRICS (only if true value is non-zero)
    # note that pkg::fun here was not working without loading simhelpers
    if (tv != 0) {
      rel_metrics <- simhelpers::calc_relative(
        data = df,
        estimates = est,
        true_param = true_param,
        criteria = c("relative bias", "relative mse", "relative rmse"),
        winz = Inf
      )
      
      rel_bias_mean <- rel_metrics$rel_bias - 1
      rel_bias_mean_pct <- rel_bias_mean * 100
      rel_bias_mcse <- rel_metrics$rel_bias_mcse
      
      rel_mse <- rel_metrics$rel_mse
      rel_rmse <- rel_metrics$rel_rmse
      rel_mse_mcse <- rel_metrics$rel_mse_mcse
      rel_rmse_mcse <- rel_metrics$rel_rmse_mcse
    } else {
      rel_bias_mean <- NA_real_
      rel_bias_mean_pct <- NA_real_
      rel_mse <- NA_real_
      rel_rmse <- NA_real_
      rel_bias_mcse <- NA_real_
      rel_mse_mcse <- NA_real_
      rel_rmse_mcse <- NA_real_
    }
    
    # COVERAGE AND CI WIDTH
    cov_metrics <- simhelpers::calc_coverage(
      data = df,
      lower_bound = lower_bound, 
      upper_bound = upper_bound, 
      true_param = true_param,
      criteria = c("coverage", "width"),
      winz = Inf
    )
    
    # REJECTION RATES
    rej_metrics <- simhelpers::calc_rejection(
      data = df,
      p_values = p_val,
      alpha = alpha
    )
    
    rejection_rate <- rej_metrics$rej_rate * 100
    rej_mcse <- rej_metrics$rej_rate_mcse
    
    # Type I Error vs Power
    if (tv == 0) {
      type_i <- rejection_rate
      power <- NA_real_
      typei_mcse <- rej_mcse
      power_mcse <- NA_real_
    } else {
      type_i <- NA_real_
      power <- rejection_rate
      typei_mcse <- NA_real_
      power_mcse <- rej_mcse
    }
    
    # ADDITIONAL METRICS
    mean_est <- mean(df$est)
    mean_se <- mean(df$se)
    se_sd_ratio <- mean_se / abs_metrics$stddev
    
    # convergence/outlier info from summary
    conv_outlier_row <- convergence_outliers_summary[
      convergence_outliers_summary$Condition == condition_num &
        convergence_outliers_summary$Method == method_name &
        convergence_outliers_summary$Parameter == p_name, ]
    
    if (nrow(conv_outlier_row) == 0) {
      # shouldn't happen, but just in case
      warning(paste("No convergence/outlier info found for", condition_num, method_name, p_name))
      next
    }
    
    # COMPILE RESULTS
    results_summary <- dplyr::bind_rows(
      results_summary,
      tibble::tibble(
        Condition = condition_num,
        Method = method_name,
        Model = model_name,
        Distribution = distribution,
        SampleSize = sample_size,
        Reliability = reliability,
        Parameter = p_name,
        TrueValue = tv,
        
        # Mean-based metrics
        MeanEstimate = mean_est,
        Bias_Mean = abs_metrics$bias,
        RelativeBias_Mean = rel_bias_mean,
        PercentRelativeBias_Mean = rel_bias_mean_pct,
        
        # Median-based metrics
        MedianEstimate = median_est,
        Bias_Median = bias_median,
        RelativeBias_Median = rel_bias_median,
        PercentRelativeBias_Median = rel_bias_median_pct,
        MAD = mad_est,
        
        # Variability and accuracy metrics
        Variance = abs_metrics$var,
        SD = abs_metrics$stddev,
        MSE_Mean = abs_metrics$mse,
        RMSE_Mean = abs_metrics$rmse,
        RMSE_Median = rmse_median,
        
        # Relative metrics
        Relative_MSE = rel_mse,
        Relative_RMSE = rel_rmse,
        
        # Standard error metrics
        MeanSE = mean_se,
        SE_SD_Ratio = se_sd_ratio,
        
        # Coverage and CI metrics
        CoverageRate = cov_metrics$coverage * 100,
        CI_Width = cov_metrics$width,
        
        # Testing metrics
        RejectionRate = rejection_rate,
        TypeI_Error = type_i,
        Power = power,
        
        # Monte Carlo standard errors
        Bias_Mean_MCSE = abs_metrics$bias_mcse,
        RelativeBias_MCSE = rel_bias_mcse,
        Variance_MCSE = abs_metrics$var_mcse,
        SD_MCSE = abs_metrics$stddev_mcse,
        MSE_Mean_MCSE = abs_metrics$mse_mcse,
        RMSE_Mean_MCSE = abs_metrics$rmse_mcse,
        Relative_MSE_MCSE = rel_mse_mcse,
        Relative_RMSE_MCSE = rel_rmse_mcse,
        CoverageRate_MCSE = cov_metrics$coverage_mcse,
        CI_Width_MCSE = cov_metrics$width_mcse,
        RejectionRate_MCSE = rej_mcse,
        TypeI_Error_MCSE = typei_mcse,
        Power_MCSE = power_mcse,
        
        # Add convergence/outlier tracking from summary
        N_Total = conv_outlier_row$N_Total,
        N_Converged_This_Method = conv_outlier_row$N_Converged_This_Method,
        N_Warnings_This_Method = conv_outlier_row$N_Warnings_This_Method,
        Convergence_Rate_This_Method = conv_outlier_row$Convergence_Rate_This_Method,
        N_After_Global_Convergence = conv_outlier_row$N_After_Global_Convergence,
        Global_Convergence_Rate = conv_outlier_row$Global_Convergence_Rate,
        N_Excluded_By_Warnings = conv_outlier_row$N_Excluded_By_Warnings,
        N_Outliers_This_Method = conv_outlier_row$N_Outliers_This_Method,
        Outlier_Rate_This_Method = conv_outlier_row$Outlier_Rate_This_Method,
        N_After_Global_Outlier = conv_outlier_row$N_After_Global_Outlier,
        Global_Outlier_Rate = conv_outlier_row$Global_Outlier_Rate,
        N_Final = conv_outlier_row$N_Final,
        N_Excluded_Convergence = conv_outlier_row$N_Excluded_Convergence,
        N_Excluded_Warnings = conv_outlier_row$N_Excluded_Warnings,
        N_Excluded_Outliers = conv_outlier_row$N_Excluded_Outliers,
        N_Total_Excluded = conv_outlier_row$N_Total_Excluded,
        Percent_Total_Excluded = conv_outlier_row$Percent_Total_Excluded
      )
    )
  }
  
  results_summary
}

load("Simulations/Study_1/Data/Results_Study_1_final.RData")

# wrapper for manuscript currently 
CalculatePerformance <- function(all_results,
                                 parameters_of_interest = c("eta1","eta2","eta1:eta2","eta1:eta1","eta2:eta2"),
                                 remove_outliers = TRUE,
                                 outlier_threshold = 3, 
                                 alpha = 0.05,
                                 min_reps = 10,
                                 exclude_warnings = FALSE,
                                 return_convergence_details = TRUE) {
  
  # convergence and outliers
  extraction_results <- ExtractConvergenceOutliers(
    all_results = all_results,
    parameters_of_interest = parameters_of_interest,
    remove_outliers = remove_outliers,
    outlier_threshold = outlier_threshold,
    min_reps = min_reps,
    exclude_warnings = exclude_warnings
  )
  
  # performance metrics
  performance_results <- CalculatePerformanceMetrics(
    filtered_data = extraction_results$filtered_data,
    convergence_outliers_summary = extraction_results$convergence_outliers_summary,
    alpha = alpha
  )
  
  if (return_convergence_details) {
    list(
      results = performance_results,
      convergence_outliers_details = extraction_results$convergence_outliers_details
    )
  } else {
    performance_results
  }
}


results_study_1 <- CalculatePerformance(
   all_results,
   parameters_of_interest = c("eta1","eta2","eta1:eta2","eta1:eta1","eta2:eta2"),
   remove_outliers = TRUE,
   outlier_threshold = 3,
   alpha = 0.05,
   min_reps = 10,
   exclude_warnings = FALSE,  # TRUE to exclude iterations with warnings
   return_convergence_details = FALSE
)
