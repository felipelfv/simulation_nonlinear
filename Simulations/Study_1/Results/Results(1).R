############################ 1. General Information ############################

# See README file for more information concerning this file. 
# This file contains the code necessary to calculate the performance results 
# reported in the manuscript. 
# The final dataset contains all the necessary metrics for plotting and
# creating tables reported in the manuscript for Simulation 1.

############################### 2. Documentation ###############################

#' extract_params: Extract Structural Parameters from Method Output
#' 
#' @param table         Data.frame. Parameter table from estimation method output
#' @param method        Character. Type of method: "lsam", "lms", "qml", or "upi"
#' 
#' @return List containing:
#'   - Estimates:       Named vector of parameter estimates
#'   - Standard Errors: Named vector of standard errors
#'   - P-values:        Named vector of p-values
#'   - CI_lower:        Named vector of confidence interval lower bounds
#'   - CI_upper:        Named vector of confidence interval upper bounds
#' 
#' @note Method-specific handling:
#'   - QML/LMS: Use "std.error" and "p.value" columns; require parameter reordering
#'   - LSAM/UPI: Use "se" and "pvalue" columns; maintain original order
#'   - Automatic detection of model complexity based on row count (3 = linear, 5 = full)
#'   
#' @details Parameter extraction:
#'   - Linear model (3 params): eta1, eta2, eta1:eta2
#'   - Full model (5 params): eta1, eta2, eta1:eta2, eta1:eta1, eta2:eta2
#'   - LMS/QML output eta1:eta1 before eta1:eta2, requiring position swap
#'
#' ExtractConvergenceOutliers: Process Convergence and Identify Outliers
#' 
#' @param all_results              List. Complete simulation results from all conditions
#' @param parameters_of_interest   Character vector. Parameters to extract 
#'                                 (default = c("eta1","eta2","eta1:eta2","eta1:eta1","eta2:eta2"))
#' @param remove_outliers          Logical. Whether to identify and exclude outliers (default = TRUE)
#' @param outlier_threshold        Numeric. IQR multiplier for outlier detection (default = 3)
#' @param min_reps                 Integer. Minimum replications required per condition (default = 10)
#' @param exclude_warnings         Logical. Whether to exclude iterations with warnings (default = FALSE)
#' 
#' @return List containing:
#'   - convergence_outliers_summary: Tibble with convergence/outlier statistics per condition/method/parameter
#'   - convergence_outliers_details: List with detailed iteration-level information
#'   - filtered_data:               List with clean data after outlier removal (for mean-based metrics)
#'   - converged_data:              List with converged data before outlier removal (for robust metrics)
#' 
#' @details Processing Pipeline:
#'   1. Check convergence for each method independently (non-missing, finite values)
#'   2. Track warnings in converged iterations only
#'   3. Identify valid replications (converged across all methods - global convergence)
#'   4. Optionally exclude iterations with warnings (if exclude_warnings = TRUE)
#'   5. Detect outliers using IQR method (Q1 - threshold*IQR, Q3 + threshold*IQR)
#'   6. Apply global synchronization (exclude if outlier in any method)
#'   7. Extract filtered data for performance metrics
#'   8. Skip conditions with fewer than min_reps valid replications
#' 
#' @details Convergence Criteria:
#'   Per-method convergence requires:
#'   - All required parameters present in output
#'   - No NA, NaN, or Inf values in estimates, SEs, p-values, or CIs
#'   - Method-specific convergence tracked independently
#'   
#'   Global convergence requires:
#'   - Convergence across ALL methods simultaneously for a given replication
#'   - Used as baseline for outlier detection
#' 
#' @details Outlier Detection Method (Boomsma, 2013):
#'   - Applied separately to each parameter
#'   - Calculate Q1 (25th percentile) and Q3 (75th percentile)
#'   - Compute IQR = Q3 - Q1
#'   - Lower bound = Q1 - threshold * IQR
#'   - Upper bound = Q3 + threshold * IQR
#'   - Flag replication as outlier if ANY parameter exceeds bounds
#'   - Global synchronization: exclude from all methods if outlier in any
#' 
#' @details Summary Statistics Tracked (per condition/method/parameter):
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
#'   - N_Excluded_Convergence:      Excluded due to non-convergence
#'   - N_Excluded_Warnings:         Excluded due to warnings
#'   - N_Excluded_Outliers:         Excluded due to outlier detection
#'   - N_Total_Excluded:            Total excluded replications
#'   - Percent_Total_Excluded:      Percentage of total excluded
#'
#' CalculatePerformanceMetrics: Calculate All Performance Metrics
#' 
#' @param filtered_data                 List. Filtered data from ExtractConvergenceOutliers (after outlier removal)
#' @param converged_data                List. Converged data from ExtractConvergenceOutliers (before outlier removal)
#' @param convergence_outliers_summary  Tibble. Summary statistics from ExtractConvergenceOutliers
#' @param alpha                         Numeric. Significance level for tests (default = 0.05)
#' 
#' @return Tibble containing all performance metrics per condition/method/parameter:
#' 
#' @details Mean-based metrics (calculated from filtered_data after outlier removal):
#'   - MeanEstimate:                Average of parameter estimates
#'   - Bias_Mean:                   Mean estimate - true value
#'   - RelativeBias_Mean:           (Bias / true value) - 1 (NA if true value = 0)
#'   - PercentRelativeBias_Mean:    Relative bias * 100
#'   - Variance:                    Variance of estimates
#'   - SD:                          Standard deviation of estimates
#'   - MSE_Mean:                    Mean squared error
#'   - RMSE_Mean:                   Root mean squared error
#'   
#' @details Median-based metrics (calculated from converged_data before outlier removal):
#'   - MedianEstimate:              Median of parameter estimates
#'   - Bias_Median:                 Median estimate - true value  
#'   - RelativeBias_Median:         Median bias / true value (NA if true value = 0)
#'   - PercentRelativeBias_Median:  Relative median bias * 100
#'   - MAD:                         Median absolute deviation (standardized with factor 1.4826)
#'   - RMSE_Median:                 sqrt(Bias_Median^2 + MAD^2)
#'   
#' @details Trimmed mean-based metrics (calculated from converged_data before outlier removal):
#'   - TrimmedEstimate:             20% trimmed mean of parameter estimates
#'   - Bias_Trimmed:                Trimmed mean estimate - true value
#'   - RMSE_Trimmed:                sqrt(Bias_Trimmed^2 + MAD^2)
#'   
#' @details Standard error metrics:
#'   - MeanSE:                      Average of standard errors (from filtered_data)
#'   - TrimmedSE:                   20% trimmed mean of standard errors (from converged_data)
#'   - SE_SD_Ratio:                 Mean SE / SD of estimates (from filtered_data)
#'   - RSB:                         Robust relative SE bias = TrimmedSE/MAD - 1 
#'   
#' @details Coverage and confidence interval metrics:
#'   - CoverageRate:                Percentage of CIs containing true value
#'   - CI_Width:                    Average width of confidence intervals
#'   
#' @details Hypothesis testing metrics:
#'   - RejectionRate:               Percentage of p-values < alpha
#'   - TypeI_Error:                 Rejection rate when true value = 0
#'   - Power:                       Rejection rate when true value not 0
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
#'   
#' @note Performance calculations use simhelpers package for Monte Carlo metrics
#' @note Robust metrics use converged data to avoid influence of outliers
#' @note Mean-based metrics use filtered data after outlier removal
#'
#' CalculatePerformance: Wrapper Function for Complete Simulation 1 Analysis
#' 
#' @param all_results                  List. Complete simulation results from all conditions
#' @param parameters_of_interest       Character vector. Parameters to extract 
#'                                     (default = c("eta1","eta2","eta1:eta2","eta1:eta1","eta2:eta2"))
#' @param remove_outliers              Logical. Whether to remove outliers (default = TRUE)
#' @param outlier_threshold            Numeric. IQR multiplier for outlier detection (default = 3)
#' @param alpha                        Numeric. Significance level for hypothesis tests (default = 0.05)
#' @param min_reps                     Integer. Minimum replications required per condition (default = 10)
#' @param exclude_warnings             Logical. Whether to exclude iterations with warnings (default = FALSE)
#' @param return_convergence_details   Logical. Whether to return detailed convergence info (default = FALSE)
#' 
#' @return If return_convergence_details = TRUE:
#'   - List containing:
#'     * results: Complete performance metrics tibble
#'     * convergence_outliers_details: Detailed iteration-level information
#'   
#'   If return_convergence_details = FALSE:
#'   - Tibble with performance metrics only
#'   
#' @description Convenience wrapper that combines extraction and calculation steps into a single call.
#'              Handles the complete analysis pipeline from raw simulation results to final metrics.
#'              
#' @note Simulation 1-Specific Features:
#'   - Single structural equation (eta3)
#'   - Four methods: LSAM, LMS, QML, UPI
#'   - Linear model (3 parameters) or Full model (5 parameters)
#'   - True values depend on model type:
#'     * Linear: c(0.316, 0.316, 0, 0, 0)
#'     * Full: c(0.316, 0.316, 0.139, 0.101, 0.101)
#'     
#' @note Dependencies:
#'   Required packages: dplyr, tibble, simhelpers
#'   
#' @examples
#' \dontrun{
#' results_study_1 <- CalculatePerformance(
#'   all_results,
#'   parameters_of_interest = c("eta1","eta2","eta1:eta2","eta1:eta1","eta2:eta2"),
#'   remove_outliers = TRUE,
#'   outlier_threshold = 3,
#'   alpha = 0.05,
#'   min_reps = 10,
#'   exclude_warnings = FALSE,
#'   return_convergence_details = FALSE
#' )
#' }

############################### 3. Functions ###################################

extract_params <- function(table, method) {
  if(is.null(table)) return(NULL)
  
  rows <- table[table$lhs == "eta3" & table$op == "~", ]
  params <- if(nrow(rows) == 3) {
    c("eta1", "eta2", "eta1:eta2")
  } else {
    c("eta1", "eta2", "eta1:eta2", "eta1:eta1", "eta2:eta2")
  }
  
  # qml and lms have different names (in contrast to sam/upi)
  se_col <- if(method %in% c("lms", "qml")) "std.error" else "se"
  pval_col <- if(method %in% c("lms", "qml")) "p.value" else "pvalue"
  
  # LMS/QML need reordering (they output eta1:eta1 before eta1:eta2)
  if(method %in% c("lms", "qml") && nrow(rows) == 5) {
    idx <- c(1, 2, 4, 3, 5)  # swap positions 3 and 4
    rows <- rows[idx, ]
  }
  
  list(
    Estimates = setNames(rows$est, params),
    "Standard Errors" = setNames(rows[[se_col]], params),
    "P-values" = setNames(rows[[pval_col]], params),
    CI_lower = setNames(rows$ci.lower, params),
    CI_upper = setNames(rows$ci.upper, params)
  )
}

# main extraction function
ExtractConvergenceOutliers <- function(all_results,
                                       parameters_of_interest = c("eta1","eta2","eta1:eta2","eta1:eta1","eta2:eta2"),
                                       remove_outliers = TRUE,
                                       outlier_threshold = 3,
                                       min_reps = 10,
                                       exclude_warnings = FALSE) {
  
  options(scipen = 999)
  methods <- c("lsam", "lms", "qml", "upi")
  # mapping from parameter names to true value positions:
  param_index <- setNames(1:5, parameters_of_interest)
  
  # null-coalescing operator:
  `%||%` <- function(x, y) if(is.null(x)) y else x
  
  convergence_outliers_summary <- tibble::tibble()
  convergence_outliers_details <- list()
  filtered_data <- list()  # after outlier removal (for mean-based metrics)
  converged_data <- list()  # after convergence only (for robust metrics)
  
  for (i in seq_along(all_results)) {
    condition <- all_results[[i]]$condition
    res_list <- all_results[[i]]$results
    
    # get true values
    true_values <- all_results[[i]]$true_parameters %||% 
      if(condition$Model_Type == "linear") c(0.316, 0.316, 0, 0, 0) 
    else c(0.316, 0.316, 0.139, 0.101, 0.101)
    
    # get max replications
    # in principle should always be 1000 but just a sanity check
    max_reps <- max(sapply(methods, function(m) 
      length(res_list[[paste0(m, "_tables")]]) %||% 0))
    
    if (max_reps == 0) next
    
    # check for convergence
    # here we create a convergence indicator for each method
    method_convergence <- setNames(lapply(methods, function(method) {
      tbls <- res_list[[paste0(method, "_tables")]]
      if (is.null(tbls)) return(rep(FALSE, max_reps))
      # then check each replication
      vapply(seq_len(max_reps), function(rep) {
        if (rep > length(tbls) || is.null(tbls[[rep]])) return(FALSE)
        extr <- extract_params(tbls[[rep]], method)
        if (is.null(extr)) return(FALSE)
        # note that we check if all parameters exist and are valid 
        # (e.g.,  no NA, NaN, or infinite values in any statistic)
        all(sapply(parameters_of_interest, function(p) {
          p %in% names(extr$Estimates) && 
            !any(is.na(c(extr$Estimates[p], extr$`Standard Errors`[p], 
                         extr$`P-values`[p], extr$CI_lower[p], extr$CI_upper[p])) |
                   is.nan(c(extr$Estimates[p], extr$`Standard Errors`[p], 
                            extr$`P-values`[p], extr$CI_lower[p], extr$CI_upper[p])) |
                   is.infinite(c(extr$Estimates[p], extr$`Standard Errors`[p], 
                                 extr$`P-values`[p], extr$CI_lower[p], extr$CI_upper[p])))
        }))
      }, logical(1))
    }), methods)
    
    # warnings processing
    # we only keep warnings from converged replications
    method_warnings <- setNames(lapply(methods, function(method) {
      warns <- res_list$warnings[[method]]
      if (is.null(warns)) return(list())
      
      conv_reps <- which(method_convergence[[method]])
      warn_list <- list()
      for(rep in conv_reps) {
        if(rep <= length(warns) && !is.null(warns[[rep]]) && length(warns[[rep]]) > 0) {
          warn_list[[as.character(rep)]] <- warns[[rep]]
        }
      }
      warn_list
    }), methods)
    
    method_warning_counts <- sapply(method_warnings, length)
    
    convergence_outliers_details[[paste0("Condition_", i)]] <- list(
      non_converged_iterations = lapply(method_convergence, function(x) which(!x)),
      warning_iterations = method_warnings,
      warning_counts = method_warning_counts,
      outlier_iterations = list()
    )
    
    # convergence stats
    method_convergence_counts <- sapply(method_convergence, sum)
    method_convergence_rates <- sapply(method_convergence, function(x) mean(x) * 100)
    # global valid reps
    # meaning a logical vector where TRUE means all methods converged
    valid_reps <- Reduce("&", method_convergence)
    
    # if exclude warnings
    n_excluded_by_warnings <- 0
    if (exclude_warnings) {
      warn_reps <- unique(unlist(lapply(method_warnings, function(w) as.numeric(names(w)))))
      if(length(warn_reps) > 0) {
        valid_reps[warn_reps] <- FALSE
        n_excluded_by_warnings <- sum(Reduce("&", method_convergence)) - sum(valid_reps)
      }
    }
    
    valid_rep_indices <- which(valid_reps)
    n_after_global_convergence <- length(valid_rep_indices)
    
    if (n_after_global_convergence < min_reps) {
      message(sprintf("Condition %d: Only %d valid replications. Skipping.", 
                      i, n_after_global_convergence))
      next
    }
    
    # store converged data (BEFORE outlier removal) for robust metrics
    model_name <- ifelse(condition$Model_Type == "linear", "Linear", "Full")
    distribution <- as.character(condition$Distribution)
    sample_size <- condition$N
    reliability <- condition$Rel
    
    for(method in methods) {
      tbls <- res_list[[paste0(method, "_tables")]]
      if (is.null(tbls)) next
      
      for(p_name in parameters_of_interest) {
        # extract converged data (no outlier removal yet)
        converged_list <- lapply(valid_rep_indices, function(rep) {
          extr <- extract_params(tbls[[rep]], method)
          if(!is.null(extr) && p_name %in% names(extr$Estimates)) {
            list(est = unname(extr$Estimates[p_name]),
                 se = unname(extr$`Standard Errors`[p_name]),
                 pv = unname(extr$`P-values`[p_name]),
                 lo = unname(extr$CI_lower[p_name]),
                 hi = unname(extr$CI_upper[p_name]))
          }
        })
        converged_list <- Filter(Negate(is.null), converged_list)
        
        if(length(converged_list) > 0) {
          tv <- true_values[param_index[p_name]]
          key <- paste(i, toupper(method), p_name, sep = "_")
          converged_data[[key]] <- list(
            condition = i,
            method = toupper(method),
            model = model_name,
            distribution = distribution,
            sample_size = sample_size,
            reliability = reliability,
            parameter = p_name,
            true_value = tv,
            estimates = sapply(converged_list, "[[", "est"),
            standard_errors = sapply(converged_list, "[[", "se"),
            p_values = sapply(converged_list, "[[", "pv"),
            ci_lower = sapply(converged_list, "[[", "lo"),
            ci_upper = sapply(converged_list, "[[", "hi")
          )
        }
      }
    }
    
    # outlier detection
    method_outliers <- setNames(lapply(methods, function(method) {
      outliers <- rep(FALSE, max_reps)
      if (!remove_outliers) return(outliers)
      
      tbls <- res_list[[paste0(method, "_tables")]]
      if (is.null(tbls)) return(outliers)
      
      # collect estimates
      param_ests <- setNames(lapply(parameters_of_interest, function(p) {
        unlist(lapply(valid_rep_indices, function(rep) {
          extr <- extract_params(tbls[[rep]], method)
          if(!is.null(extr) && p %in% names(extr$Estimates)) extr$Estimates[p]
        }))
      }), parameters_of_interest)
      
      # bounds used in lonati et al., 2024 
      # based on boomsma, (2013). "Reporting Monte Carlo studies in SEM (...)"
      bounds <- lapply(param_ests, function(ests) {
        if(length(ests) == 0) return(c(-Inf, Inf))
        q <- quantile(ests, c(0.25, 0.75), na.rm = TRUE)
        iqr <- diff(q)
        c(q[1] - outlier_threshold * iqr, q[2] + outlier_threshold * iqr)
      })
      
      # mark as outlier if ANY parameter is outside bounds
      for(rep in valid_rep_indices) {
        extr <- extract_params(tbls[[rep]], method)
        if(!is.null(extr)) {
          outliers[rep] <- any(sapply(parameters_of_interest, function(p) {
            p %in% names(extr$Estimates) && 
              !is.na(extr$Estimates[p]) && !is.nan(extr$Estimates[p]) && !is.infinite(extr$Estimates[p]) &&
              (extr$Estimates[p] < bounds[[p]][1] || extr$Estimates[p] > bounds[[p]][2])
          }))
        }
      }
      outliers
    }), methods)
    
    method_outlier_counts <- sapply(method_outliers, function(x) sum(x[valid_rep_indices]))
    convergence_outliers_details[[paste0("Condition_", i)]]$outlier_iterations <- 
      lapply(method_outliers, which)
    
    # final valid indices
    # remove replications that are outliers in ANY method (global outlier removal)
    valid_after_outliers <- valid_reps & !Reduce("|", method_outliers)
    valid_final_indices <- which(valid_after_outliers)
    n_after_global_outlier <- length(valid_final_indices)
    
    if (n_after_global_outlier < min_reps) {
      message(sprintf("Condition %d: Only %d replications after outliers. Skipping.", 
                      i, n_after_global_outlier))
      next
    }
    
    # extract data and build summary
    model_name <- ifelse(condition$Model_Type == "linear", "Linear", "Full")
    distribution <- as.character(condition$Distribution)
    sample_size <- condition$N
    reliability <- condition$Rel
    
    for(method in methods) {
      tbls <- res_list[[paste0(method, "_tables")]]
      if (is.null(tbls)) next
      
      for(p_name in parameters_of_interest) {
        # extract data
        data_list <- lapply(valid_final_indices, function(rep) {
          extr <- extract_params(tbls[[rep]], method)
          if(!is.null(extr) && p_name %in% names(extr$Estimates)) {
            list(est = unname(extr$Estimates[p_name]),
                 se = unname(extr$`Standard Errors`[p_name]),
                 pv = unname(extr$`P-values`[p_name]),
                 lo = unname(extr$CI_lower[p_name]),
                 hi = unname(extr$CI_upper[p_name]))
          }
        })
        data_list <- Filter(Negate(is.null), data_list)
        
        if(length(data_list) == 0) next
        
        tv <- true_values[param_index[p_name]]
        
        # summary row
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
            N_Total = max_reps,
            N_Converged_This_Method = method_convergence_counts[method],
            N_Warnings_This_Method = method_warning_counts[method],
            Convergence_Rate_This_Method = method_convergence_rates[method],
            N_After_Global_Convergence = n_after_global_convergence,
            Global_Convergence_Rate = n_after_global_convergence / max_reps * 100,
            N_Excluded_By_Warnings = n_excluded_by_warnings,
            N_Outliers_This_Method = method_outlier_counts[method],
            Outlier_Rate_This_Method = round(method_outlier_counts[method] / n_after_global_convergence * 100, 2),
            N_After_Global_Outlier = n_after_global_outlier,
            Global_Outlier_Rate = n_after_global_outlier / max_reps * 100,
            N_Final = length(data_list),
            N_Excluded_Convergence = max_reps - sum(Reduce("&", method_convergence)),
            N_Excluded_Warnings = n_excluded_by_warnings,
            N_Excluded_Outliers = n_after_global_convergence - n_after_global_outlier,
            N_Total_Excluded = max_reps - length(data_list),
            Percent_Total_Excluded = round((max_reps - length(data_list)) / max_reps * 100, 2)
          )
        )
        
        # filtered data
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
          estimates = sapply(data_list, "[[", "est"),
          standard_errors = sapply(data_list, "[[", "se"),
          p_values = sapply(data_list, "[[", "pv"),
          ci_lower = sapply(data_list, "[[", "lo"),
          ci_upper = sapply(data_list, "[[", "hi")
        )
      }
    }
  }
  
  list(
    convergence_outliers_summary = convergence_outliers_summary,
    convergence_outliers_details = convergence_outliers_details,
    filtered_data = filtered_data,
    converged_data = converged_data
  )
}

# performance metrics calculation
CalculatePerformanceMetrics <- function(filtered_data,
                                        converged_data,
                                        convergence_outliers_summary,
                                        alpha = 0.05) {
  
  # %||% if not available
  `%||%` <- function(x, y) if(is.null(x)) y else x
  
  results_summary <- tibble::tibble()
  
  for (key in names(filtered_data)) {
    d <- filtered_data[[key]]
    tv <- d$true_value
    
    df <- tibble::tibble(
      est = d$estimates,
      se = d$standard_errors,
      p_val = d$p_values,
      lower_bound = d$ci_lower,
      upper_bound = d$ci_upper,
      true_param = tv
    )
    
    if (nrow(df) == 0) next
    
    # metrics using simhelpers
    abs_metrics <- simhelpers::calc_absolute(df, est, true_param, 
                                             c("bias", "variance", "stddev", "mse", "rmse"), winz = Inf)
    cov_metrics <- simhelpers::calc_coverage(df, lower_bound, upper_bound, true_param, 
                                             c("coverage", "width"), winz = Inf)
    rej_metrics <- simhelpers::calc_rejection(df, p_val, alpha)
    
    # relative metrics if true value non-zero
    if(tv != 0) {
      rel_metrics <- simhelpers::calc_relative(df, est, true_param, 
                                               c("relative bias", "relative mse", "relative rmse"), winz = Inf)
      rel_bias_mean <- rel_metrics$rel_bias - 1
      rel_bias_mean_pct <- rel_bias_mean * 100
    } else {
      rel_metrics <- list(rel_mse = NA, rel_rmse = NA, 
                          rel_bias_mcse = NA, rel_mse_mcse = NA, rel_rmse_mcse = NA)
      rel_bias_mean <- NA
      rel_bias_mean_pct <- NA
    }
    
    # ROBUST METRICS: use converged data (BEFORE outlier removal)
    d_conv <- converged_data[[key]]
    if (is.null(d_conv)) {
      # fallback to filtered if converged not available
      # for the sake of double checking
      d_conv <- d
      warning(paste("Converged data not found for", key, "- using filtered data"))
    }
    
    df_conv <- tibble::tibble(
      est = d_conv$estimates,
      se = d_conv$standard_errors
    )
    
    # median metrics (from converged data)
    median_est <- median(df_conv$est)
    bias_median <- median_est - tv
    mad_est <- mad(df_conv$est)  # uses default constant = 1.4826 to scale to SD
    
    # trimmed mean metrics (20% trimmed mean as per Wilcox, 2016; Lai & Hsiao, 2022 from converged data)
    trimmed_est <- mean(df_conv$est, trim = 0.20)
    bias_trimmed <- trimmed_est - tv
    rmse_trimmed <- sqrt(bias_trimmed^2 + mad_est^2)
    
    # robust relative SE bias (RSB, from converged data)
    # 20% trimmed mean of SEs divided by MAD minus 1
    trimmed_se <- mean(df_conv$se, trim = 0.20)
    rsb <- trimmed_se / mad_est - 1
    
    # get convergence info
    conv_row <- convergence_outliers_summary[
      convergence_outliers_summary$Condition == d$condition &
        convergence_outliers_summary$Method == d$method &
        convergence_outliers_summary$Parameter == d$parameter, ]
    
    if (nrow(conv_row) == 0) {
      warning(paste("No convergence info for", d$condition, d$method, d$parameter))
      next
    }
    
    # results row
    results_summary <- dplyr::bind_rows(
      results_summary,
      tibble::tibble(
        Condition = d$condition,
        Method = d$method,
        Model = d$model,
        Distribution = d$distribution,
        SampleSize = d$sample_size,
        Reliability = d$reliability,
        Parameter = d$parameter,
        TrueValue = tv,
        
        # estimates
        MeanEstimate = mean(df$est),
        MedianEstimate = median_est,
        TrimmedEstimate = trimmed_est,
        
        # bias metrics
        Bias_Mean = abs_metrics$bias,
        Bias_Median = bias_median,
        Bias_Trimmed = bias_trimmed,
        RelativeBias_Mean = rel_bias_mean,
        RelativeBias_Median = if(tv != 0) bias_median / tv else NA,
        PercentRelativeBias_Mean = rel_bias_mean_pct,
        PercentRelativeBias_Median = if(tv != 0) bias_median / tv * 100 else NA,
        
        # variability
        MAD = mad_est,
        Variance = abs_metrics$var,
        SD = abs_metrics$stddev,
        MSE_Mean = abs_metrics$mse,
        RMSE_Mean = abs_metrics$rmse,
        RMSE_Median = sqrt(bias_median^2 + mad_est^2),
        RMSE_Trimmed = rmse_trimmed,
        Relative_MSE = rel_metrics$rel_mse %||% NA,
        Relative_RMSE = rel_metrics$rel_rmse %||% NA,
        
        # SE metrics
        MeanSE = mean(df$se),
        TrimmedSE = trimmed_se,
        SE_SD_Ratio = mean(df$se) / abs_metrics$stddev,
        RSB = rsb,
        
        # coverage
        CoverageRate = cov_metrics$coverage * 100,
        CI_Width = cov_metrics$width,
        
        # rejection rates
        RejectionRate = rej_metrics$rej_rate * 100,
        TypeI_Error = if(tv == 0) rej_metrics$rej_rate * 100 else NA,
        Power = if(tv != 0) rej_metrics$rej_rate * 100 else NA,
        
        # MCSEs
        Bias_Mean_MCSE = abs_metrics$bias_mcse,
        RelativeBias_MCSE = rel_metrics$rel_bias_mcse %||% NA,
        Variance_MCSE = abs_metrics$var_mcse,
        SD_MCSE = abs_metrics$stddev_mcse,
        MSE_Mean_MCSE = abs_metrics$mse_mcse,
        RMSE_Mean_MCSE = abs_metrics$rmse_mcse,
        Relative_MSE_MCSE = rel_metrics$rel_mse_mcse %||% NA,
        Relative_RMSE_MCSE = rel_metrics$rel_rmse_mcse %||% NA,
        CoverageRate_MCSE = cov_metrics$coverage_mcse,
        CI_Width_MCSE = cov_metrics$width_mcse,
        RejectionRate_MCSE = rej_metrics$rej_rate_mcse,
        TypeI_Error_MCSE = if(tv == 0) rej_metrics$rej_rate_mcse else NA,
        Power_MCSE = if(tv != 0) rej_metrics$rej_rate_mcse else NA,
        
        # convergence and outliers columns
        N_Total = conv_row$N_Total,
        N_Converged_This_Method = conv_row$N_Converged_This_Method,
        N_Warnings_This_Method = conv_row$N_Warnings_This_Method,
        Convergence_Rate_This_Method = conv_row$Convergence_Rate_This_Method,
        N_After_Global_Convergence = conv_row$N_After_Global_Convergence,
        Global_Convergence_Rate = conv_row$Global_Convergence_Rate,
        N_Excluded_By_Warnings = conv_row$N_Excluded_By_Warnings,
        N_Outliers_This_Method = conv_row$N_Outliers_This_Method,
        Outlier_Rate_This_Method = conv_row$Outlier_Rate_This_Method,
        N_After_Global_Outlier = conv_row$N_After_Global_Outlier,
        Global_Outlier_Rate = conv_row$Global_Outlier_Rate,
        N_Final = conv_row$N_Final,
        N_Excluded_Convergence = conv_row$N_Excluded_Convergence,
        N_Excluded_Warnings = conv_row$N_Excluded_Warnings,
        N_Excluded_Outliers = conv_row$N_Excluded_Outliers,
        N_Total_Excluded = conv_row$N_Total_Excluded,
        Percent_Total_Excluded = conv_row$Percent_Total_Excluded
      )
    )
  }
  
  results_summary
}

# wrapper function
# this actually not necessary (personal preference)
CalculatePerformance <- function(all_results,
                                 parameters_of_interest = c("eta1","eta2","eta1:eta2","eta1:eta1","eta2:eta2"),
                                 remove_outliers = TRUE,
                                 outlier_threshold = 3,
                                 alpha = 0.05,
                                 min_reps = 10,
                                 exclude_warnings = FALSE,
                                 return_convergence_details = FALSE) {
  
  extraction_results <- ExtractConvergenceOutliers(
    all_results = all_results,
    parameters_of_interest = parameters_of_interest,
    remove_outliers = remove_outliers,
    outlier_threshold = outlier_threshold,
    min_reps = min_reps,
    exclude_warnings = exclude_warnings
  )
  
  performance_results <- CalculatePerformanceMetrics(
    filtered_data = extraction_results$filtered_data,
    converged_data = extraction_results$converged_data,
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
