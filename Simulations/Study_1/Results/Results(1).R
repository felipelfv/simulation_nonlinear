############################ 1. General Information ############################

# This file contains flexible convergence checking with optional strictness levels:
# - Basic (default): NA/NaN/Inf checks only (original)
# - Strict: Adds positive SEs, positive variances, positive definite factor cov matrix

############################### 3. Functions ###################################

# helper functions
`%||%` <- function(x, y) if (is.null(x)) y else x

is_valid <- function(x) !any(is.na(x) | is.nan(x) | is.infinite(x))

extract_data_list <- function(tbls, rep_indices, method, p_name) {
  out <- lapply(rep_indices, function(rep) {
    extr <- extract_params_study_1(tbls[[rep]], method)
    if (!is.null(extr) && p_name %in% names(extr$Estimates)) {
      list(est = unname(extr$Estimates[p_name]),
           se  = unname(extr$`Standard Errors`[p_name]),
           pv  = unname(extr$`P-values`[p_name]),
           lo  = unname(extr$CI_lower[p_name]),
           hi  = unname(extr$CI_upper[p_name]))
    }
  })
  Filter(Negate(is.null), out)
}

# parameter estimates table has a different order depending on the method:
extract_params_study_1 <- function(table, method) {
  if (is.null(table)) return(NULL)
  
  rows <- table[table$lhs == "eta3" & table$op == "~", ]
  params <- if (nrow(rows) == 3) {
    c("eta1", "eta2", "eta1:eta2")
  } else {
    c("eta1", "eta2", "eta1:eta2", "eta1:eta1", "eta2:eta2")
  }
  
  se_col <- if (method %in% c("lms", "qml")) "std.error" else "se"
  pval_col <- if (method %in% c("lms", "qml")) "p.value" else "pvalue"
  
  if (method %in% c("lms", "qml") && nrow(rows) == 5) {
    rows <- rows[c(1, 2, 4, 3, 5), ]
  }
  
  list(
    Estimates = setNames(rows$est, params),
    "Standard Errors" = setNames(rows[[se_col]], params),
    "P-values" = setNames(rows[[pval_col]], params),
    CI_lower = setNames(rows$ci.lower, params),
    CI_upper = setNames(rows$ci.upper, params)
  )
}

#' Check for proper solution with configurable strictness
#' @param table Data.frame. Parameter table from estimation method output
#' @param method Character. Type of method: "lsam", "lms", "qml", or "upi"
#' @param check_positive_se Logical. Check that all estimated SEs > 0 (default FALSE)
#' @param check_heywood Logical. Check for negative variances (default FALSE)
#' @param check_pd Logical. Check positive definiteness of factor cov matrix (default FALSE)
#' @return Logical. TRUE if passes all enabled checks, FALSE otherwise.
check_proper_solution <- function(table, method,
                                  check_positive_se = FALSE,
                                  check_heywood = FALSE,
                                  check_pd = FALSE) {
  
  # no additional checks requested, basic check done elsewhere
  if (!check_positive_se && !check_heywood && !check_pd) {
    return(TRUE)
  }
  
  if (is.null(table) || nrow(table) == 0) return(FALSE)
  
  se_col <- if (method %in% c("lms", "qml")) "std.error" else "se"
  
  if (!se_col %in% names(table)) return(FALSE)
  
  # check 1: positive SEs 
  if (check_positive_se) {
    se_values <- table[[se_col]]
    # dont count the fixed parameters (SE = 0 for LSAM/UPI, SE = NA for LMS/QML)
    estimable_rows <- !is.na(se_values) & se_values != 0
    
    if (sum(estimable_rows) == 0) return(FALSE)
    
    se_vals <- se_values[estimable_rows]
    if (any(se_vals < 0, na.rm = TRUE)) return(FALSE)
  }
  
  # check 2: positive variances (Heywood cases)
  if (check_heywood) {
    if (all(c("lhs", "rhs", "op") %in% names(table))) {
      var_rows <- table$lhs == table$rhs & table$op == "~~"
      if (any(var_rows)) {
        var_estimates <- table$est[var_rows]
        if (any(var_estimates < 0, na.rm = TRUE)) return(FALSE)
      }
    }
  }
  
  # check 3: positive definite factor covariance matrix 
  if (check_pd) {
    latents <- c("eta1", "eta2", "eta3")
    
    if (all(c("lhs", "rhs", "op") %in% names(table))) {
      cov_rows <- table[table$op == "~~" & 
                          table$lhs %in% latents & 
                          table$rhs %in% latents, ]
      
      if (nrow(cov_rows) > 0) {
        # build covariance matrix as we didn't store it from sim runs
        cov_matrix <- matrix(0, 3, 3, dimnames = list(latents, latents))
        
        for (i in seq_len(nrow(cov_rows))) {
          row_lhs <- cov_rows$lhs[i]
          row_rhs <- cov_rows$rhs[i]
          est <- cov_rows$est[i]
          
          cov_matrix[row_lhs, row_rhs] <- est
          cov_matrix[row_rhs, row_lhs] <- est
        }
        
        # positive definiteness through eig. values
        eigenvalues <- eigen(cov_matrix, only.values = TRUE)$values
        if (!all(eigenvalues > 0)) return(FALSE)
      }
    }
  }
  
  TRUE
}

#' Extract Convergence and Outliers with Configurable Checks
#' @param all_results List. Complete simulation results
#' @param parameters_of_interest Character vector. Parameters to extract
#' @param remove_outliers Logical. Whether to remove outliers (default TRUE)
#' @param outlier_threshold Numeric. IQR multiplier for outlier detection (default 3)
#' @param min_reps Integer. Minimum replications required (default 10)
#' @param exclude_warnings Logical. Whether to exclude warnings (default FALSE)
#' @param check_positive_se Logical. Check positive SEs (default FALSE)
#' @param check_heywood Logical. Check for Heywood cases (default FALSE)
#' @param check_pd Logical. Check positive definiteness (default FALSE)
#' 
#' @details Convergence Check Levels:
#'   - Basic (all FALSE): Only NA/NaN/Inf checks on parameters of interest
#'   - Strict (all TRUE): Adds positive SEs, no Heywood cases, PD factor cov matrix
ExtractConvergenceOutliers <- function(all_results,
                                       parameters_of_interest = c("eta1","eta2","eta1:eta2","eta1:eta1","eta2:eta2"),
                                       remove_outliers = TRUE,
                                       outlier_threshold = 3,
                                       min_reps = 10,
                                       exclude_warnings = FALSE,
                                       check_positive_se = FALSE,
                                       check_heywood = FALSE,
                                       check_pd = FALSE) {
  
  options(scipen = 999)
  methods <- c("lsam", "lms", "qml", "upi")
  param_index <- setNames(1:5, parameters_of_interest)
  
  # which checks are enabled
  cat("Convergence checks enabled:\n")
  cat("\nBasic (NA/NaN/Inf): TRUE (always)\n")
  cat("\nPositive SEs:", check_positive_se, "\n")
  cat("\nNo Heywood cases:", check_heywood, "\n")
  cat("\nPositive definite factor cov:", check_pd, "\n\n")
  
  summary_list <- list()
  convergence_outliers_details <- list()
  filtered_data <- list()
  converged_data <- list()
  
  for (i in seq_along(all_results)) {
    condition <- all_results[[i]]$condition
    res_list <- all_results[[i]]$results
    
    # condition variables 
    model_name <- ifelse(condition$Model_Type == "linear", "Linear", "Full")
    distr_exo <- as.character(condition$Distr_Exo)
    distr_epsilon <- as.character(condition$Distr_Epsilon)
    distr_zeta <- as.character(condition$Distr_Zeta)
    sample_size <- condition$N
    reliability <- condition$Rel
    
    true_values <- all_results[[i]]$true_parameters %||% 
      if (condition$Model_Type == "linear") c(0.316, 0.316, 0, 0, 0) 
    else c(0.316, 0.316, 0.139, 0.101, 0.101)
    
    max_reps <- max(sapply(methods, function(m) 
      length(res_list[[paste0(m, "_tables")]]) %||% 0))
    
    if (max_reps == 0) next
    
    # check convergence with the configurable checks
    method_convergence <- setNames(lapply(methods, function(method) {
      tbls <- res_list[[paste0(method, "_tables")]]
      if (is.null(tbls)) return(rep(FALSE, max_reps))
      
      vapply(seq_len(max_reps), function(rep) {
        if (rep > length(tbls) || is.null(tbls[[rep]])) return(FALSE)
        
        # additional checks (if any enabled)
        if (!check_proper_solution(tbls[[rep]], method,
                                   check_positive_se = check_positive_se,
                                   check_heywood = check_heywood,
                                   check_pd = check_pd)) {
          return(FALSE)
        }
        
        # basic check: parameters of interest have valid values
        extr <- extract_params_study_1(tbls[[rep]], method)
        if (is.null(extr)) return(FALSE)
        
        all(sapply(parameters_of_interest, function(p) {
          p %in% names(extr$Estimates) && 
            is_valid(c(extr$Estimates[p], extr$`Standard Errors`[p], 
                       extr$`P-values`[p], extr$CI_lower[p], extr$CI_upper[p]))
        }))
      }, logical(1))
    }), methods)
    
    # warnings processing
    method_warnings <- setNames(lapply(methods, function(method) {
      warns <- res_list$warnings[[method]]
      if (is.null(warns)) return(list())
      
      conv_reps <- which(method_convergence[[method]])
      warn_list <- list()
      for (rep in conv_reps) {
        if (rep <= length(warns) && !is.null(warns[[rep]]) && length(warns[[rep]]) > 0) {
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
    
    method_convergence_counts <- sapply(method_convergence, sum)
    method_convergence_rates <- sapply(method_convergence, function(x) mean(x) * 100)
    valid_reps <- Reduce("&", method_convergence)
    
    n_excluded_by_warnings <- 0
    if (exclude_warnings) {
      warn_reps <- unique(unlist(lapply(method_warnings, function(w) as.numeric(names(w)))))
      if (length(warn_reps) > 0) {
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
    
    # store converged data using helper
    for (method in methods) {
      tbls <- res_list[[paste0(method, "_tables")]]
      if (is.null(tbls)) next
      
      for (p_name in parameters_of_interest) {
        converged_list <- extract_data_list(tbls, valid_rep_indices, method, p_name)
        
        if (length(converged_list) > 0) {
          key <- paste(i, toupper(method), p_name, sep = "_")
          converged_data[[key]] <- list(
            condition = i, method = toupper(method), model = model_name,
            distr_exo = distr_exo, distr_epsilon = distr_epsilon, distr_zeta = distr_zeta,
            sample_size = sample_size, reliability = reliability,
            parameter = p_name, true_value = true_values[param_index[p_name]],
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
      
      param_ests <- setNames(lapply(parameters_of_interest, function(p) {
        unlist(lapply(valid_rep_indices, function(rep) {
          extr <- extract_params_study_1(tbls[[rep]], method)
          if (!is.null(extr) && p %in% names(extr$Estimates)) extr$Estimates[p]
        }))
      }), parameters_of_interest)
      
      bounds <- lapply(param_ests, function(ests) {
        if (length(ests) == 0) return(c(-Inf, Inf))
        q <- quantile(ests, c(0.25, 0.75), na.rm = TRUE)
        iqr <- diff(q)
        c(q[1] - outlier_threshold * iqr, q[2] + outlier_threshold * iqr)
      })
      
      for (rep in valid_rep_indices) {
        extr <- extract_params_study_1(tbls[[rep]], method)
        if (!is.null(extr)) {
          outliers[rep] <- any(sapply(parameters_of_interest, function(p) {
            p %in% names(extr$Estimates) && is_valid(extr$Estimates[p]) &&
              (extr$Estimates[p] < bounds[[p]][1] || extr$Estimates[p] > bounds[[p]][2])
          }))
        }
      }
      outliers
    }), methods)
    
    method_outlier_counts <- sapply(method_outliers, function(x) sum(x[valid_rep_indices]))
    convergence_outliers_details[[paste0("Condition_", i)]]$outlier_iterations <- 
      lapply(method_outliers, which)
    
    valid_after_outliers <- valid_reps & !Reduce("|", method_outliers)
    valid_final_indices <- which(valid_after_outliers)
    n_after_global_outlier <- length(valid_final_indices)
    
    if (n_after_global_outlier < min_reps) {
      message(sprintf("Condition %d: Only %d replications after outliers. Skipping.", 
                      i, n_after_global_outlier))
      next
    }
    
    # extract filtered data and build summary using helper
    for (method in methods) {
      tbls <- res_list[[paste0(method, "_tables")]]
      if (is.null(tbls)) next
      
      for (p_name in parameters_of_interest) {
        data_list <- extract_data_list(tbls, valid_final_indices, method, p_name)
        if (length(data_list) == 0) next
        
        tv <- true_values[param_index[p_name]]
        key <- paste(i, toupper(method), p_name, sep = "_")
        
        # accumulate summary rows in list
        summary_list[[length(summary_list) + 1]] <- data.frame(
          Condition = i, Method = toupper(method), Model = model_name,
          Distr_Exo = distr_exo, Distr_Epsilon = distr_epsilon, Distr_Zeta = distr_zeta,
          SampleSize = sample_size, Reliability = reliability,
          Parameter = p_name, TrueValue = tv, N_Total = max_reps,
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
          Percent_Total_Excluded = round((max_reps - length(data_list)) / max_reps * 100, 2),
          stringsAsFactors = FALSE
        )
        
        filtered_data[[key]] <- list(
          condition = i, method = toupper(method), model = model_name,
          distr_exo = distr_exo, distr_epsilon = distr_epsilon, distr_zeta = distr_zeta,
          sample_size = sample_size, reliability = reliability,
          parameter = p_name, true_value = tv,
          estimates = sapply(data_list, "[[", "est"),
          standard_errors = sapply(data_list, "[[", "se"),
          p_values = sapply(data_list, "[[", "pv"),
          ci_lower = sapply(data_list, "[[", "lo"),
          ci_upper = sapply(data_list, "[[", "hi")
        )
      }
    }
  }
  
  # combine all summary rows at once
  convergence_outliers_summary <- do.call(rbind, summary_list)
  
  list(
    convergence_outliers_summary = convergence_outliers_summary,
    convergence_outliers_details = convergence_outliers_details,
    filtered_data = filtered_data, # for measures without outliers
    converged_data = converged_data # for measures with outliers
  )
}

CalculatePerformanceMetrics <- function(filtered_data,
                                        converged_data,
                                        convergence_outliers_summary,
                                        alpha = 0.05) {
  
  results_list <- list()
  
  for (key in names(filtered_data)) {
    d <- filtered_data[[key]]
    tv <- d$true_value
    
    df <- data.frame(
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
    
    if (tv != 0) {
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
    
    # robust metrics from converged data
    d_conv <- converged_data[[key]] %||% d
    df_conv <- data.frame(est = d_conv$estimates, se = d_conv$standard_errors)
    
    median_est <- median(df_conv$est)
    bias_median <- median_est - tv
    mad_est <- mad(df_conv$est)
    trimmed_est <- mean(df_conv$est, trim = 0.20)
    bias_trimmed <- trimmed_est - tv
    rmse_trimmed <- sqrt(bias_trimmed^2 + mad_est^2)
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
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Condition = d$condition, Method = d$method, Model = d$model,
      Distr_Exo = d$distr_exo, Distr_Epsilon = d$distr_epsilon, Distr_Zeta = d$distr_zeta,
      SampleSize = d$sample_size, Reliability = d$reliability,
      Parameter = d$parameter, TrueValue = tv,
      MeanEstimate = mean(df$est), MedianEstimate = median_est, TrimmedEstimate = trimmed_est,
      Bias_Mean = abs_metrics$bias, Bias_Median = bias_median, Bias_Trimmed = bias_trimmed,
      RelativeBias_Mean = rel_bias_mean,
      RelativeBias_Median = if (tv != 0) bias_median / tv else NA,
      PercentRelativeBias_Mean = rel_bias_mean_pct,
      PercentRelativeBias_Median = if (tv != 0) bias_median / tv * 100 else NA,
      MAD = mad_est, Variance = abs_metrics$var, SD = abs_metrics$stddev,
      MSE_Mean = abs_metrics$mse, RMSE_Mean = abs_metrics$rmse,
      RMSE_Median = sqrt(bias_median^2 + mad_est^2), RMSE_Trimmed = rmse_trimmed,
      Relative_MSE = rel_metrics$rel_mse %||% NA, Relative_RMSE = rel_metrics$rel_rmse %||% NA,
      MeanSE = mean(df$se), TrimmedSE = trimmed_se,
      SE_SD_Ratio = mean(df$se) / abs_metrics$stddev, RSB = rsb,
      CoverageRate = cov_metrics$coverage * 100, CI_Width = cov_metrics$width,
      RejectionRate = rej_metrics$rej_rate * 100,
      TypeI_Error = if (tv == 0) rej_metrics$rej_rate * 100 else NA,
      Power = if (tv != 0) rej_metrics$rej_rate * 100 else NA,
      Bias_Mean_MCSE = abs_metrics$bias_mcse, RelativeBias_MCSE = rel_metrics$rel_bias_mcse %||% NA,
      Variance_MCSE = abs_metrics$var_mcse, SD_MCSE = abs_metrics$stddev_mcse,
      MSE_Mean_MCSE = abs_metrics$mse_mcse, RMSE_Mean_MCSE = abs_metrics$rmse_mcse,
      Relative_MSE_MCSE = rel_metrics$rel_mse_mcse %||% NA,
      Relative_RMSE_MCSE = rel_metrics$rel_rmse_mcse %||% NA,
      CoverageRate_MCSE = cov_metrics$coverage_mcse, CI_Width_MCSE = cov_metrics$width_mcse,
      RejectionRate_MCSE = rej_metrics$rej_rate_mcse,
      TypeI_Error_MCSE = if (tv == 0) rej_metrics$rej_rate_mcse else NA,
      Power_MCSE = if (tv != 0) rej_metrics$rej_rate_mcse else NA,
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
      Percent_Total_Excluded = conv_row$Percent_Total_Excluded,
      stringsAsFactors = FALSE
    )
  }
  
  do.call(rbind, results_list)
}


# ============================================================================
# USAGE EXAMPLES
# ============================================================================

# --- BASIC CHECKS (for main paper) ---
# extraction_basic <- ExtractConvergenceOutliers(
#   all_results = all_results,
#   parameters_of_interest = c("eta1", "eta2", "eta1:eta2", "eta1:eta1", "eta2:eta2"),
#   remove_outliers = TRUE,
#   outlier_threshold = 3,
#   min_reps = 0,
#   check_positive_se = FALSE,
#   check_heywood = FALSE,
#   check_pd = FALSE
# )
# 
# results_basic <- CalculatePerformanceMetrics(
#   filtered_data = extraction_basic$filtered_data,
#   converged_data = extraction_basic$converged_data,
#   convergence_outliers_summary = extraction_basic$convergence_outliers_summary
# )

# --- STRICT CHECKS (for supplementary materials) ---
# extraction_strict <- ExtractConvergenceOutliers(
#   all_results = all_results,
#   parameters_of_interest = c("eta1", "eta2", "eta1:eta2", "eta1:eta1", "eta2:eta2"),
#   remove_outliers = TRUE,
#   outlier_threshold = 3,
#   min_reps = 0,
#   exclude_warnings = TRUE, # note this also
#   check_positive_se = TRUE,
#   check_heywood = TRUE,
#   check_pd = TRUE
# )
# 
# results_strict <- CalculatePerformanceMetrics(
#   filtered_data = extraction_strict$filtered_data,
#   converged_data = extraction_strict$converged_data,
#   convergence_outliers_summary = extraction_strict$convergence_outliers_summary
# )