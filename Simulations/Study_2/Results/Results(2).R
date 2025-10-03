############################### Study 2 Results Processing (5-Factor Model) ###############################

# helper function for Study 2 parameter extraction with QML ordering fix
extract_study2_parameters <- function(table, method_type, equation) {
  if(is.null(table)) return(NULL)
  
  if(method_type == "dblcent") {
    rows <- table[table$lhs == equation & table$op == "~", ]
    
    # map parameters based on equation (UPDATED FOR NEW STRUCTURE)
    if(equation == "eta4") {
      params <- if(nrow(rows) == 3) {
        c("eta1", "eta2", "eta3")
      } else {
        c("eta1", "eta2", "eta3", "eta1:eta2", "eta1:eta3", "eta1:eta1", "eta2:eta2")
      }
    } else if(equation == "eta5") {
      params <- if(nrow(rows) == 4) {
        c("eta4", "eta1", "eta2", "eta3")
      } else {
        c("eta4", "eta1", "eta2", "eta3", "eta1:eta4", "eta2:eta4", "eta1:eta1", "eta3:eta3")
      }
    }
    
    RESULTS <- list(
      "Estimates" = setNames(rows$est, params),
      "Standard Errors" = setNames(rows$se, params),
      "P-values" = setNames(rows$pvalue, params),
      "CI_lower" = setNames(rows$ci.lower, params),
      "CI_upper" = setNames(rows$ci.upper, params)
    )
    
  } else if(method_type == "sam") {
    coefs <- table[table$lhs == equation & table$op == "~", ]
    
    if(equation == "eta4") {
      params <- if(nrow(coefs) == 3) {
        c("eta1", "eta2", "eta3")
      } else {
        c("eta1", "eta2", "eta3", "eta1:eta2", "eta1:eta3", "eta1:eta1", "eta2:eta2")
      }
    } else if(equation == "eta5") {
      params <- if(nrow(coefs) == 4) {
        c("eta4", "eta1", "eta2", "eta3")
      } else {
        c("eta4", "eta1", "eta2", "eta3", "eta1:eta4", "eta2:eta4", "eta1:eta1", "eta3:eta3")
      }
    }
    
    rows <- table$lhs == equation & table$op == "~" & table$rhs %in% params
    
    ests <- setNames(table$est[rows], table$rhs[rows])
    ses  <- setNames(table$se[rows], table$rhs[rows])
    pval <- setNames(table$pvalue[rows], table$rhs[rows])
    ci_lower <- setNames(table$ci.lower[rows], table$rhs[rows])
    ci_upper <- setNames(table$ci.upper[rows], table$rhs[rows])
    
    RESULTS <- list(
      "Estimates"        = ests[params],
      "Standard Errors"  = ses[params],
      "P-values"         = pval[params],
      "CI_lower"         = ci_lower[params],
      "CI_upper"         = ci_upper[params]
    )
    
  } else if(method_type == "qml") {
    rows <- table[table$lhs == equation & table$op == "~", ]
    
    # QML specific ordering issues (UPDATED FOR NEW INTERACTIONS)
    if(equation == "eta4" && nrow(rows) == 7) {
      # For eta4 with new interactions: eta1:eta2, eta1:eta3, eta1:eta1, eta2:eta2
      # Check if any ordering issues exist and handle them
      RESULTS <- list(
        "Estimates" = setNames(rows$est, rows$rhs),
        "Standard Errors" = setNames(rows$std.error, rows$rhs),
        "P-values" = setNames(rows$p.value, rows$rhs),
        "CI_lower" = setNames(rows$ci.lower, rows$rhs),
        "CI_upper" = setNames(rows$ci.upper, rows$rhs)
      )
    } else if(equation == "eta5" && nrow(rows) == 8) {
      # For eta5 with new interactions: eta1:eta4, eta2:eta4, eta1:eta1, eta3:eta3
      # QML might have ordering issues with quadratics vs interactions
      
      # Check for potential swaps needed
      idx <- 1:nrow(rows)
      
      # Find positions of interactions and quadratics
      eta1eta4_idx <- which(rows$rhs == "eta1:eta4")
      eta2eta4_idx <- which(rows$rhs == "eta2:eta4")
      eta1eta1_idx <- which(rows$rhs == "eta1:eta1")
      eta3eta3_idx <- which(rows$rhs == "eta3:eta3")
      
      # Apply any necessary reordering here if QML outputs in wrong order
      # For now, keep original order unless specific issues are identified
      
      RESULTS <- list(
        "Estimates" = setNames(rows$est[idx], rows$rhs[idx]),
        "Standard Errors" = setNames(rows$std.error[idx], rows$rhs[idx]),
        "P-values" = setNames(rows$p.value[idx], rows$rhs[idx]),
        "CI_lower" = setNames(rows$ci.lower[idx], rows$rhs[idx]),
        "CI_upper" = setNames(rows$ci.upper[idx], rows$rhs[idx])
      )
    } else {
      # cases without interaction/quadratic terms, no reordering needed
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


ExtractConvergenceOutliers_Study2 <- function(all_results,
                                              parameters_of_interest = list(
                                                eta4 = c("eta1","eta2","eta3","eta1:eta2","eta1:eta3","eta1:eta1","eta2:eta2"),
                                                eta5 = c("eta4","eta1","eta2","eta3","eta1:eta4","eta2:eta4","eta1:eta1","eta3:eta3")
                                              ),
                                              remove_outliers = TRUE,
                                              outlier_threshold = 3,
                                              min_reps = 10,
                                              exclude_warnings = FALSE) {
  
  options(scipen = 999)

  methods <- c("sam", "qml", "dblcent") 
  convergence_outliers_summary <- dplyr::tibble()
  convergence_outliers_details <- list()
  filtered_data <- list()
  
  for (i in seq_along(all_results)) {
    condition <- all_results[[i]]$condition
    res_list  <- all_results[[i]]$results
    
    # true values from stored parameters
    true_params <- all_results[[i]]$true_parameters
    
    model_name <- ifelse(condition$Model_Type == "linear", "Linear", "Full")
    distribution <- as.character(condition$Distribution)
    sample_size <- condition$N
    reliability <- condition$Rel
    
    # max replications
    max_reps <- 0
    for (method in methods) {
      tbls <- res_list[[paste0(method, "_tables")]]
      if (!is.null(tbls)) {
        max_reps <- max(max_reps, length(tbls))
      }
    }
    
    if (max_reps == 0) next
    
    # track convergence for each method
    method_convergence <- list()
    for (method in methods) {
      method_convergence[[method]] <- rep(FALSE, max_reps)  
    }
    
    # check convergence for each method and equation (NOW ONLY eta4 and eta5)
    for (method in methods) {
      tbls <- res_list[[paste0(method, "_tables")]]
      if (is.null(tbls)) next
      
      for (rep in seq_len(max_reps)) {
        if (rep > length(tbls) || is.null(tbls[[rep]])) next
        
        converged <- TRUE
        # check both equations (eta4 and eta5 only)
        for (eq in names(parameters_of_interest)) {
          extr <- extract_study2_parameters(tbls[[rep]], method, eq)
          if (is.null(extr)) {
            converged <- FALSE
            break
          }
          
          # check each parameter in this equation
          for (p in parameters_of_interest[[eq]]) {
            if (p %in% names(extr$Estimates)) {
              if (is.na(extr$Estimates[p]) || is.nan(extr$Estimates[p]) || is.infinite(extr$Estimates[p]) ||
                  is.na(extr$`Standard Errors`[p]) || is.nan(extr$`Standard Errors`[p]) || is.infinite(extr$`Standard Errors`[p]) ||
                  is.na(extr$`P-values`[p]) || is.nan(extr$`P-values`[p]) || is.infinite(extr$`P-values`[p]) ||
                  is.na(extr$CI_lower[p]) || is.nan(extr$CI_lower[p]) || is.infinite(extr$CI_lower[p]) ||
                  is.na(extr$CI_upper[p]) || is.nan(extr$CI_upper[p]) || is.infinite(extr$CI_upper[p])) {
                converged <- FALSE
                break
              }
            } else if (!grepl(":", p)) {  # main effects should always be present
              converged <- FALSE
              break
            }
          }
          if (!converged) break
        }
        method_convergence[[method]][rep] <- converged
      }
    }
    
    # track warnings for converged iterations
    method_warnings <- list()
    method_warning_counts <- list()
    
    for (method in methods) {
      method_warnings[[method]] <- list()
      method_warning_counts[[method]] <- 0
      
      method_warns <- res_list[["warnings"]][[method]]
      if (!is.null(method_warns)) {
        for (rep in which(method_convergence[[method]])) {
          if (rep <= length(method_warns) && !is.null(method_warns[[rep]]) && length(method_warns[[rep]]) > 0) {
            method_warnings[[method]][[as.character(rep)]] <- method_warns[[rep]]
            method_warning_counts[[method]] <- method_warning_counts[[method]] + 1
          }
        }
      }
    }
    
    # store non-converged iterations
    non_converged_iterations <- list()
    for (method in methods) {
      non_converged_iterations[[method]] <- which(!method_convergence[[method]])
    }
    
    convergence_outliers_details[[paste0("Condition_", i)]] <- list(
      non_converged_iterations = non_converged_iterations,
      warning_iterations = method_warnings,
      warning_counts = method_warning_counts,
      outlier_iterations = list()
    )
    
    # individual method convergence statistics
    method_convergence_counts <- list()
    method_convergence_rates <- list()
    
    for (method in methods) {
      method_convergence_counts[[method]] <- sum(method_convergence[[method]])
      method_convergence_rates[[method]] <- sum(method_convergence[[method]]) / max_reps * 100
    }
    
    # valid replications across all methods
    valid_reps <- rep(TRUE, max_reps)
    for (method in methods) {
      valid_reps <- valid_reps & method_convergence[[method]]
    }
    
    # optionally exclude warnings
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
    
    # track outliers for each method independently
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
        
        # collect all estimates from valid replications
        all_estimates <- list()
        for (eq in names(parameters_of_interest)) {
          for (p in parameters_of_interest[[eq]]) {
            all_estimates[[paste0(eq, "_", p)]] <- numeric()
          }
        }
        
        for (rep in valid_rep_indices) {
          for (eq in names(parameters_of_interest)) {
            extr <- extract_study2_parameters(tbls[[rep]], method, eq)
            if (!is.null(extr)) {
              for (p in names(extr$Estimates)) {
                key <- paste0(eq, "_", p)
                if (key %in% names(all_estimates)) {
                  all_estimates[[key]] <- c(all_estimates[[key]], extr$Estimates[p])
                }
              }
            }
          }
        }
        
        # calculate outlier bounds per parameter
        outlier_bounds <- list()
        for (key in names(all_estimates)) {
          if (length(all_estimates[[key]]) > 0) {
            q1 <- stats::quantile(all_estimates[[key]], 0.25, na.rm = TRUE)
            q3 <- stats::quantile(all_estimates[[key]], 0.75, na.rm = TRUE)
            iqr <- q3 - q1
            outlier_bounds[[key]] <- list(
              lower = q1 - outlier_threshold * iqr,
              upper = q3 + outlier_threshold * iqr
            )
          }
        }
        
        # mark outliers for each replication
        for (rep in valid_rep_indices) {
          is_outlier <- FALSE
          for (eq in names(parameters_of_interest)) {
            extr <- extract_study2_parameters(tbls[[rep]], method, eq)
            if (!is.null(extr)) {
              for (p in parameters_of_interest[[eq]]) {
                if (p %in% names(extr$Estimates)) {
                  key <- paste0(eq, "_", p)
                  if (key %in% names(outlier_bounds)) {
                    est_val <- extr$Estimates[p]
                    if (!is.na(est_val) && !is.nan(est_val) && !is.infinite(est_val)) {
                      if (est_val < outlier_bounds[[key]]$lower || 
                          est_val > outlier_bounds[[key]]$upper) {
                        is_outlier <- TRUE
                        break
                      }
                    }
                  }
                }
              }
              if (is_outlier) break
            }
          }
          method_outliers[[method]][rep] <- is_outlier
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
    
    # global outlier synchronization
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
    
    # extract filtered data for each method and equation
    for (method in methods) {
      tbls <- res_list[[paste0(method, "_tables")]]
      if (is.null(tbls)) next
      
      n_converged_this_method <- method_convergence_counts[[method]]
      convergence_rate_this_method <- method_convergence_rates[[method]]
      n_outliers_this_method <- method_outlier_counts[[method]]
      outlier_rate_this_method <- round(n_outliers_this_method / n_after_global_convergence * 100, 2)
      
      # process each equation
      for (eq in names(parameters_of_interest)) {
        for (p_name in parameters_of_interest[[eq]]) {
          # collect data from globally valid replications
          ests <- numeric()
          ses <- numeric()
          pvs <- numeric()
          los <- numeric()
          his <- numeric()
          
          for (rep in valid_final_indices) {
            tab <- tbls[[rep]]
            if (is.null(tab)) next
            extr <- extract_study2_parameters(tab, method, eq)
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
          
          # get true value (UPDATED FOR NEW PARAMETER NAMES)
          param_key <- paste0(eq, "_", gsub(":", "", p_name))
          tv <- if (!is.null(true_params[[param_key]])) {
            true_params[[param_key]]
          } else {
            0  # for interaction/quadratic terms in linear model
          }
          
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
              Equation = eq,
              Parameter = p_name,
              TrueValue = tv,
              
              # Sample size tracking
              N_Total = n_total,
              N_Converged_This_Method = n_converged_this_method,
              N_Warnings_This_Method = method_warning_counts[[method]],
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
          key <- paste(i, toupper(method), eq, p_name, sep = "_")
          filtered_data[[key]] <- list(
            condition = i,
            method = toupper(method),
            model = model_name,
            distribution = distribution,
            sample_size = sample_size,
            reliability = reliability,
            equation = eq,
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
  }
  
  list(
    convergence_outliers_summary = convergence_outliers_summary,
    convergence_outliers_details = convergence_outliers_details,
    filtered_data = filtered_data
  )
}

# performance Metrics for Study 2 (unchanged)
CalculatePerformanceMetrics_Study2 <- function(filtered_data, 
                                               convergence_outliers_summary,
                                               alpha = 0.05) {
  
  library(simhelpers)
  results_summary <- dplyr::tibble()
  
  for (key in names(filtered_data)) {
    data_item <- filtered_data[[key]]
    
    condition_num <- data_item$condition
    method_name <- data_item$method
    model_name <- data_item$model
    distribution <- data_item$distribution
    sample_size <- data_item$sample_size
    reliability <- data_item$reliability
    eq <- data_item$equation
    p_name <- data_item$parameter
    tv <- data_item$true_value
    
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
    abs_metrics <- calc_absolute(
      data = df,
      estimates = est, 
      true_param = true_param,
      criteria = c("bias", "variance", "stddev", "mse", "rmse"),
      winz = Inf
    )
    
    # RELATIVE METRICS (only if true value is non-zero)
    if (tv != 0) {
      rel_metrics <- calc_relative(
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
    cov_metrics <- calc_coverage(
      data = df,
      lower_bound = lower_bound, 
      upper_bound = upper_bound, 
      true_param = true_param,
      criteria = c("coverage", "width"),
      winz = Inf
    )
    
    # REJECTION RATES
    rej_metrics <- calc_rejection(
      data = df,
      p_values = p_val,
      alpha = alpha
    )
    
    rejection_rate <- rej_metrics$rej_rate * 100
    rej_mcse <- rej_metrics$rej_rate_mcse
    
    # type I Error vs Power
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
        convergence_outliers_summary$Equation == eq &
        convergence_outliers_summary$Parameter == p_name, ]
    
    if (nrow(conv_outlier_row) == 0) {
      warning(paste("No convergence/outlier info found for", 
                    condition_num, method_name, eq, p_name))
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
        Equation = eq,
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

# wrapper function for Study 2 (UPDATED PARAMETERS)
CalculatePerformance_Study2 <- function(all_results,
                                        parameters_of_interest = list(
                                          eta4 = c("eta1","eta2","eta3","eta1:eta2","eta1:eta3","eta1:eta1","eta2:eta2"),
                                          eta5 = c("eta4","eta1","eta2","eta3","eta1:eta4","eta2:eta4","eta1:eta1","eta3:eta3")
                                        ),
                                        remove_outliers = TRUE,
                                        outlier_threshold = 3,
                                        alpha = 0.05,
                                        min_reps = 10,
                                        exclude_warnings = FALSE,
                                        return_convergence_details = TRUE) {
  
  # convergence and outliers
  extraction_results <- ExtractConvergenceOutliers_Study2(
    all_results = all_results,
    parameters_of_interest = parameters_of_interest,
    remove_outliers = remove_outliers,
    outlier_threshold = outlier_threshold,
    min_reps = min_reps,
    exclude_warnings = exclude_warnings
  )
  
  # performance metrics
  performance_results <- CalculatePerformanceMetrics_Study2(
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

# RUN THE ANALYSIS FOR STUDY 2 (UPDATED PARAMETERS)
#load("Simulations/Study_2/Data/Results_Study_2_final.RData")

# results with the new 5-factor structure
#results_study_2 <- CalculatePerformance_Study2(
#  all_results,
#  parameters_of_interest = list(
#    eta4 = c("eta1", "eta2", "eta3", "eta1:eta2", "eta1:eta3", "eta1:eta1", "eta2:eta2"),
#    eta5 = c("eta4", "eta1", "eta2", "eta3", "eta1:eta4", "eta2:eta4", "eta1:eta1", "eta3:eta3")
#  ),
#  remove_outliers = TRUE,
#  outlier_threshold = 3,
#  alpha = 0.05,
#  min_reps = 10,
#  exclude_warnings = FALSE,
#  return_convergence_details = TRUE
#)