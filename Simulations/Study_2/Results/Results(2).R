# parameter extraction 
extract_study2_params <- function(table, method, equation) {
  if(is.null(table)) return(NULL)
  
  rows <- table[table$lhs == equation & table$op == "~", ]
  
  # expected parameter order based on sam 
  params <- if(equation == "eta4") {
    if(nrow(rows) == 3) c("eta1", "eta2", "eta3")
    else c("eta1", "eta2", "eta3", "eta1:eta2", "eta1:eta3", "eta1:eta1", "eta2:eta2")
  } else { # eta5
    if(nrow(rows) == 4) c("eta4", "eta1", "eta2", "eta3")
    else c("eta4", "eta1", "eta2", "eta3", "eta1:eta4", "eta2:eta4", "eta1:eta1", "eta3:eta3")
  }
  
  se_col <- if(method == "qml") "std.error" else "se"
  pval_col <- if(method == "qml") "p.value" else "pvalue"
  
  if(method == "qml") {
    matched_rows <- match(params, rows$rhs)
    rows <- rows[matched_rows, ]
  }
  
  list(
    Estimates = setNames(rows$est, params),
    "Standard Errors" = setNames(rows[[se_col]], params),
    "P-values" = setNames(rows[[pval_col]], params),
    CI_lower = setNames(rows$ci.lower, params),
    CI_upper = setNames(rows$ci.upper, params)
  )
}

# main extraction function for sim 2
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
  
  `%||%` <- function(x, y) if(is.null(x)) y else x
  
  convergence_outliers_summary <- dplyr::tibble()
  convergence_outliers_details <- list()
  filtered_data <- list()
  
  for (i in seq_along(all_results)) {
    condition <- all_results[[i]]$condition
    res_list <- all_results[[i]]$results
    true_params <- all_results[[i]]$true_parameters
    
    # max replications
    max_reps <- max(sapply(methods, function(m) 
      length(res_list[[paste0(m, "_tables")]]) %||% 0))
    
    if (max_reps == 0) next
    
    # check convergence 
    method_convergence <- setNames(lapply(methods, function(method) {
      tbls <- res_list[[paste0(method, "_tables")]]
      if (is.null(tbls)) return(rep(FALSE, max_reps))
      
      vapply(seq_len(max_reps), function(rep) {
        if (rep > length(tbls) || is.null(tbls[[rep]])) return(FALSE)
        
        # check both equations (!)
        all(sapply(names(parameters_of_interest), function(eq) {
          extr <- extract_study2_params(tbls[[rep]], method, eq)
          if (is.null(extr)) return(FALSE)
          
          all(sapply(parameters_of_interest[[eq]], function(p) {
            if (!(p %in% names(extr$Estimates))) {
              return(!grepl(":", p))  # main effects must exist, interactions "optional"
            }
            vals <- c(extr$Estimates[p], extr$`Standard Errors`[p], 
                      extr$`P-values`[p], extr$CI_lower[p], extr$CI_upper[p])
            !any(is.na(vals) | is.nan(vals) | is.infinite(vals))
          }))
        }))
      }, logical(1))
    }), methods)
    
    # warnings
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
    
    # convergence details
    convergence_outliers_details[[paste0("Condition_", i)]] <- list(
      non_converged_iterations = lapply(method_convergence, function(x) which(!x)),
      warning_iterations = method_warnings,
      warning_counts = method_warning_counts,
      outlier_iterations = list()
    )
    
    method_convergence_counts <- sapply(method_convergence, sum)
    method_convergence_rates <- sapply(method_convergence, function(x) mean(x) * 100)
    
    # global valid reps
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
    
    # outlier 
    method_outliers <- setNames(lapply(methods, function(method) {
      outliers <- rep(FALSE, max_reps)
      if (!remove_outliers) return(outliers)
      
      tbls <- res_list[[paste0(method, "_tables")]]
      if (is.null(tbls)) return(outliers)
      
      # collect estimates by equation and parameter
      all_estimates <- list()
      for(eq in names(parameters_of_interest)) {
        for(p in parameters_of_interest[[eq]]) {
          key <- paste0(eq, "_", p)
          all_estimates[[key]] <- unlist(lapply(valid_rep_indices, function(rep) {
            extr <- extract_study2_params(tbls[[rep]], method, eq)
            if(!is.null(extr) && p %in% names(extr$Estimates)) extr$Estimates[p]
          }))
        }
      }
      
      # same as lonati et al. 2024 again
      bounds <- lapply(all_estimates, function(ests) {
        if(length(ests) == 0) return(c(-Inf, Inf))
        q <- quantile(ests, c(0.25, 0.75), na.rm = TRUE)
        iqr <- diff(q)
        c(q[1] - outlier_threshold * iqr, q[2] + outlier_threshold * iqr)
      })
      
      # flag all outliers
      for(rep in valid_rep_indices) {
        is_outlier <- FALSE
        for(eq in names(parameters_of_interest)) {
          extr <- extract_study2_params(tbls[[rep]], method, eq)
          if(!is.null(extr)) {
            for(p in parameters_of_interest[[eq]]) {
              if(p %in% names(extr$Estimates)) {
                key <- paste0(eq, "_", p)
                est_val <- extr$Estimates[p]
                if(!is.na(est_val) && !is.nan(est_val) && !is.infinite(est_val) &&
                   (est_val < bounds[[key]][1] || est_val > bounds[[key]][2])) {
                  is_outlier <- TRUE
                  break
                }
              }
            }
          }
          if(is_outlier) break
        }
        outliers[rep] <- is_outlier
      }
      outliers
    }), methods)
    
    method_outlier_counts <- sapply(method_outliers, function(x) sum(x[valid_rep_indices]))
    convergence_outliers_details[[paste0("Condition_", i)]]$outlier_iterations <- 
      lapply(method_outliers, which)
    
    # fianl valid indices
    valid_after_outliers <- valid_reps & !Reduce("|", method_outliers)
    valid_final_indices <- which(valid_after_outliers)
    n_after_global_outlier <- length(valid_final_indices)
    
    if (n_after_global_outlier < min_reps) {
      message(sprintf("Condition %d: Only %d replications after outliers. Skipping.", 
                      i, n_after_global_outlier))
      next
    }
    
    # data and build summary
    model_name <- ifelse(condition$Model_Type == "linear", "Linear", "Full")
    distribution <- as.character(condition$Distribution)
    sample_size <- condition$N
    reliability <- condition$Rel
    
    for(method in methods) {
      tbls <- res_list[[paste0(method, "_tables")]]
      if (is.null(tbls)) next
      
      for(eq in names(parameters_of_interest)) {
        for(p_name in parameters_of_interest[[eq]]) {
          # extract data
          data_list <- lapply(valid_final_indices, function(rep) {
            extr <- extract_study2_params(tbls[[rep]], method, eq)
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
          
          # true values
          param_key <- paste0(eq, "_", gsub(":", "", p_name))
          tv <- true_params[[param_key]] %||% 0
          
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
              Equation = eq,
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
            estimates = sapply(data_list, "[[", "est"),
            standard_errors = sapply(data_list, "[[", "se"),
            p_values = sapply(data_list, "[[", "pv"),
            ci_lower = sapply(data_list, "[[", "lo"),
            ci_upper = sapply(data_list, "[[", "hi")
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

# performance metrics calculation for Study 2
CalculatePerformanceMetrics_Study2 <- function(filtered_data, 
                                               convergence_outliers_summary,
                                               alpha = 0.05) {
  
  library(simhelpers)
  `%||%` <- function(x, y) if(is.null(x)) y else x
  
  results_summary <- dplyr::tibble()
  
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
    
    # metrics with simhelpers
    abs_metrics <- calc_absolute(df, est, true_param, 
                                 c("bias", "variance", "stddev", "mse", "rmse"), winz = Inf)
    cov_metrics <- calc_coverage(df, lower_bound, upper_bound, true_param, 
                                 c("coverage", "width"), winz = Inf)
    rej_metrics <- calc_rejection(df, p_val, alpha)
    
    # relative metrics
    if(tv != 0) {
      rel_metrics <- calc_relative(df, est, true_param, 
                                   c("relative bias", "relative mse", "relative rmse"), winz = Inf)
      rel_bias_mean <- rel_metrics$rel_bias - 1
      rel_bias_mean_pct <- rel_bias_mean * 100
    } else {
      rel_metrics <- list(rel_mse = NA, rel_rmse = NA,
                          rel_bias_mcse = NA, rel_mse_mcse = NA, rel_rmse_mcse = NA)
      rel_bias_mean <- NA
      rel_bias_mean_pct <- NA
    }
    
    # nedian metrics
    median_est <- median(df$est)
    bias_median <- median_est - tv
    mad_est <- mad(df$est, constant = 1)
    
    # convergence info
    conv_row <- convergence_outliers_summary[
      convergence_outliers_summary$Condition == d$condition &
        convergence_outliers_summary$Method == d$method &
        convergence_outliers_summary$Equation == d$equation &
        convergence_outliers_summary$Parameter == d$parameter, ]
    
    if (nrow(conv_row) == 0) {
      warning(paste("No convergence info for", d$condition, d$method, d$equation, d$parameter))
      next
    }
    
    results_summary <- dplyr::bind_rows(
      results_summary,
      tibble::tibble(
        Condition = d$condition,
        Method = d$method,
        Model = d$model,
        Distribution = d$distribution,
        SampleSize = d$sample_size,
        Reliability = d$reliability,
        Equation = d$equation,
        Parameter = d$parameter,
        TrueValue = tv,
        
        # Estimates
        MeanEstimate = mean(df$est),
        MedianEstimate = median_est,
        
        # Bias metrics
        Bias_Mean = abs_metrics$bias,
        Bias_Median = bias_median,
        RelativeBias_Mean = rel_bias_mean,
        RelativeBias_Median = if(tv != 0) bias_median / tv else NA,
        PercentRelativeBias_Mean = rel_bias_mean_pct,
        PercentRelativeBias_Median = if(tv != 0) bias_median / tv * 100 else NA,
        
        # Variability
        MAD = mad_est,
        Variance = abs_metrics$var,
        SD = abs_metrics$stddev,
        MSE_Mean = abs_metrics$mse,
        RMSE_Mean = abs_metrics$rmse,
        RMSE_Median = sqrt(bias_median^2 + mad_est^2),
        Relative_MSE = rel_metrics$rel_mse %||% NA,
        Relative_RMSE = rel_metrics$rel_rmse %||% NA,
        
        # SE metrics
        MeanSE = mean(df$se),
        SE_SD_Ratio = mean(df$se) / abs_metrics$stddev,
        
        # Coverage
        CoverageRate = cov_metrics$coverage * 100,
        CI_Width = cov_metrics$width,
        
        # Testing
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
        
        # Convergence and outliers
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

# Wrapper function for Study 2
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
  
  extraction_results <- ExtractConvergenceOutliers_Study2(
    all_results = all_results,
    parameters_of_interest = parameters_of_interest,
    remove_outliers = remove_outliers,
    outlier_threshold = outlier_threshold,
    min_reps = min_reps,
    exclude_warnings = exclude_warnings
  )
  
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