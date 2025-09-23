library(dplyr); library(tidyr); library(purrr); library(stringr); library(kableExtra); library(glue); library(simhelpers)

load("Simulations/Study_1/Data/Results_Study_1_final.RData")

CalculatePerformance <- function(all_results,
                                     parameters_of_interest = c("eta1","eta2","eta1:eta2","eta1:eta1","eta2:eta2"),
                                     remove_outliers = TRUE,
                                     outlier_threshold = 3, 
                                     alpha = 0.05,
                                     min_reps = 10,
                                     return_convergence_details = TRUE) {  
  options(scipen = 999)
  if (is.null(all_results) || length(all_results) == 0) {
    stop("all_results is NULL or empty. Please load your data first.")
  }
  methods <- c("sam", "lms", "qml", "dblcent")
  results_summary <- dplyr::tibble()
  param_index <- c("eta1"=1, "eta2"=2, "eta1:eta2"=3, "eta1:eta1"=4, "eta2:eta2"=5)
  
  convergence_details <- list()
  
  for (i in seq_along(all_results)) {
    condition <- all_results[[i]]$condition
    res_list  <- all_results[[i]]$results
    
    # true values (fallbacks)
    true_values <- all_results[[i]]$true_parameters
    if (is.null(true_values)) {
      if (!is.null(condition$Model_Type)) {
        if (condition$Model_Type == "null") {
          true_values <- c(0.316, 0.316, 0, 0, 0)
        } else if (condition$Model_Type == "alternative") {
          true_values <- c(0.316, 0.316, 0.139, 0.101, 0.101)
        }
      } else {
        true_values <- c(0.316, 0.316, 0.139, 0.101, 0.101)
      }
    }
    model_name <- if (!is.null(condition$Model_Type)) {
      ifelse(condition$Model_Type == "null", "Linear", "Full")
    } else "Full"
    distribution <- as.character(condition$Distribution)
    if (is.null(distribution) || distribution == "") distribution <- "normal"
    sample_size <- condition$N
    reliability <- condition$Rel

    # TRACK CONVERGENCE FOR EACH METHOD SEPARATEL
    
    # number of replications
    max_reps <- 0
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
        # this method doesn't exist, all remain false
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
            # required parameter is missing (?)
            converged <- FALSE
            break
          }
        }
        method_convergence[[method]][rep] <- converged
      }
    }
    
    # sconvergence details for this condition
    convergence_details[[paste0("Condition_", i)]] <- list(
      condition_info = list(
        Model = model_name,
        Distribution = distribution,
        SampleSize = sample_size,
        Reliability = reliability
      ),
      convergence_matrix = do.call(cbind, method_convergence),
      converged_iterations = lapply(method_convergence, which),  # which iterations converged for each method
      convergence_counts = sapply(method_convergence, sum),
      convergence_rates = sapply(method_convergence, sum) / max_reps * 100
    )
    colnames(convergence_details[[paste0("Condition_", i)]]$convergence_matrix) <- methods
    
    # individual method convergence statistics
    method_convergence_counts <- list()
    method_convergence_rates <- list()
    method_outlier_counts <- list()  
    
    for (method in methods) {
      method_convergence_counts[[method]] <- sum(method_convergence[[method]])
      method_convergence_rates[[method]] <- sum(method_convergence[[method]]) / max_reps * 100
      method_outlier_counts[[method]] <- 0 
    }
    
    # TRACK OUTLIERS FOR EACH METHOD INDEPENDENTLY (BEFORE COMPLETE CASE)
    
    if (remove_outliers) {
      for (method in methods) {
        tbls <- res_list[[paste0(method, "_tables")]]
        if (is.null(tbls)) next
        
        # all estimates for this method (only from converged cases)
        all_estimates <- list()
        for (p in parameters_of_interest) {
          all_estimates[[p]] <- numeric()
        }
        
        for (rep in which(method_convergence[[method]])) {
          extr <- extract_eta3_parameters(tbls[[rep]], method)
          if (!is.null(extr)) {
            for (p in names(extr$Estimates)) {
              if (p %in% parameters_of_interest) {
                all_estimates[[p]] <- c(all_estimates[[p]], extr$Estimates[p])
              }
            }
          }
        }
        
        # outliers across all parameters for this method
        total_outliers <- 0
        for (p in parameters_of_interest) {
          if (length(all_estimates[[p]]) > 0) {
            q1 <- stats::quantile(all_estimates[[p]], 0.25, na.rm = TRUE)
            q3 <- stats::quantile(all_estimates[[p]], 0.75, na.rm = TRUE)
            iqr <- q3 - q1
            lb <- q1 - outlier_threshold * iqr
            ub <- q3 + outlier_threshold * iqr
            outliers <- sum(all_estimates[[p]] < lb | all_estimates[[p]] > ub, na.rm = TRUE)
            total_outliers <- total_outliers + outliers
          }
        }
        method_outlier_counts[[method]] <- total_outliers
      }
    }

    # IDENTIFY VALID REPLICATIONS ACROSS ALL METHODS
    
    # replication is valid only if all methods converged
    valid_reps <- rep(TRUE, max_reps)
    for (method in methods) {
      valid_reps <- valid_reps & method_convergence[[method]]
    }
    
    # indices of valid replications
    valid_rep_indices <- which(valid_reps)
    n_total <- max_reps
    n_valid_complete_case <- length(valid_rep_indices)
    complete_case_convergence_rate <- n_valid_complete_case / n_total * 100
    
    # add complete case info to convergence details
    convergence_details[[paste0("Condition_", i)]]$complete_case_iterations <- valid_rep_indices
    convergence_details[[paste0("Condition_", i)]]$complete_case_count <- n_valid_complete_case
    convergence_details[[paste0("Condition_", i)]]$complete_case_rate <- complete_case_convergence_rate
    
    if (n_valid_complete_case < min_reps) {
      message(paste("Condition", i, ": Only", n_valid_complete_case, 
                    "valid replications across all methods. Skipping."))
      next
    }
    
    # PROCESS EACH METHOD USING ONLY VALID REPLICATIONS
    
    for (method in methods) {
      tbls <- res_list[[paste0(method, "_tables")]]
      if (is.null(tbls)) next
      
      # convergence stats for this specific method
      n_converged_this_method <- method_convergence_counts[[method]]
      convergence_rate_this_method <- method_convergence_rates[[method]]
      n_outliers_this_method <- method_outlier_counts[[method]]
      
      # accumulators by parameter - now only for valid reps
      params <- parameters_of_interest
      ests <- setNames(vector("list", length(params)), params)
      ses  <- setNames(vector("list", length(params)), params)
      pvs  <- setNames(vector("list", length(params)), params)
      los  <- setNames(vector("list", length(params)), params)
      his  <- setNames(vector("list", length(params)), params)
      
      # gather only from valid replications
      for (rep in valid_rep_indices) {
        tab <- tbls[[rep]]
        if (is.null(tab)) next 
        extr <- extract_eta3_parameters(tab, method)
        if (is.null(extr)) next 
        
        have <- intersect(names(extr$Estimates), params)
        if (length(have) == 0) next
        for (p in have) {
          ests[[p]] <- c(ests[[p]], unname(extr$Estimates[p]))
          ses [[p]] <- c(ses [[p]], unname(extr$`Standard Errors`[p]))
          pvs [[p]] <- c(pvs [[p]], unname(extr$`P-values`[p]))
          los [[p]] <- c(los [[p]], unname(extr$CI_lower[p]))
          his [[p]] <- c(his [[p]], unname(extr$CI_upper[p]))
        }
      }
      
      # compute summaries per parameter
      for (p_name in params) {
        df <- tibble::tibble(
          est = unlist(ests[[p_name]]),
          se  = unlist(ses [[p_name]]),
          p_val = unlist(pvs [[p_name]]),
          lower_bound = unlist(los [[p_name]]),
          upper_bound = unlist(his [[p_name]])
        )
        if (nrow(df) == 0) next
        
        tv <- true_values[param_index[p_name]]
        
        # all data should already be complete cases now (?)
        dfv <- df

        # CALCULATE MEDIAN-BASED METRICS BEFORE OUTLIER REMOVAL

        median_est <- stats::median(dfv$est)
        bias_median <- median_est - tv
        mad_est <- stats::mad(dfv$est, constant = 1)
        rmse_median <- sqrt(bias_median^2 + mad_est^2)
        
        if (tv != 0) {
          rel_bias_median <- bias_median / tv
          rel_bias_median_pct <- rel_bias_median * 100
        } else {
          rel_bias_median <- NA_real_
          rel_bias_median_pct <- NA_real_
        }
        
        # REMOVE OUTLIERS FOR MEAN-BASED METRICS
        
        n_out <- 0
        if (remove_outliers && nrow(dfv) > 0) {
          q1 <- stats::quantile(dfv$est, 0.25, na.rm = TRUE)
          q3 <- stats::quantile(dfv$est, 0.75, na.rm = TRUE)
          iqr <- q3 - q1
          lb <- q1 - outlier_threshold * iqr
          ub <- q3 + outlier_threshold * iqr
          keep <- dfv$est >= lb & dfv$est <= ub
          n_out <- sum(!keep, na.rm = TRUE)
          dfv_clean <- dfv[keep, , drop = FALSE]  # cleaned dataset
        } else {
          dfv_clean <- dfv  # no outlier removal
        }
        
        # check if we still have enough data after outlier removal
        if (nrow(dfv_clean) < min_reps) next
        
        # add true parameter column for simhelpers functions
        dfv_clean$true_param <- tv
        
        # CALCULATE ALL MEAN-BASED METRICS USING SIMHELPERS ON CLEANED DATA

        # ABSOLUTE METRICS (bias, variance, MSE, RMSE) without winsorization 
        abs_metrics <- calc_absolute(
          data = dfv_clean,
          estimates = est, 
          true_param = true_param,
          criteria = c("bias", "variance", "stddev", "mse", "rmse"),
          winz = Inf  # no winsorization, we already removed outliers
        )
        
        # RELATIVE METRICS (only if true value is non-zero)
        if (tv != 0) {
          rel_metrics <- calc_relative(
            data = dfv_clean,
            estimates = est,
            true_param = true_param,
            criteria = c("relative bias", "relative mse", "relative rmse"),
            winz = Inf  # No winsorization
          )
          
          # simhelpers returns E(T)/θ (ratio)
          rel_bias_mean <- rel_metrics$rel_bias - 1
          rel_bias_mean_pct <- rel_bias_mean * 100       # convert to percentage
          rel_bias_mcse <- rel_metrics$rel_bias_mcse     # MCSE for the ratio
          
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
          data = dfv_clean,
          lower_bound = lower_bound, 
          upper_bound = upper_bound, 
          true_param = true_param,
          criteria = c("coverage", "width"),
          winz = Inf  # No winsorization
        )
        
        # REJECTION RATES (Type I Error / Power)
        rej_metrics <- calc_rejection(
          data = dfv_clean,
          p_values = p_val,
          alpha = alpha
        )
        
        # convert rejection rate to percentage but keep MCSE as proportion
        rejection_rate <- rej_metrics$rej_rate * 100
        rej_mcse <- rej_metrics$rej_rate_mcse  # as proportion
        
        # determine Type I Error vs Power
        if (tv == 0) {
          type_i <- rejection_rate
          power <- NA_real_
          typei_mcse <- rej_mcse  # as proportion
          power_mcse <- NA_real_
        } else {
          type_i <- NA_real_
          power <- rejection_rate
          typei_mcse <- NA_real_
          power_mcse <- rej_mcse  # as proportion
        }

        # ADDITIONAL METRICS
        
        # mean estimate (on cleaned data)
        mean_est <- mean(dfv_clean$est)
        
        # standard error metrics (on cleaned data)
        mean_se <- mean(dfv_clean$se)
        se_sd_ratio <- mean_se / abs_metrics$stddev
        
        # COMPILE RESULTS
        
        results_summary <- dplyr::bind_rows(
          results_summary,
          tibble::tibble(
            Condition = i,
            Method = toupper(method),
            Model = model_name,
            Distribution = distribution,
            SampleSize = sample_size,
            Reliability = reliability,
            Parameter = p_name,
            TrueValue = tv,
            
            # mean-based metrics from simhelpers (on cleaned data)
            MeanEstimate = mean_est,
            Bias_Mean = abs_metrics$bias,
            RelativeBias_Mean = rel_bias_mean,
            PercentRelativeBias_Mean = rel_bias_mean_pct,
            
            # median-based metrics (calculated before outlier removal)
            MedianEstimate = median_est,
            Bias_Median = bias_median,
            RelativeBias_Median = rel_bias_median,
            PercentRelativeBias_Median = rel_bias_median_pct,
            MAD = mad_est,
            
            # variability and accuracy metrics from simhelpers (on cleaned data)
            Variance = abs_metrics$var,
            SD = abs_metrics$stddev,
            MSE_Mean = abs_metrics$mse,
            RMSE_Mean = abs_metrics$rmse,
            RMSE_Median = rmse_median,  # on pre-outlier-removal median
            
            # relative metrics from simhelpers (on cleaned data)
            Relative_MSE = rel_mse,
            Relative_RMSE = rel_rmse,
            
            # standard error metrics (on cleaned data)
            MeanSE = mean_se,
            SE_SD_Ratio = se_sd_ratio,
            
            # coverage and CI metrics from simhelpers (on cleaned data)
            CoverageRate = cov_metrics$coverage * 100,
            CI_Width = cov_metrics$width,
            
            # testing metrics (on cleaned data)
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
            
            # sample size tracking
            N_Total = n_total,
            N_Converged_This_Method = n_converged_this_method,  # how many converged for THIS method
            Convergence_Rate_This_Method = convergence_rate_this_method,  # % converged for THIS method
            N_Used_After_Convergence = n_valid_complete_case,  # what's actually used after complete-case deletion
            Complete_Case_Convergence_Rate = complete_case_convergence_rate,  # % converged for ALL methods
            N_Outliers_Removed = n_out,  # outliers removed from this method/parameter
            Percent_Outliers_Removed = round(n_out / n_valid_complete_case * 100, 2),  # % of complete cases that were outliers
            N_Final_After_Outliers = nrow(dfv_clean),  # final N after convergence AND outlier removal
            
            # TOTAL EXCLUSION SUMMARY
            N_Total_Excluded = n_total - nrow(dfv_clean),
            Percent_Total_Excluded = round((n_total - nrow(dfv_clean)) / n_total * 100, 2)
          )
        )
      } # per-parameter
    } # per-method
  } # per-condition
  
  if (return_convergence_details) {
    list(
      results = results_summary,
      convergence_details = convergence_details
    )
  } else {
    results_summary
  }
}

results_with_mcse <- CalculatePerformance(
  all_results,
  parameters_of_interest = c("eta1","eta2","eta1:eta2","eta1:eta1","eta2:eta2"),
  remove_outliers = TRUE,
  outlier_threshold = 3,
  alpha = 0.05,
  min_reps = 10,
  return_convergence_details = FALSE  
)


